#!/usr/bin/env python
"""Generate kmer specific, position specific, and modification specific accuracy results"""
########################################################################
# File: plot_multiple_variant_accuracy.py
#  executable: plot_multiple_variant_accuracy.py
#
# Author: Andrew Bailey
# History: Created 04/02/20
########################################################################

# py3helpers
from py3helpers.seq_tools import ReferenceHandler
from py3helpers.classification import ClassificationMetrics
from py3helpers.utils import load_json, create_dot_dict, list_dir, merge_lists
# other libs
import pandas as pd
import numpy as np
# std libs
import os
import sys
import platform
from argparse import ArgumentParser
from timeit import default_timer as timer
# matplotlib
import matplotlib as mpl

if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
if platform.system() == "Darwin":
    mpl.use("macosx")


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--config', '-c', action='store',
                        dest='config', required=True, type=str, default=None,
                        help="Path to json config file")

    args = parser.parse_args()
    return args


def load_positions_file(path_to_positions_file):
    assert os.path.exists(path_to_positions_file), "Path to positions file does not exist: {}".format(
        path_to_positions_file)
    return pd.read_csv(path_to_positions_file, sep="\t", names=["contig", "reference_index", "strand", "base", "label"])


def load_sa2bed_variant_data(path_to_variant_data):
    assert os.path.exists(path_to_variant_data), "Path to variant data does not exist: {}".format(path_to_variant_data)
    return pd.read_csv(path_to_variant_data)


def get_prob_and_label(variants):
    # create probability columns
    variant_strings = list(set(variants["variants"]))
    assert len(variant_strings) > 0, "No variant data passed into function get_prob_and_label"
    n_variants = len(variant_strings[0])
    lengths = [len(x) for x in variant_strings]
    if len(variant_strings) > 1:

        assert np.sum(lengths) == lengths[0] * len(lengths), "All modifications must have the " \
                                                                                  "same number of possible variants"
        prob = variants[["prob"+str(x+1) for x in range(n_variants)]].rename(
            columns={"prob" + str(i + 1): value for i, value in enumerate([str(i) for i in range(n_variants)])})
        label = variants[["prob"+str(x+1)+"_label" for x in range(n_variants)]].rename(
            columns={"prob" + str(i + 1)+"_label": value for i, value in enumerate([str(i) for i in range(n_variants)])})

    else:
        prob = variants[["prob"+str(x+1) for x in range(n_variants)]].rename(
            columns={"prob" + str(i + 1): value for i, value in enumerate(variants["variants"].iloc[0])})
        label = variants[["prob"+str(x+1)+"_label" for x in range(n_variants)]].rename(
            columns={"prob" + str(i + 1)+"_label": value for i, value in enumerate(variants["variants"].iloc[0])})

    label_ids = ["_".join([str(y) for y in x]) for x in zip(variants["read_id"],
                                                            variants["contig"],
                                                            variants["reference_index"],
                                                            variants["strand"])]
    return label, prob, label_ids


def create_master_table(positions, variants):
    n_variants = len(variants.columns)-5
    label_data = variants[["prob"+str(x+1) for x in range(n_variants)]] * 0
    label_data.columns = ["prob"+str(x+1)+"_label" for x in range(n_variants)]
    variants2 = pd.concat([variants, label_data], axis=1)
    labelled_dfs = []
    for x, y in variants2.groupby(['contig', 'reference_index', "strand", "variants"], as_index=False):
        label_row = positions[(positions['contig'] == x[0])
                              & (positions['reference_index'] == x[1])
                              & (positions['strand'] == x[2])]
        if len(label_row) == 0:
            continue

        first_row = y.iloc[0]
        label = label_row.iloc[0]["label"]
        index = first_row["variants"].find(label)
        assert index != -1, "Variant label is not in variants at this position. Check model file and positions file"
        y.loc[:, "prob"+str(index+1)+"_label"] = 1
        labelled_dfs.append(y)
    complete_table = pd.concat(labelled_dfs)
    return complete_table


def plot_variant_data(labels, probs, label_ids, output_dir, name):
    if len(labels) != 0:
        n_variants = len(probs.columns)
        class_n = probs.columns[-1]
        cm = ClassificationMetrics(labels,
                                   probs,
                                   label_ids=label_ids)
        plot = cm.plot_confusion_matrix(threshold=0.5,
                                        title="Confusion Matrix " + name,
                                        save_fig_path=os.path.join(output_dir, name + "_confusion_matrix.png"),
                                        class_n=None)
        plot.close()
        if n_variants > 2:
            plot = cm.plot_multiclass_precision_recall(save_fig_path=os.path.join(output_dir,
                                                                                  name + "_precision_recall.png"),
                                                       title='Multi-class precision-recall ' + name)
            plot.close()

            plot = cm.plot_multiclass_roc(save_fig_path=os.path.join(output_dir,
                                                                     name + "_roc.png"),
                                          title='Multi-class ROC ' + name)
            plot.close()
            for class_n in probs.columns:
                plot = cm.plot_probability_hist(class_n,
                                                save_fig_path=os.path.join(output_dir, name + "_" +
                                                                           class_n + "_probability_hist.png"),
                                                bins=None,
                                                normalize=True)
                plot.close()
        else:
            plot = cm.plot_probability_hist(class_n,
                                            save_fig_path=os.path.join(output_dir, name + "_probability_hist.png"),
                                            bins=None,
                                            normalize=True)
            plot.close()

            plot = cm.plot_roc(class_n, save_fig_path=os.path.join(output_dir, name + "_roc.png"),
                               title="Receiver operating characteristic " + name)
            plot.close()

            plot = cm.plot_precision_recall(class_n, save_fig_path=os.path.join(output_dir,
                                                                                name + "_precision_recall.png"),
                                            title="Precision Recall curve" + name)
            plot.close()


def main(config=None):
    start = timer()
    if config is None:
        args = parse_args()
        # load model files
        assert os.path.exists(args.config), "Config file does not exist: {}".format(args.config)
        config = load_json(args.config)

    config = create_dot_dict(config)
    assert os.path.isdir(config.output_dir), "Output directory does not exist. {}".format(config.output_dir)
    print("config.output_dir: {}".format(config.output_dir))
    samples = config.samples
    # aggregate and label all data
    all_data = []
    for sample in samples:
        variants = load_sa2bed_variant_data(sample.variant_csv)
        positions = load_positions_file(sample.positions_file)
        assert len(positions) == len("".join(positions["label"])), "Positions file must have only one value for label"
        sample_data = create_master_table(positions, variants)
        if sample_data is not None:
            all_data.append(sample_data)

    assert len(all_data) > 0, \
        "No variants overlap with positions within all the samples. Check positions file and variants file"
    all_data_df = pd.concat(all_data)

    # all variant accuracy for 2, 3 ... variants
    output_dir = os.path.join(config.output_dir, "all_variants")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    possible_number_of_variants = list(set(all_data_df['variants'].str.len()))
    for x in possible_number_of_variants:
        labels, probs, label_ids = get_prob_and_label(all_data_df[all_data_df['variants'].str.len() == x])
        plot_variant_data(labels, probs, label_ids, output_dir, "variants_length_{}".format(x))

    # all variant accuracy for each variant type
    output_dir = os.path.join(config.output_dir, "per_variant")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    possible_variants = list(set(all_data_df['variants']))
    for x in possible_variants:
        labels, probs, label_ids = get_prob_and_label(all_data_df[all_data_df['variants'] == x])
        plot_variant_data(labels, probs, label_ids, output_dir, "variants_{}".format(x))

    # # all variant accuracy for each position
    output_dir = os.path.join(config.output_dir, "per_position")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for x, y in all_data_df.groupby(['contig', 'reference_index', "strand", "variants"], as_index=False):
        labels, probs, label_ids = get_prob_and_label(y)
        plot_variant_data(labels, probs, label_ids, output_dir, "_".join([str(i) for i in x]))

    stop = timer()
    print("Running Time = {} seconds".format(stop - start), file=sys.stderr)


if __name__ == '__main__':
    main()
