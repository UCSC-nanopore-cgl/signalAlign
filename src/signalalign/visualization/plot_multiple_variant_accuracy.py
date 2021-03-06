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
import shutil
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
    pd.set_option('mode.chained_assignment', None)
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


def plot_variant_data(labels, probs, label_ids, output_dir, name, threshold=0.5):
    if len(labels) != 0:
        n_variants = len(probs.columns)
        class_n = probs.columns[-1]
        cm = ClassificationMetrics(labels,
                                   probs,
                                   label_ids=label_ids)
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
            plot = cm.plot_calibration_curve(class_n,
                                             save_fig_path=os.path.join(output_dir, name + "_calibration_curve.png"),
                                             title="Calibration curve " + name,
                                             n_bins=20)
            plot.close()
            plot = cm.plot_confusion_matrix(threshold=threshold,
                                            title="Confusion Matrix " + name + " threshold:"+str(threshold),
                                            save_fig_path=os.path.join(output_dir, name + "_" +
                                                                       str(threshold) + "_confusion_matrix.png"),
                                            class_n=class_n)
            plot.close()

            plot = cm.plot_roc(class_n, save_fig_path=os.path.join(output_dir, name + "_roc.png"),
                               title="Receiver operating characteristic " + name)
            plot.close()

            plot = cm.plot_precision_recall(class_n, save_fig_path=os.path.join(output_dir,
                                                                                name + "_precision_recall.png"),
                                            title="Precision Recall curve " + name)
            plot.close()
        return cm
    return None


def plot_multiple_variant_accuracy(output_dir, samples, threshold=0.5):
    # aggregate and label all data
    all_data = []
    for sample in samples:
        sample = create_dot_dict(sample)
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
    all_variants_output_dir = os.path.join(output_dir, "all_variants")
    if not os.path.exists(all_variants_output_dir):
        os.makedirs(all_variants_output_dir)
    possible_number_of_variants = list(set(all_data_df['variants'].str.len()))
    for x in possible_number_of_variants:
        print("variants_length_{}".format(x))
        labels, probs, label_ids = get_prob_and_label(all_data_df[all_data_df['variants'].str.len() == x])
        plot_variant_data(labels, probs, label_ids, all_variants_output_dir, "variants_length_{}".format(x), threshold=threshold)

    # all variant accuracy for each variant type
    per_variant_output_dir = os.path.join(output_dir, "per_variant")
    if not os.path.exists(per_variant_output_dir):
        os.makedirs(per_variant_output_dir)
    possible_variants = list(set(all_data_df['variants']))
    for x in possible_variants:
        print("variants_{}".format(x))
        labels, probs, label_ids = get_prob_and_label(all_data_df[all_data_df['variants'] == x])
        plot_variant_data(labels, probs, label_ids, per_variant_output_dir, "variants_{}".format(x), threshold=threshold)

    # # all variant accuracy for each position
    per_position_output_dir = os.path.join(output_dir, "per_position")
    if not os.path.exists(per_position_output_dir):
        os.makedirs(per_position_output_dir)
    with open(os.path.join(per_position_output_dir, "per_position_data_"+str(threshold)+".csv"), 'w') as fh:
        print(",".join(['contig', 'reference_index', "strand", "variants",
                        "accuracy",
                        "precision",
                        "negative_predictive_value",
                        "recall",
                        "specificity",
                        "positive_likelihood_ratio",
                        "negative_likelihood_ratio",
                        "diagnostic_odds_ratio",
                        "f1_score",
                        "prevalence",
                        "aucroc", "avg_precision", "brier_score"]),
              file=fh)
        for x, y in all_data_df.groupby(['contig', 'reference_index', "strand", "variants"], as_index=False):
            print("_".join([str(i) for i in x]))
            labels, probs, label_ids = get_prob_and_label(y)
            cm = plot_variant_data(labels, probs, label_ids, per_position_output_dir, "_".join([str(i) for i in x]),
                                   threshold=threshold)
            if cm is not None:
                class_n = probs.columns[-1]
                rocauc = cm.roc_auc[class_n]
                avg_precision = cm.get_average_precision(class_n)
                brier_score = cm.brier_score[class_n]
                accuracy = cm.accuracy(class_n, threshold=threshold)
                precision = cm.precision(class_n, threshold=threshold)
                negative_predictive_value = cm.negative_predictive_value(class_n, threshold=threshold)
                recall = cm.recall(class_n, threshold=threshold)
                specificity = cm.specificity(class_n, threshold=threshold)
                positive_likelihood_ratio = cm.positive_likelihood_ratio(class_n, threshold=threshold)
                negative_likelihood_ratio = cm.negative_likelihood_ratio(class_n, threshold=threshold)
                diagnostic_odds_ratio = cm.diagnostic_odds_ratio(class_n, threshold=threshold)
                f1_score = cm.f1_score(class_n, threshold=threshold)
                prevalence = cm.prevalence(class_n)
                line = [str(i) for i in x] + [
                    str(round(accuracy, 4)),
                    str(round(precision, 4)),
                    str(round(negative_predictive_value, 4)),
                    str(round(recall, 4)),
                    str(round(specificity, 4)),
                    str(round(positive_likelihood_ratio, 4)),
                    str(round(negative_likelihood_ratio, 4)),
                    str(round(diagnostic_odds_ratio, 4)),
                    str(round(f1_score, 4)),
                    str(round(prevalence, 4)),
                    str(round(rocauc, 4)),
                    str(round(avg_precision, 4)),
                    str(round(brier_score, 4))]
                print(",".join(line), file=fh)


def main():
    start = timer()
    args = parse_args()
    # load model files
    assert os.path.exists(args.config), "Config file does not exist: {}".format(args.config)
    config = create_dot_dict(load_json(args.config))

    assert os.path.isdir(config.output_dir), "Output directory does not exist. {}".format(config.output_dir)
    print("output_dir: {}".format(config.output_dir))
    if args.config != os.path.join(config.output_dir, os.path.basename(args.config)):
        shutil.copyfile(args.config, os.path.join(config.output_dir, os.path.basename(args.config)))

    output_dir = config.output_dir
    samples = config.samples
    threshold = 0.5
    if config.threshold:
        threshold = config.threshold
    plot_multiple_variant_accuracy(output_dir, samples, threshold=threshold)

    stop = timer()
    print("Running Time = {} seconds".format(stop - start), file=sys.stderr)


if __name__ == '__main__':
    main()
