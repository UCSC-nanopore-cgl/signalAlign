#!/usr/bin/env python
"""Plot ROC curve of variant called data"""
########################################################################
# File: plot_variant_accuracy.py
#  executable: plot_variant_accuracy.py
#
# Author: Andrew Bailey
# History: Created 01/07/19
########################################################################

from argparse import ArgumentParser
import pandas as pd
import pickle
import os
import sys
import platform
import matplotlib as mpl

if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
if platform.system() == "Darwin":
    mpl.use("macosx")
from py3helpers.utils import list_dir
from py3helpers.classification import ClassificationMetrics
from py3helpers.utils import load_json, create_dot_dict
from signalalign.filter_reads import find_fast5s_from_ids_readdb, write_readdb, copy_files_from_readdb
from signalalign.variantCaller import AggregateOverReadsFull
from signalalign.utils.sequenceTools import CustomAmbiguityPositions
from timeit import default_timer as timer


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--config', '-c', action='store',
                        dest='config', required=True, type=str, default=None,
                        help="Path to json config file")

    args = parser.parse_args()
    return args


def plot_roc_from_config(config):
    """Plotting function to handle logic of the config file. Mainly created to test function"""
    config = create_dot_dict(config)

    variants = config.variants
    samples = config.samples
    if isinstance(config.threshold, float):
        threshold = config.threshold
    else:
        threshold = 0.500000001

    if isinstance(config.jobs, int):
        n_processes = config.jobs
    else:
        n_processes = 2

    save_fig_dir = config.save_fig_dir

    assert len(samples) > 0, "Must include samples in order to do comparison"
    aor_handles = []
    gwa_lables_list = []
    per_site_label_list = []
    plot_per_read = False
    plot_genome_position_aggregate = False
    plot_per_call = False

    # process samples
    for sample in samples:
        tsvs = sample.full_tsvs
        positions = sample.positions_file
        label = sample.label
        aor_h = AggregateOverReadsFull(tsvs, variants, verbose=True, processes=n_processes)
        aor_h.marginalize_over_all_reads()
        aor_handles.append(aor_h)
        assert positions or label, "Must provide either a label: {} or a positions file: {}".format(label,
                                                                                                    positions)
        # use character as label if given
        if label:
            plot_genome_position_aggregate = True
            plot_per_call = True
            plot_per_read = True
            aor_h.aggregate_position_probs = aor_h.generate_labels2(predicted_data=aor_h.aggregate_position_probs,
                                                                    true_char=label)
            aor_h.per_read_data = aor_h.generate_labels2(predicted_data=aor_h.per_read_data,
                                                         true_char=label)
            aor_h.per_position_data = aor_h.generate_labels2(predicted_data=aor_h.per_position_data,
                                                             true_char=label)

        # if positions file is given, check accuracy from that
        elif positions:
            plot_genome_position_aggregate = True
            plot_per_call = True

            genome_position_labels = CustomAmbiguityPositions.parseAmbiguityFile(positions)
            aor_h.aggregate_position_probs = aor_h.generate_labels(labelled_positions=genome_position_labels,
                                                                   predicted_data=aor_h.aggregate_position_probs)
            aor_h.per_position_data = aor_h.generate_labels(labelled_positions=genome_position_labels,
                                                            predicted_data=aor_h.per_position_data)

    # plot per read ROC curve
    if plot_per_read:
        all_per_read_labels = pd.concat([x.per_read_data for x in aor_handles], ignore_index=True)
        data_type_name = "per_read"
        plot_all_roc_curves(all_per_read_labels, variants, save_fig_dir, data_type_name, threshold=threshold)

    # plot per call ROC curve
    if plot_per_call:
        all_site_labels = pd.concat([x.per_position_data for x in aor_handles], ignore_index=True)
        data_type_name = "per_site_per_read"
        plot_all_roc_curves(all_site_labels, variants, save_fig_dir, data_type_name, threshold=threshold)

    # plot genome position calls
    if plot_genome_position_aggregate:
        all_genome_positions_labels = pd.concat([x.aggregate_position_probs for x in aor_handles], ignore_index=True)
        data_type_name = "per_genomic_site"
        plot_all_roc_curves(all_genome_positions_labels, variants, save_fig_dir, data_type_name, label_key="contig",
                            threshold=threshold)

    return 0


def plot_roc_and_precision_and_save_data(per_read_labels_only, per_read_probs_only, name, variants, save_fig_dir,
                                         label_ids=None, threshold=0.5):
    roc_h = ClassificationMetrics(per_read_labels_only, per_read_probs_only, label_ids=label_ids)
    for variant in variants:
        roc_path = None
        precision_recall_path = None
        confusion_recall_path = None

        if save_fig_dir:
            roc_path = os.path.join(save_fig_dir, "{}_roc_{}".format(name, variant))
            precision_recall_path = os.path.join(save_fig_dir, "{}_pr_{}".format(name, variant))
            confusion_recall_path = os.path.join(save_fig_dir, "{}_confusion_{}".format(name, variant))

        roc_h.plot_roc(variant, title="{} ROC for {}".format(name, variant), save_fig_path=roc_path)
        roc_h.plot_precision_recall(variant, title="{} Precison Recall for {}".format(name, variant),
                                    save_fig_path=precision_recall_path)
        roc_h.plot_confusion_matrix(title="{} Confusion Matrix for {}".format(name, variant),
                                    save_fig_path=confusion_recall_path, threshold=threshold, class_n=variant)

    print("{} confusion matrix".format(name))
    print(roc_h.confusion_matrix())
    # save pickle of classification metrics class
    if save_fig_dir:
        path = os.path.join(save_fig_dir, "{}_classificationMetrics.pkl".format(name))
        with open(path, "wb") as f:
            pickle.dump(roc_h, f)
    return 0


def plot_all_roc_curves(all_labels, variants, save_fig_dir, data_type_name, label_key="read_name", threshold=0.5):
    all_per_read_labels_template = all_labels[all_labels["strand"] == 't']
    all_per_read_labels_complement = all_labels[all_labels["strand"] == 'c']

    names = ["{}_template".format(data_type_name), "{}_complement".format(data_type_name),
             "{}_total".format(data_type_name)]
    for name, data in zip(names, [all_per_read_labels_template, all_per_read_labels_complement, all_labels]):
        per_read_labels_only = data[[x + "_label" for x in variants]]
        per_read_probs_only = data[list(variants)]
        label_ids = list(data[label_key])
        per_read_labels_only.columns = list(variants)
        if len(per_read_labels_only) > 0 and len(per_read_probs_only) > 0:
            plot_roc_and_precision_and_save_data(per_read_labels_only, per_read_probs_only, name, variants,
                                                 save_fig_dir, label_ids=label_ids, threshold=threshold)

    if save_fig_dir is not None:
        all_labels.to_csv(os.path.join(save_fig_dir, data_type_name + ".tsv"), sep='\t', index=False)


def log_tp_fn_overlap(tp_classifications, fn_classifications, class_name):
    """Write out ids that overlap between true positives and false negatves"""
    tp_ids = set(tp_classifications.get_tp_ids(class_name))
    fn_ids = set(fn_classifications.get_fn_ids(class_name))
    ids = tp_ids & fn_ids
    return ids


def load_classifcation_metrics_pkl(classification_pkl):
    """Load a classificationMetrics pickle file"""
    with open(classification_pkl, 'rb') as fh:
        cm_h = pickle.load(fh)
    return cm_h


def write_tp_fn_overlap_readdb(readdb, tp_pkl, fn_pkl, class_n, out_path, read_dirs, recursive=False):
    """Write a readdb file of reads which were true positives in one experiment and false negatives
    in another experiment
    :param readdb: read db file with ids for all data in both experiments
    :param tp_pkl: path to ClassificationMetrics pkl data where  true positives are going to be gathered
    :param fn_pkl: path to ClassificationMetrics pkl data where false negatives are going to be gathered
    :param class_n: name of the class to inspect
    :param out_path: output path for readdb file
    """
    tp_metrics = load_classifcation_metrics_pkl(tp_pkl)
    fn_metrics = load_classifcation_metrics_pkl(fn_pkl)
    overlap_ids = [x.split(".")[0] for x in log_tp_fn_overlap(tp_metrics, fn_metrics, class_n)]
    data = [[id_name, f5_path] for id_name, f5_path in find_fast5s_from_ids_readdb(readdb, overlap_ids,
                                                                                   read_dirs, recursive=recursive)]
    write_readdb(data, out_path)
    print("{} tp in {} and fp in {}".format(len(overlap_ids), tp_pkl, fn_pkl))
    return len(overlap_ids)


def main(config=None):
    start = timer()
    if config is None:
        args = parse_args()
        # load model files
        assert os.path.exists(args.config), "Config file does not exist: {}".format(args.config)
        config = load_json(args.config)

    plot_roc_from_config(config)
    stop = timer()
    print("Running Time = {} seconds".format(stop - start), file=sys.stderr)


if __name__ == '__main__':
    main()
