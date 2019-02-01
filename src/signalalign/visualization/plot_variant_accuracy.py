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
import matplotlib as mpl

if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
from py3helpers.utils import list_dir
from py3helpers.classification import ClassificationMetrics
from py3helpers.utils import load_json, create_dot_dict
from signalalign.variantCaller import AggregateOverReads
from signalalign.utils.sequenceTools import CustomAmbiguityPositions


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

    save_fig_dir = config.save_fig_dir

    assert len(samples) > 0, "Must include samples in order to do comparison"
    aor_handles = []
    gwa_lables_list = []
    per_site_label_list = []
    plot_per_read = False
    plot_genome_position_aggregate = False
    plot_per_call = False

    # process samples
    for strand in ("t", "c"):
        for sample in samples:
            tsvs = sample.variant_tsvs
            positions = sample.positions_file
            label = sample.label
            aor_h = AggregateOverReads(tsvs, variants)
            aor_h.marginalize_over_all_reads()
            aor_handles.append(aor_h)
            assert positions or label, "Must provide either a label: {} or a positions file: {}".format(label,
                                                                                                        positions)
            # use character as label if given
            if label:
                plot_genome_position_aggregate = True
                plot_per_call = True
                plot_per_read = True
                for nuc in variants:
                    if nuc == label:
                        # set
                        aor_h.per_read_data.loc[:, "{}_label".format(nuc)] = pd.Series(1,
                                                                                       index=aor_h.per_read_data.index)
                    else:
                        aor_h.per_read_data.loc[:, "{}_label".format(nuc)] = pd.Series(0,
                                                                                       index=aor_h.per_read_data.index)
                genome_wide_aggregate_label = aor_h.generate_labels2(predicted_data=aor_h.aggregate_position_probs,
                                                                     true_char=label)
                gwa_lables_list.append(genome_wide_aggregate_label)

                per_site_label = aor_h.generate_labels2(predicted_data=aor_h.per_position_data, true_char=label)
                per_site_label_list.append(per_site_label)

            # if positions file is given, check accuracy from that
            elif positions:
                plot_genome_position_aggregate = True
                plot_per_call = True

                genome_position_labels = CustomAmbiguityPositions.parseAmbiguityFile(positions)
                genome_wide_aggregate_label = aor_h.generate_labels(labelled_positions=genome_position_labels,
                                                                    predicted_data=aor_h.aggregate_position_probs)
                gwa_lables_list.append(genome_wide_aggregate_label)

                per_site_label = aor_h.generate_labels(labelled_positions=genome_position_labels,
                                                       predicted_data=aor_h.per_position_data)
                per_site_label_list.append(per_site_label)

        # plot per read ROC curve
        if plot_per_read:
            all_per_read_labels = pd.concat([x.per_read_data for x in aor_handles])
            data_type_name = "per_read"
            plot_all_roc_curves(all_per_read_labels, variants, save_fig_dir, data_type_name)

        # plot per call ROC curve
        if plot_per_call:
            all_site_labels = pd.concat([x for x in per_site_label_list])
            data_type_name = "per_site_per_read"
            plot_all_roc_curves(all_site_labels, variants, save_fig_dir, data_type_name)

        # plot genome position calls
        if plot_genome_position_aggregate:
            all_genome_positions_labels = pd.concat([x for x in gwa_lables_list])
            data_type_name = "per_genomic_site"
            plot_all_roc_curves(all_genome_positions_labels, variants, save_fig_dir, data_type_name)

        return 0


def plot_roc_and_save_data(per_read_labels_only, per_read_probs_only, name, variants, save_fig_dir):
    roc_h = ClassificationMetrics(per_read_labels_only, per_read_probs_only)
    for variant in variants:
        path = None
        if save_fig_dir:
            path = os.path.join(save_fig_dir, "{}_roc_{}".format(name, variant))

        roc_h.plot_roc(variant, title="{} ROC for {}".format(name, variant), save_fig_path=path)
    print("{} confusion matrix".format(name))
    print(roc_h.confusion_matrix())
    # save pickle of classification metrics class
    path = os.path.join(save_fig_dir, "{}_classificationMetrics.pkl".format(name))
    with open(path, "wb") as f:
        pickle.dump(roc_h, f)


def plot_all_roc_curves(all_labels, variants, save_fig_dir, data_type_name):
    all_per_read_labels_template = all_labels[all_labels["strand"] == 't']
    all_per_read_labels_complement = all_labels[all_labels["strand"] == 'c']

    names = ["{}_template".format(data_type_name), "{}_complement".format(data_type_name), "{}_total".format(data_type_name)]
    for name, data in zip(names, [all_per_read_labels_template, all_per_read_labels_complement, all_labels]):

        per_read_labels_only = data[[x + "_label" for x in variants]]
        per_read_probs_only = data[list(variants)]

        per_read_labels_only.columns = list(variants)

        plot_roc_and_save_data(per_read_labels_only, per_read_probs_only, name, variants, save_fig_dir)


def main(config=None):
    if config is None:
        args = parse_args()
        # load model files
        assert os.path.exists(args.config), "Config file does not exist: {}".format(args.config)
        config = load_json(args.config)

    plot_roc_from_config(config)


if __name__ == '__main__':
    main()
