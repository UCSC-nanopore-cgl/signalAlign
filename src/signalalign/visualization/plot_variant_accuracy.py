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
import os
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
    for sample in samples:
        tsvs = sample.variant_tsvs
        positions = sample.positions_file
        label = sample.label
        aor_h = AggregateOverReads(tsvs, variants)
        aor_h.marginalize_over_all_reads()
        aor_handles.append(aor_h)
        # if positions file is given, check accuracy from that
        if positions:
            plot_genome_position_aggregate = True
            plot_per_call = True

            genome_position_labels = CustomAmbiguityPositions.parseAmbiguityFile(positions)
            genome_wide_aggregate_label = aor_h.generate_labels(labelled_positions=genome_position_labels,
                                                                predicted_data=aor_h.aggregate_position_probs)
            gwa_lables_list.append(genome_wide_aggregate_label)

            per_site_label = aor_h.generate_labels(labelled_positions=genome_position_labels,
                                                   predicted_data=aor_h.per_position_data)
            per_site_label_list.append(per_site_label)

        # use character as label if given
        if label:
            plot_per_read = True
            for nuc in variants:
                if nuc == label:
                    # set
                    aor_h.per_read_data.loc[:, "{}_label".format(nuc)] = pd.Series(1, index=aor_h.per_read_data.index)
                else:
                    aor_h.per_read_data.loc[:, "{}_label".format(nuc)] = pd.Series(0, index=aor_h.per_read_data.index)

    # plot per read ROC curve
    if plot_per_read:
        all_per_read_labels = pd.concat([x.per_read_data for x in aor_handles])

        per_read_labels_only = all_per_read_labels[[x+"_label" for x in variants]]
        per_read_probs_only = all_per_read_labels[list(variants)]
        per_read_labels_only.columns = list(variants)

        roc_h = ClassificationMetrics(per_read_labels_only, per_read_probs_only)
        for variant in variants:
            path = None
            if save_fig_dir:
                path = os.path.join(save_fig_dir, "per_read_roc_{}".format(variant))

            roc_h.plot_roc(variant, title="Per read ROC for {}".format(variant), save_fig_path=path)

    # plot per call ROC curve
    if plot_per_call:
        all_site_labels = pd.concat([x for x in per_site_label_list])
        gw_labels_only = all_site_labels[[x+"_label" for x in variants]]
        gw_probs_only = all_site_labels[list(variants)]
        gw_labels_only.columns = list(variants)
        roc_h = ClassificationMetrics(gw_labels_only, gw_probs_only)
        for variant in variants:
            path = None
            if save_fig_dir:
                path = os.path.join(save_fig_dir, "per_site_per_read_roc_{}".format(variant))

            roc_h.plot_roc(variant, title="Per site per read ROC for {}".format(variant), save_fig_path=path)

    # plot genome position calls
    if plot_genome_position_aggregate:
        all_genome_positions_labels = pd.concat([x for x in gwa_lables_list])
        gw_labels_only = all_genome_positions_labels[[x+"_label" for x in variants]]
        gw_probs_only = all_genome_positions_labels[list(variants)]
        gw_labels_only.columns = list(variants)
        roc_h = ClassificationMetrics(gw_labels_only, gw_probs_only)
        for variant in variants:
            path = None
            if save_fig_dir:
                path = os.path.join(save_fig_dir, "per_site_genomic_roc_{}".format(variant))

            roc_h.plot_roc(variant, title="Per site on genome ROC for {}".format(variant), save_fig_path=path)

    return 0


def main(config=None):
    if config is None:
        args = parse_args()
        # load model files
        assert os.path.exists(args.config), "Config file does not exist: {}".format(args.config)
        config = load_json(args.config)

    plot_roc_from_config(config)


if __name__ == '__main__':
    main()
