#!/usr/bin/env python
"""Plot ROC curve of variant called data"""
########################################################################
# File: plot_variant_accuracy.py
#  executable: plot_variant_accuracy.py
#
# Author: Andrew Bailey
# History: Created 01/07/19
########################################################################

import pandas as pd
from py3helpers.utils import list_dir
from py3helpers.classification import ClassificationMetrics
from signalalign.variantCaller import AggregateOverReads
from signalalign.utils.sequenceTools import CustomAmbiguityPositions


def main():
    # variant_caller_tsvs = "/Users/andrewbailey/data/m6a_rna_data/tempFiles_alignment/m6a"
    # all_calls = list_dir(variant_caller_tsvs, ext="tsv")
    # variants = "AF"
    # labeled_positions_file = "/Users/andrewbailey/data/m6a_rna_data/positions_file.tsv"
    # data = CustomAmbiguityPositions.parseAmbiguityFile(labeled_positions_file)
    # aor_h = AggregateOverReads(variant_caller_tsvs, variants)
    # data2 = aor_h.marginalize_over_all_reads()
    # data3 = aor_h.aggregate_all_variantcalls()
    #
    # data4 = aor_h.generate_labels(data, data2)
    # data5 = aor_h.generate_labels(data, data3)
    #
    # labels1 = data4[[x+"_label" for x in variants]]
    # probs1 = data4[list(variants)]
    #
    # labels2 = data5[[x+"_label" for x in variants]]
    # probs2 = data5[list(variants)]
    # labels1.columns = list(variants)
    # labels2.columns = list(variants)
    #
    # # roc_h = PlotROC(labels1, probs1)
    # # roc_h.plot_multiclass_roc()
    # roc_h = ClassificationMetrics(labels2, probs2)
    # # roc_h.plot_multiclass_roc()
    # roc_h.plot_roc("A")
    # # roc_h.plot_roc("F")
    # # print(roc_h.confusion_matrix())
    # # roc_h.plot_calibration_curve("A")

    canonical_tsvs = "/Users/andrewbailey/data/ccwgg_test_variant_calls/canonical"
    mC_tsvs = "/Users/andrewbailey/data/ccwgg_test_variant_calls/5-mC"

    variants = "CE"
    canonical_positions_file = "/Users/andrewbailey/data/references/ecoli/CCWGG_ecoli_k12_mg1655_C_C.positions"
    mC_positions_file = "/Users/andrewbailey/data/references/ecoli/CCWGG_ecoli_k12_mg1655.positions"

    canonical_labels = CustomAmbiguityPositions.parseAmbiguityFile(canonical_positions_file)
    mC_labels = CustomAmbiguityPositions.parseAmbiguityFile(mC_positions_file)

    canonical_aor_h = AggregateOverReads(canonical_tsvs, variants)
    mC_aor_h = AggregateOverReads(mC_tsvs, variants)

    #################
    # canonical_genome_wide_predictions = canonical_aor_h.marginalize_over_all_reads()
    # canonical_gw_labels = canonical_aor_h.generate_labels(labelled_positions=canonical_labels, predicted_data=canonical_genome_wide_predictions)
    #
    # mC_genome_wide_predictions = mC_aor_h.marginalize_over_all_reads()
    # mC_gw_labels = canonical_aor_h.generate_labels(labelled_positions=mC_labels, predicted_data=mC_genome_wide_predictions)
    #
    #
    # all_gw_labels = pd.concat([canonical_gw_labels, mC_gw_labels])
    #
    #
    # gw_labels_only = all_gw_labels[[x+"_label" for x in variants]]
    # gw_probs_only = all_gw_labels[list(variants)]
    # gw_labels_only.columns = list(variants)
    # roc_h = ClassificationMetrics(gw_labels_only, gw_probs_only)
    # roc_h.plot_roc("C")

    #################
    # canonical_per_read_predictions = canonical_aor_h.per_position_data()
    # mC_per_read_predictions = mC_aor_h.per_position_data
    #
    # canonical_per_read_labels = canonical_aor_h.generate_labels(labelled_positions=canonical_labels, predicted_data=canonical_per_read_predictions)
    # mC_per_read_labels = canonical_aor_h.generate_labels(labelled_positions=mC_labels, predicted_data=mC_per_read_predictions)
    #
    # all_per_read_labels = pd.concat([canonical_per_read_labels, mC_per_read_labels])
    #
    # per_read_labels_only = all_per_read_labels[[x+"_label" for x in variants]]
    # per_read_probs_only = all_per_read_labels[list(variants)]
    # per_read_labels_only.columns = list(variants)
    #
    # roc_h = ClassificationMetrics(per_read_labels_only, per_read_probs_only)
    # roc_h.plot_roc("C", title="Per read DCM methylase ROC (300 reads)")


    #################
    canonical_per_read_predictions = canonical_aor_h.per_read_data
    mC_per_read_predictions = mC_aor_h.per_read_data

    # set
    mC_per_read_predictions.loc[:, "C_label"] = pd.Series(0, index=mC_per_read_predictions.index)
    canonical_per_read_predictions.loc[:, "C_label"] = pd.Series(1, index=canonical_per_read_predictions.index)
    mC_per_read_predictions.loc[:, "E_label"] = pd.Series(1, index=mC_per_read_predictions.index)
    canonical_per_read_predictions.loc[:, "E_label"] = pd.Series(0, index=canonical_per_read_predictions.index)

    all_per_read_labels = pd.concat([canonical_per_read_predictions, mC_per_read_predictions])

    per_read_labels_only = all_per_read_labels[[x+"_label" for x in variants]]
    per_read_probs_only = all_per_read_labels[list(variants)]
    per_read_labels_only.columns = list(variants)

    roc_h = ClassificationMetrics(per_read_labels_only, per_read_probs_only)
    roc_h.plot_roc("C", title="Per read DCM methylase ROC (300 reads)")


if __name__ == '__main__':
    main()
