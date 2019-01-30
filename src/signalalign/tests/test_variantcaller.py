#!/usr/bin/env python
"""Test trainModels.py"""
########################################################################
# File: test_variantcaller.py
#  executable: test_variantcaller.py
#
# Author: Andrew Bailey
# History: 12/13/18 Created
########################################################################

import os
import unittest
import tempfile
from py3helpers.utils import list_dir, load_json
from signalalign.variantCaller import *
from signalalign.visualization.plot_variant_accuracy import plot_roc_from_config

class TestVariantCaller(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestVariantCaller, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        # cls.positions_file = os.path.join(cls.HOME, "tests/test_sequences/pUC19_SspI_Zymo.fa")
        cls.f5_files = os.path.join(cls.HOME, "tests/minion_test_reads/canonical_ecoli_R9")
        cls.runSA_config = os.path.join(cls.HOME, "tests/test_variantCalled_files/runSignalAlign-config_tmp.json")
        cls.ecoli_positions = os.path.join(cls.HOME, "tests/test_position_files/ecoli_positions_CCAGG.tsv")
        cls.variant_files = os.path.join(cls.HOME, "tests/test_variantCalled_files")
        cls.plot_variants_config = os.path.join(cls.HOME, "tests/test_variantCalled_files/plot_variants_config.json")

    def test_aggregate_all_variantcalls(self):
        aor_h = AggregateOverReads(self.variant_files, "CE")
        for i, data in aor_h.per_position_data.iterrows():
            self.assertEqual(data["contig"], "gi_ecoli")
            self.assertLessEqual(data["C"]+data["E"], 1)

    def test_marginalize_over_all_reads(self):
        aor_h = AggregateOverReads(self.variant_files, "CE")
        all_data = aor_h.marginalize_over_all_reads()
        for i, data in all_data.iterrows():
            self.assertEqual(data[0], "gi_ecoli")
            self.assertLessEqual(data["C"]+data["E"], 1)

    def test_write_data(self):
        aor_h = AggregateOverReads(self.variant_files, "CE")
        all_data = aor_h.marginalize_over_all_reads()

        with tempfile.TemporaryDirectory() as tempdir:
            new_file = os.path.join(tempdir, "test.txt")
            aor_h.write_data(new_file)
            data = pd.read_csv(new_file, delimiter="\t")
            self.assertTrue(data.equals(all_data))

    def test_normalize_all_data(self):
        aor_h = AggregateOverReads(self.variant_files, "CE")
        data_gen = aor_h._normalize_all_data(aor_h.per_position_data)
        for data in data_gen:
            self.assertEqual(data[0], "gi_ecoli")
            self.assertLessEqual(data[3]+data[4], 1)

    def test_MarginalizeVariants(self):
        for test_file in list_dir(self.variant_files, ext="tsv"):
            mv_h = MarginalizeVariants(test_file, variants="CE")
            position_probs = mv_h.get_data()
            for i, data in position_probs.iterrows():
                # self.contig, pos, self.strand], nuc_data
                self.assertEqual(data["contig"], "gi_ecoli")
                self.assertLessEqual(data["E"] + data["C"], 1)
            for i, data in mv_h.per_read_calls.iterrows():
                self.assertLessEqual(data["E"] + data["C"], 1)
                self.assertEqual(data["contig"], "gi_ecoli")

    def test_create_labels_from_positions_file(self):
        labels = create_labels_from_positions_file(self.ecoli_positions, "CE")
        for i, label in labels.iterrows():
            if label["change_to"] == "C":
                self.assertTrue(label["C"] == 1)
                self.assertTrue(label["E"] == 0)
            if label["change_to"] == "E":
                self.assertTrue(label["E"] == 1)
                self.assertTrue(label["C"] == 0)

    def test_get_true_character(self):
        labels = create_labels_from_positions_file(self.ecoli_positions, "CE")
        for i, label in labels.iterrows():
            self.assertEqual(get_true_character(labels, "gi_ecoli", label["strand"], label["position"]), label["change_to"])

    def test_generate_labels(self):
        labels = create_labels_from_positions_file(self.ecoli_positions, "CE")
        aor_h = AggregateOverReads(self.variant_files, "CE")
        data = aor_h.generate_labels(labels, aor_h.per_position_data)
        for i, row in data.iterrows():
            if row["C"] > row["E"]:
                self.assertEqual(1, row["label_E"])
                self.assertEqual(0, row["label_C"])

            if row["C"] > row["E"]:
                self.assertEqual(0, row["label_E"])
                self.assertEqual(1, row["label_C"])

    # def test_plot_variants(self):
#     with tempfile.TemporaryDirectory() as tempdir:
#         config_dict = load_json(self.plot_variants_config)
#         config_dict["save_fig_dir"] = tempdir
#         retcode = plot_roc_from_config(config_dict)
#         self.assertEqual(retcode, 0)


if __name__ == '__main__':
    unittest.main()
