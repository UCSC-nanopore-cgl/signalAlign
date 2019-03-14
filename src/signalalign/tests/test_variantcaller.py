#!/usr/bin/env python
"""Test variantcaller.py"""
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
from signalalign.train.trainModels import read_in_alignment_file


class TestVariantCaller(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestVariantCaller, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        # cls.positions_file = os.path.join(cls.HOME, "tests/test_sequences/pUC19_SspI_Zymo.fa")
        cls.f5_files = os.path.join(cls.HOME, "tests/minion_test_reads/canonical_ecoli_R9")
        cls.runSA_config = os.path.join(cls.HOME, "tests/test_variantCalled_files/runSignalAlign-config_tmp.json")
        cls.ecoli_positions = os.path.join(cls.HOME, "tests/test_position_files/CCWGG_ecoli_k12_mg1655_CC.positions")
        cls.variant_files = os.path.join(cls.HOME, "tests/test_variantCalled_files/canonical")
        cls.rna_variant_files = os.path.join(cls.HOME, "tests/test_variantCalled_files/rna")
        cls.plot_variants_config = os.path.join(cls.HOME, "tests/test_variantCalled_files/plot_variants_config.json")
        cls.rna_positions_file = os.path.join(cls.HOME, "tests/test_position_files/rna_atg_ftg_fake_ref.positions")

    def test_aggregate_all_variantcalls(self):
        aor_h = AggregateOverReadsFull(self.variant_files, "CE")
        for i, data in aor_h.per_position_data.iterrows():
            self.assertEqual(data["contig"], "gi_ecoli")
            self.assertAlmostEqual(data["C"] + data["E"], 1)

    def test_marginalize_over_all_reads(self):
        aor_h = AggregateOverReadsFull(self.variant_files, "CE")
        all_data = aor_h.marginalize_over_all_reads()
        for i, data in all_data.iterrows():
            self.assertEqual(data[0], "gi_ecoli")
            self.assertLessEqual(data["C"] + data["E"], 1)

    def test_write_data(self):
        aor_h = AggregateOverReadsFull(self.variant_files, "CE")
        all_data = aor_h.marginalize_over_all_reads()

        with tempfile.TemporaryDirectory() as tempdir:
            new_file = os.path.join(tempdir, "test.txt")
            aor_h.write_data(new_file)
            data = pd.read_csv(new_file, delimiter="\t")
            self.assertTrue(data["position"].equals(all_data["position"]))

    def test_normalize_all_data(self):
        aor_h = AggregateOverReadsFull(self.variant_files, "CE")
        data_gen = aor_h._normalize_all_data(aor_h.per_position_data)
        for data in data_gen:
            self.assertEqual(data[0], "gi_ecoli")
            self.assertLessEqual(data[4] + data[5], 1)

    def test_MarginalizeFullVariants(self):
        forward_mapped_files = list_dir(self.variant_files, ext="forward.tsv")
        for test_file in forward_mapped_files:
            read_name = os.path.basename(test_file)
            full_data = read_in_alignment_file(test_file)
            mv_h = MarginalizeFullVariants(full_data, variants="CE", read_name=read_name, forward_mapped=True)
            position_probs = mv_h.get_data()
            for i, data in position_probs.iterrows():
                self.assertEqual(data["contig"], "gi_ecoli")
                self.assertAlmostEqual(data["E"] + data["C"], 1)
            for i, data in mv_h.per_read_calls.iterrows():
                self.assertAlmostEqual(data["E"] + data["C"], 1)
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
            if i > 100:
                break
            self.assertEqual(get_true_character(labels, "gi_ecoli", label["strand"], label["position"]),
                             label["change_to"])

    def test_generate_labels(self):
        labels = create_labels_from_positions_file(self.ecoli_positions, "CE")
        aor_h = AggregateOverReadsFull(self.variant_files, "CE")
        data = aor_h.generate_labels(labels, aor_h.per_position_data)
        for i, row in data.iterrows():
            if row["read_name"] == '6e520d79-dcc5-4af3-a69f-cf5134a4c563.sm.tsv':
                self.assertEqual(1, row["E_label"])
                self.assertEqual(0, row["C_label"])
            else:
                self.assertEqual(0, row["E_label"])
                self.assertEqual(1, row["C_label"])

    def test_aggregate_all_variantcalls_full(self):
        aor_h = AggregateOverReadsFull(self.variant_files, "CE")
        for i, data in aor_h.per_position_data.iterrows():
            self.assertEqual(data["contig"], "gi_ecoli")
            self.assertAlmostEqual(data["C"] + data["E"], 1)

    def test_marginalize_over_all_reads_full(self):
        aor_h = AggregateOverReadsFull(self.variant_files, "CE")
        all_data = aor_h.marginalize_over_all_reads()
        for i, data in all_data.iterrows():
            self.assertEqual(data[0], "gi_ecoli")
            self.assertLessEqual(data["C"] + data["E"], 1)

    def test_write_data_full(self):
        aor_h = AggregateOverReadsFull(self.variant_files, "CE")
        all_data = aor_h.marginalize_over_all_reads()

        with tempfile.TemporaryDirectory() as tempdir:
            new_file = os.path.join(tempdir, "test.txt")
            aor_h.write_data(new_file)
            data = pd.read_csv(new_file, delimiter="\t")
            self.assertTrue(data["position"].equals(all_data["position"]))

    def test_plot_variants(self):
        with tempfile.TemporaryDirectory() as tempdir:
            config_dict = load_json(self.plot_variants_config)
            config_dict["save_fig_dir"] = tempdir
            retcode = plot_roc_from_config(config_dict)
            self.assertEqual(retcode, 0)

    def test_aggregate_variant_call_rna(self):
        aor_h = AggregateOverReadsFull(self.rna_variant_files, "AF")
        for i, data in aor_h.per_position_data.iterrows():
            self.assertEqual(data["contig"], "rna_fake")
            self.assertAlmostEqual(data["A"] + data["F"], 1)


if __name__ == '__main__':
    unittest.main()
