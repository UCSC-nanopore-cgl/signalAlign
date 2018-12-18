#!/usr/bin/env python
"""Test trainModels.py"""
########################################################################
# File: test_variantcaller.py
#  executable: test_variantcaller.py
#
# Author: Andrew Bailey
# History: 12/13/18 Created
########################################################################


import sys
import os
import unittest
import tempfile
from shutil import copyfile
from py3helpers.utils import list_dir
from signalalign.variantCaller import *


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

    def test_aggregate_all_variantcalls(self):
        aor_h = AggregateOverReads(self.variant_files, "CE")
        all_data = aor_h.aggregate_all_variantcalls()
        for i, data in all_data.iterrows():
            self.assertEqual(data[0], "gi_ecoli")
            self.assertLessEqual(data[3]+data[4], 1)

    def test_normalize_all_data(self):
        aor_h = AggregateOverReads(self.variant_files, "CE")
        all_data = aor_h.aggregate_all_variantcalls()
        data_gen = aor_h._normalize_all_data(all_data)
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


if __name__ == '__main__':
    unittest.main()