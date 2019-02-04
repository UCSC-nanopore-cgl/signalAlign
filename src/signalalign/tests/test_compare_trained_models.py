#!/usr/bin/env python
"""Test compare_trained_models.py"""
########################################################################
# File: test_compare_trained_models.py
#  executable: test_compare_trained_models.py
#
# Author: Andrew Bailey
# History: 01/24/18 Created
########################################################################

import os
import unittest
import tempfile
from py3helpers.utils import list_dir, load_json
from signalalign.visualization.compare_trained_models import MultipleModelHandler, main


class TestCompareTrainedModels(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestCompareTrainedModels, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.test_compare_train_config = load_json(os.path.join(cls.HOME, "tests/test_compare_trained_models/compare_trained_models.config.json"))

    def test_compare_trained_models(self):
        with tempfile.TemporaryDirectory() as tempdir:
            self.test_compare_train_config["save_fig_dir"] = tempdir
            main(self.test_compare_train_config)
            plotted_files = list_dir(tempdir)
            self.assertEqual(len(plotted_files), 5)

    def test_read_and_write_kmer_distribution_comparison_logfile(self):
        with tempfile.TemporaryDirectory() as tempdir:
            outfile = os.path.join(tempdir, "test.tsv")
            d1 = list(range(10))
            d2 = list(range(9))
            d3 = list(range(11))

            d1.append(None)
            d2.extend([None, None])

            kmers = list("ISATEST!!!r")
            MultipleModelHandler.write_kmer_distribution_comparison_logfile(kmers, d1, d2, d3, outfile)
            data = MultipleModelHandler.read_kmer_distribution_comparison_logfile(outfile)

            correct_data = [[k, d1, d2, d3] for k, d1, d2, d3 in zip(kmers, d1, d2, d3) if d1 is not None]
            correct_data.sort(key=lambda x: x[1], reverse=True)
            correct_data.append(["r", None, None, 10])
            self.assertEqual(data, correct_data)


if __name__ == '__main__':
    unittest.main()
