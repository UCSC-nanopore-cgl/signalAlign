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
from signalalign.visualization.compare_trained_models import main


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
            self.assertEqual(len(plotted_files), 2)


if __name__ == '__main__':
    unittest.main()
