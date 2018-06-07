#!/usr/bin/env python
"""Test BuildModels.py"""
########################################################################
# File: test_BuildModels.py
#  executable: test_BuildModels.py
#
# Author: Andrew Bailey
# History: 5/21/18 Created
########################################################################


import sys
import os
import numpy as np
import unittest
from collections import defaultdict
from scipy import sparse
from signalalign.train.BuildModels import make_positions_file, make_gatc_position_file, find_gatc_motifs


class TrainSignalAlignTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TrainSignalAlignTest, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.reference = os.path.join(cls.HOME, "tests/test_sequences/pUC19_SspI_Zymo.fa")
        # cls.reference = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/tests/test_sequences/pUC19_SspI_Zymo.fa"

    # def test_make_positions_file(self):
    #     with tempfile.TemporaryDirectory() as tempdir:
    #         path = os.path.join(tempdir, "test.positions")
    #
    #         os.remove(path)
    #
    #     degenerate = "adenosine"
    #     degenerate = "cytosine"
    #     positions_file = make_positions_file(fasta=self.reference,
    #                                          degenerate=degenerate,
    #                                          outfile=path)

    def test_make_gatc_position_file(self):
        # with tempfile.TemporaryDirectory() as tempdir:
            # path = os.path.join(tempdir, "test.positions")
        path = os.path.join(self.HOME, "src/signalalign/tests/test.adenosine.positions")

        positions_file = make_gatc_position_file(fasta=self.reference,
                                                 outfile=path)
        os.remove(path)

    # def test_find_gatc_motifs(self):
    #     find_gatc_motifs(seq, )


if __name__ == '__main__':
    unittest.main()
