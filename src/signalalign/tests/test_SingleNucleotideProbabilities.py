#!/usr/bin/env python3

import sys
import unittest
import glob
import os
import shutil
import pandas as pd
import numpy as np
from subprocess import call
import signalalign.singleNucleotideProbabilities as singleNuclProb

SIGNALALIGN_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
SIGNALMACHINE_EXE = os.path.join(SIGNALALIGN_ROOT, "bin/signalMachine")
TEMPLATE_MODEL=os.path.join(SIGNALALIGN_ROOT, "models/testModelR9_5mer_acgt_template.model")
# ZYMO_C_READS = os.path.join(SIGNALALIGN_ROOT, "tests/minion_test_reads/C/")
# ZYMO_REFERENCE = os.path.join(SIGNALALIGN_ROOT, "tests/test_sequences/zymo_sequence.fasta")




class SingleNuclProbsTest(unittest.TestCase):
    WORK_DIR=os.path.abspath("./singleNuclProb_unittest/")

    def setUp(self):
        if os.path.exists(SingleNuclProbsTest.WORK_DIR):
            shutil.rmtree(SingleNuclProbsTest.WORK_DIR)
        os.makedirs(SingleNuclProbsTest.WORK_DIR)
        shutil.copy(SIGNALMACHINE_EXE, os.path.join(SingleNuclProbsTest.WORK_DIR, "signalMachine"))
        os.chdir(SingleNuclProbsTest.WORK_DIR)

    def tearDown(self):
        shutil.rmtree(SingleNuclProbsTest.WORK_DIR)

    def run_single_nucl_prob(self, fast5_dir, reference_location):
        output_tmp_dir = os.path.join(os.path.abspath(SingleNuclProbsTest.WORK_DIR), "output")
        args = ['-d', fast5_dir,
                '-r', reference_location,
                '-T', TEMPLATE_MODEL,
                '-o', output_tmp_dir, ]

        singleNuclProb.main(args)

        in_file_count = len(glob.glob(os.path.join(fast5_dir, "*.fast5")))
        output_files = glob.glob(os.path.join(output_tmp_dir, "*.tsv"))
        out_file_count = len(output_files)
        assert out_file_count == in_file_count, "Expected {} output files, got {}".format(in_file_count, out_file_count)


    def test_1D_reads(self):
        oneD_reads = os.path.join(SIGNALALIGN_ROOT, "tests/minion_test_reads/1D/")
        ecoli_ref = os.path.join(SIGNALALIGN_ROOT, "tests/test_sequences/E.coli_K12.fasta")
        self.run_single_nucl_prob(oneD_reads, ecoli_ref)


def main():
    testSuite = unittest.TestSuite()
    testSuite.addTest(SingleNuclProbsTest('test_1D_reads'))

    testRunner = unittest.TextTestRunner(verbosity=2)
    testRunner.run(testSuite)

if __name__ == '__main__':
    main()
