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

    def run_single_nucl_prob(self, fast5_glob, reference_location):
        output_tmp_dir = os.path.join(os.path.abspath(SingleNuclProbsTest.WORK_DIR), "output")
        args = ['-g', fast5_glob,
                '-r', reference_location,
                '-T', TEMPLATE_MODEL,
                '-o', output_tmp_dir,
                '--step_size', '5']

        singleNuclProb.main(args)

        in_file_count = len(glob.glob(fast5_glob))
        output_files = glob.glob(os.path.join(output_tmp_dir, "*.tsv"))
        out_file_count = len(output_files)
        assert out_file_count == in_file_count, "Expected {} output files, got {}".format(in_file_count, out_file_count)


    def test_1D_reads(self):
        oneD_reads = os.path.join(SIGNALALIGN_ROOT, "tests/minion_test_reads/1D/LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch112_read108_strand.fast5")
        ecoli_ref = os.path.join(SIGNALALIGN_ROOT, "tests/test_sequences/E.coli_K12.fasta")
        self.run_single_nucl_prob(oneD_reads, ecoli_ref)


def main():
    testSuite = unittest.TestSuite()
    testSuite.addTest(SingleNuclProbsTest('test_1D_reads'))

    testRunner = unittest.TextTestRunner(verbosity=2)
    testRunner.run(testSuite)

if __name__ == '__main__':
    main()
