#!/usr/bin/env python3

import sys
import unittest
import glob
import os
import shutil
import pandas as pd
import numpy as np
from subprocess import call
from signalalign.utils.sequenceTools import getFastaDictionary

SIGNALALIGN_ROOT = '/'.join(os.path.abspath(__file__).split("/")[:-4])
ZYMO_C_READS = os.path.join(SIGNALALIGN_ROOT, "tests/minion_test_reads/C/")
ZYMO_REFERENCE = os.path.join(SIGNALALIGN_ROOT, "tests/test_sequences/zymo_sequence.fasta")
BIN_PATH = os.path.join(SIGNALALIGN_ROOT, "bin")


def parse_alignment_full(alignment_file):
    data = pd.read_table(alignment_file, usecols=(1, 2, 4, 5, 9, 12, 13),
                         dtype={'ref_pos': np.int64,
                                'ref_kmer': np.str,
                                'strand': np.str,
                                'event_index': np.int64,
                                'kmer': np.str,
                                'posterior_prob': np.float64,
                                'event_mean': np.float64},
                         header=None,
                         names=['ref_pos', 'ref_kmer', 'strand', 'event_index', 'kmer', 'posterior_prob', 'event_mean'])
    return data


class LibTest(unittest.TestCase):
    def test_signalAlign_library(self):
        current_wd = os.getcwd()
        command = os.path.join(BIN_PATH, "signalAlignLibTests")
        os.chdir(BIN_PATH)
        result = call(command, shell=True, bufsize=-1, stdout=sys.stdout, stderr=sys.stderr)
        self.assertTrue(result == 0, "signalAlign Library Tests Fail")
        os.chdir(current_wd)


class signalAlignLibTests(unittest.TestCase):
    def setUp(self):
        self.work_dir = "./signalAlign_pylibTest/"
        os.makedirs(self.work_dir)

    def tearDown(self):
        shutil.rmtree(self.work_dir)


class SignalAlignAlignmentTest(unittest.TestCase):
    def setUp(self):
        self.current_wd = os.getcwd()
        os.chdir(BIN_PATH)
        if os.path.exists("./signalAlign_unittest/"):
            shutil.rmtree("./signalAlign_unittest/")
        os.makedirs("./signalAlign_unittest/")

    def tearDown(self):
        shutil.rmtree("./signalAlign_unittest/")
        os.chdir(self.current_wd)

    def check_alignments(self, true_alignments, reads, reference, kmer_length, contig_name, extra_args=None):

        def get_kmer(start):
            kmer = referece_sequence[start:start + kmer_length]
            if type(kmer) is str:
                return kmer
            else:
                return bytes.decode(kmer)

        assert len(glob.glob(reads + "*.fast5")) > 0, "Didn't find test MinION reads"
        assert os.path.isfile(reference), "Didn't find reference sequence"
        run_signal_align = os.path.join(BIN_PATH, "runSignalAlign")
        alignment_command = "{runsignalalign} -d={reads} --bwa_reference={ref} -smt=threeState -o={testDir} " \
                            "".format(runsignalalign=run_signal_align, reads=reads, ref=reference,
                                      testDir="./signalAlign_unittest/")
        if extra_args is not None:
            alignment_command += extra_args

        # null_output = open(os.devnull, 'w')
        # result = call(alignment_command, shell=True, bufsize=-1, stdout=null_output, stderr=null_output)

        result = call(alignment_command, shell=True, bufsize=-1)

        self.assertTrue(result == 0, "error running signalAlign alignments command was {}"
                                     "".format(alignment_command))
        test_alignments = glob.glob("./signalAlign_unittest/tempFiles_alignment/*.tsv")

        referece_sequence = getFastaDictionary(reference)[contig_name]

        self.assertTrue(len(test_alignments) == len(glob.glob(true_alignments + "*.tsv")),
                        "Didn't make all alignments got {got} should be {should}".format(got=len(test_alignments),
                                                                                         should=len(glob.glob(
                                                                                             true_alignments + "*.tsv"))))

        for alignment in test_alignments:
            alignment_file = alignment.split("/")[-1]
            expected = parse_alignment_full(os.path.join(true_alignments, alignment_file))
            obs = parse_alignment_full(alignment)
            self.assertEqual(len(obs), len(expected))
            for row in obs.itertuples():
                ref_pos = row[1]
                obs_kmer = row[2]
                strand = row[3]
                exp_kmer = get_kmer(ref_pos)
                self.assertEqual(obs_kmer, exp_kmer, msg="kmer at index {idx} on strand {strand} is {obs} should be "
                                                         "{exp}, file {f}".format(idx=ref_pos,
                                                                                  strand=strand,
                                                                                  obs=obs_kmer,
                                                                                  exp=exp_kmer,
                                                                                  f=alignment))

    # def test_zymo_reads(self):
    #     zymo_true_alignments = os.path.join(SIGNALALIGN_ROOT,
    #                                         'tests/test_alignments/zymo_C_test_alignments_sm3/tempFiles_alignment/')
    #     self.check_alignments(true_alignments=zymo_true_alignments,
    #                           reads=ZYMO_C_READS,
    #                           reference=ZYMO_REFERENCE,
    #                           kmer_length=6,
    #                           contig_name="ZYMO",
    #                           extra_args="--2d --debug")

    def test_pUC_r9_reads_5mer(self):
        pUC_true_alignments = os.path.join(SIGNALALIGN_ROOT, "tests/test_alignments/pUC_5mer_tempFiles_alignment/")
        pUC_reads = os.path.join(SIGNALALIGN_ROOT, "tests/minion_test_reads/pUC/")
        pUC_reference = os.path.join(SIGNALALIGN_ROOT, "tests/test_sequences/pUC19_SspI.fa")
        self.check_alignments(true_alignments=pUC_true_alignments,
                              reads=pUC_reads,
                              reference=pUC_reference,
                              kmer_length=5,
                              contig_name="pUC19",
                              extra_args="-T=../models/testModelR9_5mer_acegot_template.model "
                                         "-C=../models/testModelR9_5mer_acegot_complement.model "
                                         "--2d ")

    def test_pUC_r9_reads_6mer(self):
        pUC_true_alignments = os.path.join(SIGNALALIGN_ROOT, "tests/test_alignments/pUC_6mer_tempFiles_alignment/")
        pUC_reads = os.path.join(SIGNALALIGN_ROOT, "tests/minion_test_reads/pUC/")
        pUC_reference = os.path.join(SIGNALALIGN_ROOT, "tests/test_sequences/pUC19_SspI.fa")
        self.check_alignments(true_alignments=pUC_true_alignments,
                              reads=pUC_reads,
                              reference=pUC_reference,
                              kmer_length=6,
                              contig_name="pUC19",
                              extra_args="--2d ")

    def test_Ecoli1D_reads_5mer(self):
        pUC_true_alignments = os.path.join(SIGNALALIGN_ROOT, "tests/test_alignments/ecoli1D_test_alignments_sm3/")
        pUC_reads = os.path.join(SIGNALALIGN_ROOT, "tests/minion_test_reads/1D/")
        pUC_reference = os.path.join(SIGNALALIGN_ROOT, "tests/test_sequences/E.coli_K12.fasta")
        self.check_alignments(true_alignments=pUC_true_alignments,
                              reads=pUC_reads,
                              reference=pUC_reference,
                              kmer_length=5,
                              contig_name="gi_ecoli",
                              extra_args="-T=../models/testModelR9p4_5mer_acegt_template.model ")

    def test_RNA_edge_alignments_reads_5mer(self):
        edge_case_true_alignments = os.path.join(SIGNALALIGN_ROOT,
                                                 "tests/test_alignments/RNA_edge_case_tempFiles_alignment/")
        edge_case_reads = os.path.join(SIGNALALIGN_ROOT, "tests/minion_test_reads/RNA_edge_cases/")
        edge_case_reference = os.path.join(SIGNALALIGN_ROOT, "tests/test_sequences/fake_rna_reversed.fa")
        self.check_alignments(true_alignments=edge_case_true_alignments,
                              reads=edge_case_reads,
                              reference=edge_case_reference,
                              kmer_length=5,
                              contig_name="rna_fake_reversed",
                              extra_args="-T=../models/testModelR9p4_5mer_acgt_RNA.model ")


class signalAlign_EM_test(unittest.TestCase):

    def setUp(self):
        self.current_wd = os.getcwd()
        os.chdir(BIN_PATH)
        if os.path.exists("./signalAlign_unittest/"):
            shutil.rmtree("./signalAlign_unittest/")
        os.makedirs("./signalAlign_unittest/")

    def tearDown(self):
        shutil.rmtree("./signalAlign_unittest/")
        os.chdir(self.current_wd)

    def test_EM(self):
        trainModels_path = os.path.join(BIN_PATH, "trainModels")
        config_path = os.path.join(SIGNALALIGN_ROOT, "tests/trainModels-config.yaml")
        em_command = "{train_models} run --config={config}".format(train_models=trainModels_path, config=config_path)
        # null_output = open(os.devnull, 'w')
        # result = call(em_command, shell=True, bufsize=-1, stdout=null_output, stderr=null_output)
        result = call(em_command, shell=True, bufsize=-1)

        self.assertTrue(result == 0, "error running signalAlign alignments command was {}"
                                     "".format(em_command))


def add_all_tests_to_Suite(test_suite, test_class):
    """Add all tests from a testCase class to testSuite
    :param test_suite: instance of unittest.TestSuite
    :param test_class: instance of instance of unittest.TestCase
    """
    # make sure they are correct instances with helpful assertion
    assert isinstance(test_suite, unittest.TestSuite), "test_Suite must be an instance of unittest.TestSuite()"
    assert issubclass(test_class, unittest.TestCase), "test_class must be an instance of unittest.TestCase()"
    # grab all functions with test
    [test_suite.addTest(test_class(func)) for func in dir(test_class) if
     callable(getattr(test_class, func)) and not func.find("test")]

    return test_suite


def main():
    testSuite = unittest.TestSuite()
    testSuite.addTest(LibTest('test_signalAlign_library'))
    testSuite.addTest(SignalAlignAlignmentTest('test_pUC_r9_reads_5mer'))
    testSuite.addTest(SignalAlignAlignmentTest('test_pUC_r9_reads_6mer'))
    testSuite.addTest(SignalAlignAlignmentTest('test_Ecoli1D_reads_5mer'))
    testSuite.addTest(SignalAlignAlignmentTest('test_RNA_edge_alignments_reads_5mer'))
    testSuite.addTest(signalAlign_EM_test('test_EM'))
    # testSuite.addTest(SignalAlignAlignmentTest('test_zymo_reads'))

    testRunner = unittest.TextTestRunner(verbosity=2)
    testRunner.run(testSuite)


if __name__ == '__main__':
    main()
