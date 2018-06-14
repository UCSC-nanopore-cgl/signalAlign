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
import tempfile
from collections import defaultdict
from shutil import copyfile
from scipy import sparse
from signalalign.signalAlignment import create_sa_sample_args
from signalalign.train.BuildModels import *
from signalalign.utils.fileHandlers import FolderHandler
from signalalign.hiddenMarkovModel import get_model
from py3helpers.utils import captured_output


class TrainSignalAlignTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TrainSignalAlignTest, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.reference = os.path.join(cls.HOME, "tests/test_sequences/pUC19_SspI_Zymo.fa")
        cls.ecoli_reference = os.path.join(cls.HOME, "tests/test_sequences/E.coli_K12.fasta")

        cls.fast5_dir = os.path.join(cls.HOME, "tests/minion_test_reads/canonical_ecoli_R9")
        cls.files = [
            "miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_strand.fast5",
            "miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch101_read456_strand.fast5",
            "miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch101_read544_strand1.fast5",
            "miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch103_read333_strand1.fast5"]
        cls.fast5_paths = [os.path.join(cls.fast5_dir, f) for f in os.listdir(cls.fast5_dir)
                           if os.path.isfile(os.path.join(cls.fast5_dir, f))]
        cls.model_file = os.path.join(cls.HOME, "models/testModelR9p4_5mer_acgt_RNA.model")
        cls.r9_complement_model_file = os.path.join(cls.HOME, "models/testModelR9_acegt_complement.model")
        cls.r9_template_model_file = os.path.join(cls.HOME, "models/testModelR9_acegt_template.model")

        cls.model = get_model(model_type="threeState", model_file=cls.model_file)
        cls.expectation_file = os.path.join(cls.HOME,
                                            "tests/test_expectation_files/4f9a316c-8bb3-410a-8cfc-026061f7e8db.template.expectations.tsv")
        cls.assignment_file = os.path.join(cls.HOME, "tests/test_assignment_files/d6160b0b-a35e-43b5-947f-adaa1abade28.sm.assignments.tsv")

        cls.path_to_bin = os.path.join(cls.HOME, "bin")
        cls.hdp_types = {
            "singleLevelFixed": 0,
            "singleLevelPrior": 1,
            "multisetFixed": 2,
            "multisetPrior": 3,
            "compFixed": 4,
            "compPrior": 5,
            "middleNtsFixed": 6,
            "middleNtsPrior": 7,
            "groupMultisetFixed": 8,
            "groupMultisetPrior": 9,
            "singleLevelPrior2": 10,
            "multisetPrior2": 11,
            "multisetPriorEcoli": 12,
            "singleLevelPriorEcoli": 13
        }


    def test_make_gatc_position_file(self):
        with tempfile.TemporaryDirectory() as tempdir:
            path = os.path.join(tempdir, "test.positions")
            positions_file = make_gatc_position_file(fasta=self.reference,
                                                    outfile=path)
            self.assertEqual(positions_file, path)

    def test_kmer_length_from_model(self):
        self.assertEqual(kmer_length_from_model(self.model_file), self.model.kmer_length)

    def test_parse_assignment_file(self):
        assignments = parse_assignment_file(self.assignment_file)
        # "kmer", "strand", "level_mean", "prob"
        self.assertEqual(len(assignments['kmer']), len(assignments['strand']))
        self.assertEqual(len(assignments['level_mean']), len(assignments['prob']))
        self.assertEqual(len(assignments['kmer']), len(assignments['prob']))

    def test_make_master_assignment_table(self):
        pandas_table = make_master_assignment_table([self.assignment_file, self.assignment_file])
        assignments1 = parse_assignment_file(self.assignment_file)
        self.assertEqual(len(pandas_table), 2*len(assignments1))

    def test_make_alignment_line(self):
        make_alignment_line(strand='t', kmer="ATGC", prob=0.1, event=23.2)
        self.assertRaises(AssertionError, make_alignment_line, strand='d', kmer="ATGC", prob=0.1, event=23.2)
        self.assertRaises(AssertionError, make_alignment_line, strand='t', kmer=1, prob=0.1, event=23.2)
        self.assertRaises(AssertionError, make_alignment_line, strand='t', kmer="ATGC", prob=1, event=23.2)
        self.assertRaises(AssertionError, make_alignment_line, strand='t', kmer="ATGC", prob=0.1, event=23)

    def test_generate_hdp_training_lines(self):
        d = dict(kmer=["ATGC", "AAAA", "ATGC", "AAAA"], strand=['t', 'c', 't', 'c'],
                 level_mean=[99.9, 89.9, 99.9, 89.9], prob=[0.01, 0.99, 0.01, 0.99])
        assignments = pd.DataFrame(d)

        line_generator = CreateHdpTrainingData._generate_hdp_training_lines(assignments, kmer_list=["ATGC", "AAAA"], max_assignments=2,
                                                     strands=['t', 'c'], min_probability=0, verbose=False)
        self.assertEqual(4, len([x for x in line_generator]))
        line_generator = CreateHdpTrainingData._generate_hdp_training_lines(assignments, kmer_list=["ATGC"], max_assignments=2,
                                                                            strands=['t', 'c'], min_probability=0, verbose=False)
        self.assertEqual(2, len([x for x in line_generator]))

        line_generator = CreateHdpTrainingData._generate_hdp_training_lines(assignments, kmer_list=["ATGC", "AAAA"], max_assignments=1,
                                                                            strands=['t', 'c'], min_probability=0, verbose=False)
        self.assertEqual(2, len([x for x in line_generator]))
        line_generator = CreateHdpTrainingData._generate_hdp_training_lines(assignments, kmer_list=["ATGC", "AAAA"], max_assignments=2,
                                                                            strands=['t'], min_probability=0, verbose=False)
        self.assertEqual(2, len([x for x in line_generator]))
        line_generator = CreateHdpTrainingData._generate_hdp_training_lines(assignments, kmer_list=["ATGC", "AAAA"], max_assignments=2,
                                                                            strands=['c'], min_probability=0, verbose=False)
        self.assertEqual(2, len([x for x in line_generator]))
        line_generator = CreateHdpTrainingData._generate_hdp_training_lines(assignments, kmer_list=["ATGC", "AAAA"], max_assignments=2,
                                                                            strands=['t', 'c'], min_probability=0.9, verbose=False)
        self.assertEqual(2, len([x for x in line_generator]))
        line_generator = CreateHdpTrainingData._generate_hdp_training_lines(assignments, kmer_list=["ATGC", "AAAA"], max_assignments=2,
                                                                            strands=[], min_probability=0, verbose=False)

        self.assertRaises(AssertionError, line_generator.__next__)

    def test_CreateHdpTrainingData(self):
        with tempfile.TemporaryDirectory() as tempdir:
            # create fast5 dir
            test_fast5 = os.path.join(tempdir, "test.fast5")
            copyfile(self.fast5_paths[0], test_fast5)
            # create fofn
            test_out = os.path.join(tempdir, "test.hdp.tsv")
            test_args = create_sa_sample_args(fast5_dirs=[tempdir],
                                              name="some_name",
                                              fw_reference=self.ecoli_reference,
                                              bwa_reference=self.ecoli_reference,
                                              number_of_kmer_assignments=1,
                                              probability_threshold= 0)
            working_folder = FolderHandler()
            working_folder.open_folder(os.path.join(tempdir, "test_dir"))
            sample = SignalAlignSample(working_folder=working_folder, **test_args)
            sample.analysis_files = [self.assignment_file, self.assignment_file]
            out_path = CreateHdpTrainingData([sample], test_out, 6, template=True, complement=False, verbose=False).write_hdp_training_file()
            n_lines = count_lines_txt_file(out_path)
            self.assertEqual(n_lines, 3182)
            # out_path = CreateHdpTrainingData([sample], test_out, 6, template=True, complement=False, verbose=False).write_hdp_training_file(bulk)
            # n_lines = count_lines_txt_file(out_path)
            # self.assertEqual(n_lines, 2*count_lines_txt_file(self.assignment_file))


# def test_build_hdp(self):
    #     with tempfile.TemporaryDirectory() as tempdir:
    #             # create fast5 dir
    #         test_fast5 = os.path.join(tempdir, "test.fast5")
    #         copyfile(self.fast5_paths[0], test_fast5)
    #         # create fofn
    #         test_out = os.path.join(tempdir, "test.hdp.tsv")
    #         test_args = create_sa_sample_args(fast5_dirs=[tempdir],
    #                                           name="some_name",
    #                                           fw_reference=self.ecoli_reference,
    #                                           bwa_reference=self.ecoli_reference)
    #         working_folder = FolderHandler()
    #         working_folder.open_folder(os.path.join(tempdir, "test_dir"))
    #         sample = SignalAlignSample(working_folder=working_folder, **test_args)
    #         sample.analysis_files = [self.assignment_file, self.assignment_file]
    #         out_path = CreateHdpTrainingData([sample], test_out, 5, template=True, complement=False, verbose=False).write_hdp_training_file()
    #         print(os.path.exists(working_folder.path), working_folder.path)
    #         hdps = build_hdp(build_alignment_path=out_path,
    #                          template_model=self.r9_template_model_file,
    #                          complement_model=self.r9_complement_model_file,
    #                          hdp_type='multiset',
    #                          outpath=working_folder.path,
    #                          samples=10,
    #                          path_to_bin=self.path_to_bin)

    def test_train_hdp(self):
        with tempfile.TemporaryDirectory() as tempdir:
                # create fast5 dir
            test_fast5 = os.path.join(tempdir, "test.fast5")
            copyfile(self.fast5_paths[0], test_fast5)
            # create fofn
            test_out = os.path.join(tempdir, "test.hdp.tsv")
            test_args = create_sa_sample_args(fast5_dirs=[tempdir],
                                              name="some_name",
                                              fw_reference=self.ecoli_reference,
                                              bwa_reference=self.ecoli_reference)
            working_folder = FolderHandler()
            working_folder.open_folder(os.path.join(tempdir, "test_dir"))
            sample = SignalAlignSample(working_folder=working_folder, **test_args)
            sample.analysis_files = [self.assignment_file, self.assignment_file]
            out_path = CreateHdpTrainingData([sample], test_out, 6, template=True, complement=True, verbose=False).write_hdp_training_file()
            args = create_hdp_training_args(build_alignment=out_path,
                                            template_model=self.r9_template_model_file,
                                            complement_model=self.r9_complement_model_file,
                                            hdp_type='singleLevelPrior2',
                                            outpath=working_folder.path,
                                            samples=10,
                                            thinning=1,
                                            number_of_assignments=2,
                                            path_to_bin=self.path_to_bin,
                                            num_alignments=10,
                                            verbose=True,
                                            twoD=False)
            hdps = train_hdp(args)


    def test_get_hdp_type(self):
        for type, integer in self.hdp_types.items():
            self.assertEqual(get_hdp_type(type), integer)

    def test_get_initial_hdp_args(self):
        args = create_hdp_training_args()
        for hdp_type in self.hdp_types:
            self.assertIsInstance(get_initial_hdp_args(args, hdp_type), str)

    def test_count_lines_txt_file(self):
        self.assertEqual(count_lines_txt_file(self.assignment_file), 17350)


if __name__ == '__main__':
    unittest.main()
