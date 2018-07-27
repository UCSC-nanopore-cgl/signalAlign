#!/usr/bin/env python
"""Test BuildModels.py"""
########################################################################
# File: test_trainModels.py
#  executable: test_trainModels.py
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
from signalalign.train.trainModels import *
from signalalign.utils.fileHandlers import FolderHandler
from signalalign.hiddenMarkovModel import HmmModel
from py3helpers.utils import captured_output, load_json, save_json, count_lines_in_file


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

        cls.model = HmmModel(model_file=cls.model_file)
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
            "singleLevelPriorEcoli": 13,
            "singleLevelFixedCanonical": 14
        }
        cls.test_hdp_training_data = os.path.join(cls.HOME, "tests/test_hdp/test_hdp_alignment.txt")
        cls.one_file_dir = os.path.join(cls.HOME, "tests/minion_test_reads/one_R9_canonical_ecoli")
        cls.config_file = os.path.join(cls.HOME, "tests/trainModels-config.json")
        cls.default_args = create_dot_dict(load_json(cls.config_file))
        cls.default_args.path_to_bin = cls.path_to_bin
        cls.default_args.output_dir = cls.path_to_bin
        cls.default_args.samples[0].fast5_dirs = [cls.one_file_dir]
        cls.default_args.samples[0].bwa_reference = cls.ecoli_reference
        cls.r9_complement_model_file_acgt = os.path.join(cls.HOME, "models/testModelR9_5mer_acgt_complement.model")
        cls.r9_template_model_file_acgt = os.path.join(cls.HOME, "models/testModelR9_5mer_acgt_template.model")

        cls.default_args.complement_hmm_model = cls.r9_complement_model_file_acgt
        cls.default_args.template_hmm_model = cls.r9_template_model_file_acgt

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
                                              probability_threshold=0,
                                              kmers_from_reference=False)
            working_folder = FolderHandler()
            working_folder.open_folder(os.path.join(tempdir, "test_dir"))
            sample = SignalAlignSample(working_folder=working_folder, **test_args)
            sample.analysis_files = [self.assignment_file, self.assignment_file]
            out_path = CreateHdpTrainingData([sample], test_out, template=True, complement=False, verbose=False).write_hdp_training_file()
            n_lines = count_lines_in_file(out_path)
            self.assertEqual(n_lines, 3182)
            with open(out_path, 'r') as fh1, open(self.test_hdp_training_data, 'r') as fh2:
                self.assertEqual(sorted(list(fh1)), sorted(list(fh2)))

    def test_get_sample_kmers(self):
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
                                              probability_threshold=0,
                                              kmers_from_reference=False)
            working_folder = FolderHandler()
            working_folder.open_folder(os.path.join(tempdir, "test_dir"))
            sample = SignalAlignSample(working_folder=working_folder, **test_args)
            sample.analysis_files = [self.assignment_file, self.assignment_file]
            hdp_data_handle = CreateHdpTrainingData([sample], test_out, template=True, complement=False, verbose=False)
            kmers = hdp_data_handle.get_sample_kmers(sample)
            self.assertEqual(kmers, {x for x in all_string_permutations("ATGC", length=6)})
            test_args = create_sa_sample_args(fast5_dirs=[tempdir],
                                              name="some_name",
                                              fw_reference=self.ecoli_reference,
                                              bwa_reference=self.ecoli_reference,
                                              number_of_kmer_assignments=1,
                                              probability_threshold=0,
                                              kmers_from_reference=False,
                                              motifs=[["ATGC", "ETGC"]])
            working_folder = FolderHandler()
            working_folder.open_folder(os.path.join(tempdir, "test_dir"))
            sample = SignalAlignSample(working_folder=working_folder, **test_args)
            sample.analysis_files = [self.assignment_file, self.assignment_file]
            hdp_data_handle = CreateHdpTrainingData([sample], test_out, template=True, complement=False, verbose=False)
            kmers = hdp_data_handle.get_sample_kmers(sample)
            self.assertEqual(kmers, get_motif_kmers(["ATGC", "ETGC"], 6, alphabet="ATGC") |
                             {x for x in all_string_permutations("ATGC", length=6)})
            test_args = create_sa_sample_args(fast5_dirs=[tempdir],
                                              name="some_name",
                                              fw_reference=self.ecoli_reference,
                                              bwa_reference=self.ecoli_reference,
                                              number_of_kmer_assignments=1,
                                              probability_threshold=0,
                                              kmers_from_reference=True,
                                              motifs=[["ATGC", "ETGC"]])
            working_folder = FolderHandler()
            working_folder.open_folder(os.path.join(tempdir, "test_dir"))
            sample = SignalAlignSample(working_folder=working_folder, **test_args)
            sample.analysis_files = [self.assignment_file, self.assignment_file]
            hdp_data_handle = CreateHdpTrainingData([sample], test_out, template=True, complement=False, verbose=False)
            kmers = hdp_data_handle.get_sample_kmers(sample)
            expected_kmers = set()
            for _, _, sequence in read_fasta(self.ecoli_reference):
                expected_kmers |= get_sequence_kmers(sequence, k=6, rev_comp=True)

            self.assertEqual(kmers, get_motif_kmers(["ATGC", "ETGC"], 6, alphabet="ATGC") | expected_kmers)


    def test_get_hdp_type(self):
        for type, integer in self.hdp_types.items():
            self.assertEqual(get_hdp_type(type), integer)

    def test_check_config(self):
        fake_args = create_dot_dict(self.default_args.copy())
        fake_args.path_to_bin = "ducked"
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)

    def test_load_hmm_models(self):
        fake_args = create_dot_dict(self.default_args.copy())
        fake_args.template_hmm_model = "ducked"
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
        fake_args = create_dot_dict(self.default_args.copy())
        fake_args.complement_hmm_model = "ducked"
        fake_args.two_d = True
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
        fake_args = create_dot_dict(self.default_args.copy())
        fake_args.complement_hmm_model = self.r9_complement_model_file
        fake_args.two_d = True
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)

    def test_create_samples(self):
        train_class = TrainSignalAlign(self.default_args)
        for x in train_class.samples:
            self.assertIsInstance(x, SignalAlignSample)

    def test_check_train_hdp_config(self):
        # wrong hdp type (acegt)
        r9_acegt = self.r9_template_model_file
        r9_acegot = os.path.join(self.HOME, "models/testModelR9_5mer_acegot_template.model")
        r9_acgt = os.path.join(self.HOME, "models/testModelR9_5mer_acgt_template.model")
        r9_acegit = os.path.join(self.HOME, "models/testModelR9_5mer_acegit_template.model")

        fake_args = create_dot_dict(self.default_args.copy())
        # test Canonical options
        fake_args.hdp_args.hdp_type = "singleLevelFixedCanonical"
        fake_args.template_hmm_model = r9_acegt
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
        fake_args.template_hmm_model = r9_acegot
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
        fake_args.template_hmm_model = r9_acegit
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
        fake_args.template_hmm_model = r9_acgt
        self.assertIsInstance(TrainSignalAlign(fake_args), TrainSignalAlign)
        # wrong hdp type (acegt)
        fake_args.hdp_args.hdp_type = "singleLevelPrior2"
        fake_args.template_hmm_model = r9_acgt
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
        fake_args.template_hmm_model = r9_acegot
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
        fake_args.template_hmm_model = r9_acegit
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
        fake_args.template_hmm_model = r9_acegt
        self.assertIsInstance(TrainSignalAlign(fake_args), TrainSignalAlign)
        # wrong hdp type (acegit) and needs to be 2D
        fake_args.hdp_args.hdp_type = "multisetPriorEcoli"
        fake_args.template_hmm_model = r9_acegit
        fake_args.complement_hmm_model = r9_acegit
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
        fake_args.two_d = True
        self.assertIsInstance(TrainSignalAlign(fake_args), TrainSignalAlign)
        fake_args.template_hmm_model = r9_acgt
        fake_args.complement_hmm_model = r9_acgt
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
        fake_args.template_hmm_model = r9_acegot
        fake_args.complement_hmm_model = r9_acegot
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
        fake_args.template_hmm_model = r9_acegt
        fake_args.complement_hmm_model = r9_acegt
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
        # wrong hdp type (acegot)
        fake_args.hdp_args.hdp_type = "singleLevelFixed"
        fake_args.template_hmm_model = r9_acgt
        fake_args.complement_hmm_model = r9_acgt
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
        fake_args.template_hmm_model = r9_acegt
        fake_args.complement_hmm_model = r9_acegt
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
        fake_args.template_hmm_model = r9_acegit
        fake_args.complement_hmm_model = r9_acegit
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
        fake_args.template_hmm_model = r9_acegot
        fake_args.complement_hmm_model = r9_acegot
        self.assertIsInstance(TrainSignalAlign(fake_args), TrainSignalAlign)
        # test built_alignments
        fake_args.hdp_args.built_alignments = "path_to_bs"
        fake_args.training.expectation_maximization = True
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
        fake_args.training.expectation_maximization = False
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
        fake_args.hdp_args.built_alignments = self.assignment_file
        self.assertIsInstance(TrainSignalAlign(fake_args), TrainSignalAlign)

    def test_check_train_transitions_config(self):
        # wrong job_count type
        fake_args = create_dot_dict(self.default_args.copy())
        fake_args.job_count = self.r9_template_model_file
        fake_args.training.transitions = True
        self.assertRaises(AssertionError, TrainSignalAlign, fake_args)

    def test_canonical_expectation_maximization_training(self):
        with captured_output() as (out, err):
            with tempfile.TemporaryDirectory() as tempdir:
                fake_args = create_dot_dict(self.default_args.copy())
                fake_args.output_dir = tempdir
                fake_args.job_count = 1
                fake_args.training.transitions = True
                fake_args.training.normal_emissions = True
                fake_args.training.hdp_emissions = True
                fake_args.training.expectation_maximization = True
                fake_args.training.em_iterations = 3
                fake_args.hdp_args.gibbs_samples = 100
                fake_args.hdp_args.burnin_multiplier = 0
                fake_args.hdp_args.number_of_assignments = 1
                # Test EM training 3 rounds
                TrainSignalAlign(fake_args).expectation_maximization_training()

    def test_hdp_and_transition_training(self):
        with captured_output() as (out, err):
            with tempfile.TemporaryDirectory() as tempdir:
                fake_args = create_dot_dict(self.default_args.copy())
                fake_args.output_dir = tempdir
                fake_args.training.transitions = True
                fake_args.training.normal_emissions = True
                fake_args.training.hdp_emissions = True
                fake_args.training.expectation_maximization = False
                fake_args.training.em_iterations = 3
                fake_args.hdp_args.gibbs_samples = 100
                fake_args.hdp_args.burnin_multiplier = 0
                fake_args.hdp_args.number_of_assignments = 1
                # Test hdp AND transition training worked
                TrainSignalAlign(fake_args).expectation_maximization_training()

    def test_transition_training(self):
        # with captured_output() as (out, err):
        with tempfile.TemporaryDirectory() as tempdir:
            fake_args = create_dot_dict(self.default_args.copy())
            fake_args.output_dir = tempdir
            fake_args.training.transitions = True
            fake_args.training.normal_emissions = True
            fake_args.training.hdp_emissions = False
            fake_args.training.expectation_maximization = False
            # with captured_output() as (out, err):
            # Test transitions training worked
            TrainSignalAlign(fake_args).expectation_maximization_training()

    def test_hdp_training(self):
        with captured_output() as (out, err):
            with tempfile.TemporaryDirectory() as tempdir:
                fake_args = create_dot_dict(self.default_args.copy())
                fake_args.output_dir = tempdir
                fake_args.training.transitions = False
                fake_args.training.normal_emissions = True
                fake_args.training.hdp_emissions = True
                fake_args.training.expectation_maximization = False
                fake_args.training.em_iterations = 3
                fake_args.hdp_args.gibbs_samples = 100
                fake_args.hdp_args.burnin_multiplier = 0
                fake_args.hdp_args.number_of_assignments = 1
                # Test hdp training worked
                TrainSignalAlign(fake_args).expectation_maximization_training()

                fake_args.training.hdp_emissions = False
                # raise error when no training is selected
                self.assertRaises(AssertionError, TrainSignalAlign(fake_args).expectation_maximization_training)

    def test_methylated_transitions_training(self):
        with captured_output() as (out, err):
            with tempfile.TemporaryDirectory() as tempdir:
                fake_args = create_dot_dict(self.default_args.copy())
                fake_args.output_dir = tempdir
                fake_args.samples.append(create_dot_dict(fake_args.samples[0]))
                fake_args.samples[1].motifs = [["CCAGG", "CEAGG"], ["CCTGG", "CETGG"]]
                fake_args.samples[1].name = "methylated"
                fake_args.complement_hmm_model = self.r9_complement_model_file
                fake_args.template_hmm_model = self.r9_template_model_file

                fake_args.job_count = 1
                fake_args.training.transitions = True
                fake_args.training.normal_emissions = False
                fake_args.training.hdp_emissions = False
                fake_args.training.expectation_maximization = False
                # Test EM training 3 rounds
                template_hmm_model_path, complement_hmm_model_path, template_hdp_model_path, complement_hdp_model_path = \
                    TrainSignalAlign(fake_args).expectation_maximization_training()

    def test_methylated_hdp_training(self):
        with captured_output() as (out, err):
            with tempfile.TemporaryDirectory() as tempdir:
                fake_args = create_dot_dict(self.default_args.copy())
                fake_args.output_dir = tempdir
                fake_args.samples.append(create_dot_dict(fake_args.samples[0]))
                fake_args.samples[1].motifs = [["CCAGG", "CEAGG"], ["CCTGG", "CETGG"]]
                fake_args.samples[1].name = "methylated"
                fake_args.complement_hmm_model = self.r9_complement_model_file
                fake_args.template_hmm_model = self.r9_template_model_file

                fake_args.job_count = 1
                fake_args.training.transitions = False
                fake_args.training.normal_emissions = False
                fake_args.training.hdp_emissions = True
                fake_args.training.expectation_maximization = False
                fake_args.training.em_iterations = 3
                fake_args.hdp_args.hdp_type = "singleLevelPrior2"
                fake_args.hdp_args.gibbs_samples = 100
                fake_args.hdp_args.burnin_multiplier = 0
                fake_args.hdp_args.number_of_assignments = 1
                fake_args.two_d = True

                # Test EM training 3 rounds
                template_hmm_model_path, complement_hmm_model_path, template_hdp_model_path, complement_hdp_model_path = \
                    TrainSignalAlign(fake_args).expectation_maximization_training()



if __name__ == '__main__':
    unittest.main()
