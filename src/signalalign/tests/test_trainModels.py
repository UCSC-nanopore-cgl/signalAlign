#!/usr/bin/env python
"""Test trainModels.py"""
########################################################################
# File: test_trainModels.py
#  executable: test_trainModels.py
#
# Author: Andrew Bailey
# History: 5/21/18 Created
########################################################################


import unittest
import tempfile
from signalalign.signalAlignment import create_sa_sample_args
from signalalign.train.trainModels import *
from signalalign.utils.fileHandlers import FolderHandler
from signalalign.hiddenMarkovModel import HmmModel
from signalalign.mixture_model import get_motif_kmer_pairs
from py3helpers.utils import captured_output, load_json, time_it
from py3helpers.seq_tools import ReferenceHandler, ReverseComplement


class TrainSignalAlignTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TrainSignalAlignTest, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.reference = os.path.join(cls.HOME, "tests/test_sequences/pUC19_SspI_Zymo.fa")
        cls.ecoli_reference = os.path.join(cls.HOME, "tests/test_sequences/E.coli_K12.fasta")
        cls.ecoli_bam = os.path.join(cls.HOME, "tests/minion_test_reads/canonical_ecoli_R9/canonical_ecoli.bam")
        cls.ecoli_readdb = os.path.join(cls.HOME, "tests/minion_test_reads/canonical_ecoli_R9/canonical_ecoli.readdb")
        cls.fast5_dir = os.path.join(cls.HOME, "tests/minion_test_reads/canonical_ecoli_R9")
        cls.files = [
            "miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_"
            "strand.fast5",
            "miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch101_"
            "read456_strand.fast5",
            "miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch101_"
            "read544_strand1.fast5",
            "miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch103_"
            "read333_strand1.fast5"]
        cls.fast5_paths = [os.path.join(cls.fast5_dir, f) for f in os.listdir(cls.fast5_dir)
                           if os.path.isfile(os.path.join(cls.fast5_dir, f))]
        cls.model_file = os.path.join(cls.HOME, "models/testModelR9p4_5mer_acgt_RNA.model")
        cls.r9_complement_model_file = os.path.join(cls.HOME, "models/testModelR9_acegt_complement.model")
        cls.r9_template_model_file = os.path.join(cls.HOME, "models/testModelR9_acegt_template.model")

        cls.model = HmmModel(ont_model_file=cls.model_file)
        cls.expectation_file = \
            os.path.join(cls.HOME,
                         "tests/test_expectation_files/4f9a316c-8bb3-410a-8cfc-026061f7e8db.template.expectations.tsv")
        cls.assignment_file = \
            os.path.join(cls.HOME,
                         "tests/test_assignment_files/d6160b0b-a35e-43b5-947f-adaa1abade28.sm.assignments.tsv")

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
            "singleLevelFixedCanonical": 14,
            "singleLevelFixedM6A": 15
        }
        cls.test_hdp_training_data = os.path.join(cls.HOME, "tests/test_hdp/test_hdp_alignment.txt")
        cls.one_file_dir = os.path.join(cls.HOME, "tests/minion_test_reads/one_R9_canonical_ecoli")
        cls.config_file = os.path.join(cls.HOME, "tests/trainModels-config.json")
        cls.default_args = create_dot_dict(load_json(cls.config_file))
        cls.default_args.path_to_bin = cls.path_to_bin
        cls.default_args.output_dir = cls.path_to_bin
        cls.default_args.debug = False
        cls.default_args.samples[0].fast5_dirs = [cls.one_file_dir]
        cls.default_args.samples[0].bwa_reference = cls.ecoli_reference
        cls.default_args.samples[0].alignment_file = cls.ecoli_bam
        cls.default_args.samples[0].readdb = cls.ecoli_readdb
        cls.r9_complement_model_file_acgt = os.path.join(cls.HOME, "models/testModelR9_5mer_acgt_complement.model")
        cls.r9_template_model_file_acgt = os.path.join(cls.HOME, "models/testModelR9_5mer_acgt_template.model")

        cls.default_args.complement_hmm_model = cls.r9_complement_model_file_acgt
        cls.default_args.template_hmm_model = cls.r9_template_model_file_acgt

    def test_parse_assignment_file(self):
        with captured_output() as (_, _):
            assignments = parse_assignment_file(self.assignment_file)
            # "kmer", "strand", "level_mean", "prob"
            self.assertEqual(len(assignments['kmer']), len(assignments['strand']))
            self.assertEqual(len(assignments['level_mean']), len(assignments['prob']))
            self.assertEqual(len(assignments['kmer']), len(assignments['prob']))

    def test_make_master_assignment_table(self):
        with captured_output() as (_, _):
            pandas_table = make_master_assignment_table([self.assignment_file, self.assignment_file])
            assignments1 = parse_assignment_file(self.assignment_file)
            self.assertEqual(len(pandas_table), 2 * len(assignments1))

    def test_multiprocess_make_kmer_assignment_tables(self):
        with captured_output() as (_, _):
            kmers = get_kmers(6, alphabet="ATGC")
            data2, time2 = time_it(multiprocess_make_kmer_assignment_tables, [self.assignment_file], kmers,
                                   set("t"), 0.0, False, False, 100000, 2)
            for x in kmers:
                kmer_data = data2.loc[data2['kmer'] == x]
                self.assertSequenceEqual(list(kmer_data["prob"]), sorted(kmer_data["prob"], reverse=True))

    def test_generate_build_alignments3(self):
        with captured_output() as (_, _):
            kmers = get_kmers(6, alphabet="ATGC")
            data_files = [self.assignment_file]

            sample_assignment_table = get_assignment_table(self.assignment_file, 0.0, False)
            data1 = generate_build_alignments(sample_assignment_table, kmers, 10, ["t"], False)

            data2, _ = time_it(multiprocess_make_kmer_assignment_tables,
                               data_files, kmers,
                               set("t"), 0.0, False, False, 10, 8)

            # get kmers associated with each sample
            self.assertEqual(len(data2), len(data1))

    def test_sort_dataframe_wrapper(self):
        with captured_output() as (_, _):
            kmers = get_kmers(6, alphabet="ATGC")
            kmers = list(kmers)
            kmer = "ATTTTT"
            index = kmers.index(kmer)
            data_table = get_assignment_kmer_tables(self.assignment_file, kmers, 0, False)
            kmer_data = data_table[index]
            data = sort_dataframe_wrapper(kmer_data, kmer, max_assignments=10, verbose=False, strands=('t', 'c'))

            self.assertSequenceEqual(list(data["prob"]), sorted(data["prob"], reverse=True))

    def test_make_alignment_line(self):
        with captured_output() as (_, _):
            make_alignment_line(strand='t', kmer="ATGC", prob=0.1, event=23.2)
            self.assertRaises(AssertionError, make_alignment_line, strand='d', kmer="ATGC", prob=0.1, event=23.2)
            self.assertRaises(AssertionError, make_alignment_line, strand='t', kmer=1, prob=0.1, event=23.2)
            self.assertRaises(AssertionError, make_alignment_line, strand='t', kmer="ATGC", prob=1, event=23.2)
            self.assertRaises(AssertionError, make_alignment_line, strand='t', kmer="ATGC", prob=0.1, event=23)

    def test_generate_build_alignments(self):
        with captured_output() as (_, _):
            d = dict(kmer=["ATGC", "AAAA", "ATGC", "AAAA"], strand=['t', 'c', 't', 'c'],
                     level_mean=[99.9, 89.9, 99.9, 89.9], prob=[0.01, 0.99, 0.01, 0.99])
            assignments = pd.DataFrame(d)

            line_generator = generate_build_alignments(assignments, kmer_list=["ATGC", "AAAA"], max_assignments=2,
                                                       strands=['t', 'c'], verbose=False)
            self.assertEqual(4, len(line_generator))
            line_generator = generate_build_alignments(assignments, kmer_list=["ATGC"], max_assignments=2,
                                                       strands=['t', 'c'], verbose=False)
            self.assertEqual(2, len(line_generator))

            line_generator = generate_build_alignments(assignments, kmer_list=["ATGC", "AAAA"], max_assignments=1,
                                                       strands=['t', 'c'], verbose=False)
            self.assertEqual(2, len(line_generator))
            line_generator = generate_build_alignments(assignments, kmer_list=["ATGC", "AAAA"], max_assignments=2,
                                                       strands=['t'], verbose=False)
            self.assertEqual(2, len(line_generator))
            line_generator = generate_build_alignments(assignments, kmer_list=["ATGC", "AAAA"], max_assignments=2,
                                                       strands=['c'], verbose=False)
            self.assertEqual(2, len(line_generator))
            line_generator = generate_build_alignments(assignments, kmer_list=["ATGC", "AAAA"], max_assignments=2,
                                                       strands=['t', 'c'], verbose=False)
            self.assertEqual(4, len(line_generator))
            self.assertRaises(AssertionError, generate_build_alignments, assignments, kmer_list=["ATGC", "AAAA"],
                              max_assignments=2,
                              strands=[], verbose=False)

    def test_CreateHdpTrainingData(self):
        with captured_output() as (_, _):
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
                out_path = CreateHdpTrainingData([sample], test_out, template=True, complement=False,
                                                 verbose=False).write_hdp_training_file(verbose=False)
                n_lines = count_lines_in_file(out_path)
                self.assertEqual(n_lines, 3182)

    def test_get_sample_kmers(self):
        with captured_output() as (_, _):
            with tempfile.TemporaryDirectory() as tempdir:
                # create fast5 dir
                test_fast5 = os.path.join(tempdir, "test.fast5")
                copyfile(self.fast5_paths[0], test_fast5)
                # create file of file names
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
                hdp_data_handle = CreateHdpTrainingData([sample], test_out, template=True, complement=False,
                                                        verbose=False)
                kmers = hdp_data_handle.get_sample_kmers(sample, kmer_len=6)
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
                hdp_data_handle = CreateHdpTrainingData([sample], test_out, template=True, complement=False,
                                                        verbose=False)
                kmers = hdp_data_handle.get_sample_kmers(sample, kmer_len=6)
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
                hdp_data_handle = CreateHdpTrainingData([sample], test_out, template=True, complement=False,
                                                        verbose=False)
                kmers = hdp_data_handle.get_sample_kmers(sample, kmer_len=6)
                expected_kmers = set()
                for _, _, sequence in read_fasta(self.ecoli_reference):
                    expected_kmers |= get_sequence_kmers(sequence, k=6, rev_comp=True)

                self.assertEqual(kmers, get_motif_kmers(["ATGC", "ETGC"], 6, alphabet="ATGC") | expected_kmers)

    def test_get_hdp_type(self):
        with captured_output() as (_, _):
            for hdp_type, integer in self.hdp_types.items():
                self.assertEqual(get_hdp_type(hdp_type), integer)

    def test_check_config(self):
        with captured_output() as (_, _):
            fake_args = create_dot_dict(self.default_args.copy())
            fake_args.path_to_bin = "ducked"
            self.assertRaises(AssertionError, TrainSignalAlign, fake_args)

    def test_load_hmm_models(self):
        with captured_output() as (_, _):
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
        with captured_output() as (_, _):
            train_class = TrainSignalAlign(self.default_args)
            for x in train_class.samples:
                self.assertIsInstance(x, SignalAlignSample)

    def test_check_train_hdp_config(self):
        with captured_output() as (_, _):
            # wrong hdp type (acegt)
            # r9_acegt = self.r9_template_model_file
            # r9_acegot = os.path.join(self.HOME, "models/testModelR9_5mer_acegot_template.model")
            # r9_acgt = os.path.join(self.HOME, "models/testModelR9_5mer_acgt_template.model")
            # r9_acegit = os.path.join(self.HOME, "models/testModelR9_5mer_acegit_template.model")

            fake_args = create_dot_dict(self.default_args.copy())
            # # test Canonical options
            # fake_args.hdp_args.hdp_type = "singleLevelFixedCanonical"
            # fake_args.template_hmm_model = r9_acegt
            # self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
            # fake_args.template_hmm_model = r9_acegot
            # self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
            # fake_args.template_hmm_model = r9_acegit
            # self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
            # fake_args.template_hmm_model = r9_acgt
            # self.assertIsInstance(TrainSignalAlign(fake_args), TrainSignalAlign)
            # # wrong hdp type (acegt)
            # fake_args.hdp_args.hdp_type = "singleLevelPrior2"
            # fake_args.template_hmm_model = r9_acgt
            # self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
            # fake_args.template_hmm_model = r9_acegot
            # self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
            # fake_args.template_hmm_model = r9_acegit
            # self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
            # fake_args.template_hmm_model = r9_acegt
            # self.assertIsInstance(TrainSignalAlign(fake_args), TrainSignalAlign)
            # # wrong hdp type (acegit) and needs to be 2D
            # fake_args.hdp_args.hdp_type = "multisetPriorEcoli"
            # fake_args.template_hmm_model = r9_acegit
            # fake_args.complement_hmm_model = r9_acegit
            # self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
            # fake_args.two_d = True
            # self.assertIsInstance(TrainSignalAlign(fake_args), TrainSignalAlign)
            # fake_args.template_hmm_model = r9_acgt
            # fake_args.complement_hmm_model = r9_acgt
            # self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
            # fake_args.template_hmm_model = r9_acegot
            # fake_args.complement_hmm_model = r9_acegot
            # self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
            # fake_args.template_hmm_model = r9_acegt
            # fake_args.complement_hmm_model = r9_acegt
            # self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
            # # wrong hdp type (acegot)
            # fake_args.hdp_args.hdp_type = "singleLevelFixed"
            # fake_args.template_hmm_model = r9_acgt
            # fake_args.complement_hmm_model = r9_acgt
            # self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
            # fake_args.template_hmm_model = r9_acegt
            # fake_args.complement_hmm_model = r9_acegt
            # self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
            # fake_args.template_hmm_model = r9_acegit
            # fake_args.complement_hmm_model = r9_acegit
            # self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
            # fake_args.template_hmm_model = r9_acegot
            # fake_args.complement_hmm_model = r9_acegot
            # self.assertIsInstance(TrainSignalAlign(fake_args), TrainSignalAlign)
            # test built_alignments
            fake_args.built_alignments = "path_to_bs"
            fake_args.training.expectation_maximization = True
            self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
            fake_args.training.expectation_maximization = False
            self.assertRaises(AssertionError, TrainSignalAlign, fake_args)
            fake_args.built_alignments = self.assignment_file
            self.assertIsInstance(TrainSignalAlign(fake_args), TrainSignalAlign)

    def test_check_train_transitions_config(self):
        with captured_output() as (_, _):
            # wrong job_count type
            fake_args = create_dot_dict(self.default_args.copy())
            fake_args.job_count = self.r9_template_model_file
            fake_args.training.transitions = True
            self.assertRaises(AssertionError, TrainSignalAlign, fake_args)

    def test_methylated_hdp_training(self):
        with captured_output() as (_, _):
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
                TrainSignalAlign(fake_args).expectation_maximization_training()

    def test_methylated_hdp_training2(self):
        with captured_output() as (_, _):
            with tempfile.TemporaryDirectory() as tempdir:
                fake_args = create_dot_dict(self.default_args.copy())
                fake_args.output_dir = tempdir
                fake_args.samples.append(create_dot_dict(fake_args.samples[0]))
                fake_args.samples[1].motifs = [["CCAGG", "CDAGG"], ["CCTGG", "CDTGG"]]
                fake_args.samples[1].name = "methylated"
                r9_complement_model_file_test = os.path.join(self.HOME,
                                                             "models/testModelR9_5mer_acdgt_complement.model")
                r9_template_model_file_test = os.path.join(self.HOME, "models/testModelR9_5mer_acdgt_template.model")
                fake_args.complement_hmm_model = r9_complement_model_file_test
                fake_args.template_hmm_model = r9_template_model_file_test

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
                TrainSignalAlign(fake_args).expectation_maximization_training()

    def test_generate_build_alignments_given_motifs(self):
        with captured_output() as (_, _):
            assignments_dir = os.path.join(self.HOME, "tests/test_alignments/ecoli1D_test_alignments_sm3")
            positions_file = os.path.join(self.HOME, "tests/test_position_files/CCWGG_ecoli_k12_mg1655.positions")
            all_data = generate_build_alignments_positions([assignments_dir, assignments_dir], positions_file,
                                                           verbose=False)
            rh = ReferenceHandler(self.ecoli_reference)
            #
            # all_data = r
            # ead_in_alignment_file("/Users/andrewbailey/data/ccwgg_new_em_trained_model/
            # mixture_models/high_sd_models/built_alignments_wide_sd_02_28_19.tsv")
            # rh = ReferenceHandler("/Users/andrewbailey/data/references/ecoli/ecoli_k12_mg1655.fa")
            # positions_file =
            # os.path.join("/Users/andrewbailey/data/references/ecoli/CCWGG_ecoli_k12_mg1655_C_C.positions")
            # positions_data = CustomAmbiguityPositions.parseAmbiguityFile(positions_file)
            # data_dir = "/Users/andrewbailey/data/ccwgg_new_em_trained_model/mixture_models/test_output/canonical"
            # # stranded_positions = positions_data[positions_data["strand"] == "t"]
            # all_data = generate_build_alignments_positions([data_dir], positions_file, verbose=True)
            all_kmer_pairs = set()
            motifs = [["CCAGG", "CEAGG"], ["CCTGG", "CETGG"]]
            for motif in motifs:
                all_kmer_pairs |= set(tuple(row) for row in get_motif_kmer_pairs(motif_pair=motif, k=5))
            all_possible_kmers = set([x[0] for x in all_kmer_pairs])

            rc = ReverseComplement()
            for i, row in all_data.iterrows():
                seq = rh.get_sequence(row["contig"], row["reference_index"] - 9, row["reference_index"] + 9)
                self.assertTrue(row["path_kmer"] in all_possible_kmers)
                if row["path_kmer"] == row["reference_kmer"]:
                    self.assertTrue("CCAGG" in seq or "CCTGG" in seq)
                else:
                    self.assertTrue("CCAGG" in rc.reverse_complement(seq) or "CCTGG" in rc.reverse_complement(seq))


if __name__ == '__main__':
    unittest.main()
