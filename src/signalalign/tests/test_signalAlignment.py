#!/usr/bin/env python
"""Test signalAlignment.py"""
########################################################################
# File: test_signalAlignment.py
#  executable: test_signalAlignment.py
#
# Author: Andrew Bailey
# History: 6/11/18 Created
########################################################################


import unittest
import tempfile
from shutil import copyfile
from subprocess import call
from signalalign.signalAlignment import *
from signalalign.fast5 import Fast5
from signalalign import parseFofn
from signalalign.filter_reads import filter_reads, filter_reads_to_string_wrapper
from signalalign.utils.fileHandlers import FolderHandler
from py3helpers.utils import captured_output, merge_dicts, list_dir, create_dot_dict, load_json, save_json
from py3helpers.seq_tools import ReverseComplement


class SignalAlignmentTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(SignalAlignmentTest, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.reference = os.path.join(cls.HOME, "tests/test_sequences/pUC19_SspI_Zymo.fa")
        cls.ecoli_reference = os.path.join(cls.HOME, "tests/test_sequences/E.coli_K12.fasta")

        cls.fast5_dir = os.path.join(cls.HOME, "tests/minion_test_reads/canonical_ecoli_R9")
        cls.files = [
            "miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_strand.fast5",
            "miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch101_read456_strand.fast5",
            "miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch101_read544_strand1.fast5",
            "miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch103_read333_strand1.fast5"]
        cls.fast5_paths = list_dir(cls.fast5_dir, ext="fast5")
        cls.fast5_bam = os.path.join(cls.HOME, "tests/minion_test_reads/canonical_ecoli_R9/canonical_ecoli.bam")
        cls.fast5_readdb = os.path.join(cls.HOME, "tests/minion_test_reads/canonical_ecoli_R9/canonical_ecoli.readdb")

        cls.template_hmm = os.path.join(cls.HOME, "models/testModelR9_acgt_template.model")
        cls.path_to_bin = os.path.join(cls.HOME, 'bin')
        cls.tmp_directory = tempfile.mkdtemp()
        cls.test_dir = os.path.join(cls.tmp_directory, "test")

        dna_dir = os.path.join(cls.HOME, "tests/minion_test_reads/1D/")
        # copy file to tmp directory
        shutil.copytree(dna_dir, cls.test_dir)
        cls.readdb = os.path.join(cls.HOME, "tests/minion_test_reads/oneD.fastq.index.readdb")
        cls.bam = os.path.join(cls.HOME, "tests/minion_test_reads/oneD.bam")

        cls.rna_bam = os.path.join(cls.HOME, "tests/minion_test_reads/RNA_edge_cases/rna_reads.bam")
        cls.rna_readdb = os.path.join(cls.HOME, "tests/minion_test_reads/RNA_edge_cases/rna_reads.readdb")
        cls.test_dir_rna = os.path.join(cls.tmp_directory, "test_rna")
        cls.rna_reference = os.path.join(cls.HOME, "tests/test_sequences/fake_rna_ref.fa")

        rna_dir = os.path.join(cls.HOME, "tests/minion_test_reads/RNA_edge_cases/")
        # copy file to tmp directory
        shutil.copytree(rna_dir, cls.test_dir_rna)
        # used to test runSignalAlign with config file
        cls.config_file = os.path.join(cls.HOME, "tests/runSignalAlign-config.json")
        cls.default_args = create_dot_dict(load_json(cls.config_file))

    def test_create_signalAlignment_args(self):
        expected_args = {"backward_reference", "forward_reference", "destination", "stateMachineType", "in_templateHmm",
                         "in_complementHmm", "in_templateHdp", "in_complementHdp", "threshold", "diagonal_expansion",
                         "constraint_trim", "target_regions", "twoD_chemistry", "alignment_file",
                         "bwa_reference",
                         'track_memory_usage', 'get_expectations', 'output_format', 'embed', 'event_table',
                         'check_for_temp_file_existance', 'path_to_bin', 'perform_kmer_event_alignment', 'filter_reads',
                         'traceBackDiagonals', 'delete_tmp'}
        args = create_signalAlignment_args()
        self.assertSetEqual(set(args.keys()), expected_args)

    def test_create_sa_sample_args(self):
        expected_args = {"fofns", "fast5_dirs", "positions_file", "motifs", "bwa_reference", "fw_reference",
                         "bw_reference", "name", "number_of_kmer_assignments", "probability_threshold",
                         "kmers_from_reference", 'alignment_file', "quality_threshold", "recursive",
                         "workers", "assignments_dir", "readdb"}
        args = create_sa_sample_args()
        self.assertSetEqual(set(args.keys()), expected_args)

    def test_parseFofn(self):
        with tempfile.TemporaryDirectory() as tempdir:
            test_out = os.path.join(tempdir, "test.fofn")
            with open(test_out, 'w+') as fofn_file:
                for path in self.fast5_paths:
                    print(path, file=fofn_file)

            files = parseFofn(test_out)
            self.assertSequenceEqual(files, self.fast5_paths)
            with open(test_out, 'w+') as fofn_file:
                for x in range(10):
                    print("Something else", file=fofn_file)
            self.assertRaises(AssertionError, parseFofn, test_out)
            self.assertRaises(AssertionError, parseFofn, "fake_file")

    def test_signal_file_and_alignment(self):
        signal_file_reads = os.path.join(self.HOME, "tests/minion_test_reads/no_event_data_1D_ecoli")
        template_model = os.path.join(self.HOME, "models/testModelR9p4_5mer_acegt_template.model")
        ecoli_reference = os.path.join(self.HOME, "tests/test_sequences/E.coli_K12.fasta")
        signal_file_guide_alignment = os.path.join(self.HOME, "tests/minion_test_reads/oneD_alignments.sam")

        with tempfile.TemporaryDirectory() as tempdir:
            new_dir = os.path.join(tempdir, "new_dir")
            working_folder = FolderHandler()
            working_folder.open_folder(os.path.join(tempdir, "test_dir"))

            shutil.copytree(signal_file_reads, new_dir)

            args = create_signalAlignment_args(alignment_file=signal_file_guide_alignment, bwa_reference=ecoli_reference,
                                               forward_reference=ecoli_reference, in_templateHmm=template_model,
                                               path_to_bin=self.path_to_bin, destination=working_folder.path)
            final_args = merge_dicts([args, dict(in_fast5=os.path.join(new_dir, "LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch6_read347_strand.fast5"))])
            handle = SignalAlignment(**final_args)
            handle.run()
            self.assertEqual(len(os.listdir(working_folder.path)), 1)
            self.assertEqual(sorted(os.listdir(working_folder.path))[0], "9e4d14b1-8167-44ef-9fdb-5c29dd0763fd.sm.backward.tsv")

    def test_embed(self):
        signal_file_reads = os.path.join(self.HOME, "tests/minion_test_reads/no_event_data_1D_ecoli")
        template_model = os.path.join(self.HOME, "models/testModelR9p4_5mer_acegt_template.model")
        ecoli_reference = os.path.join(self.HOME, "tests/test_sequences/E.coli_K12.fasta")
        signal_file_guide_alignment = os.path.join(self.HOME, "tests/minion_test_reads/oneD_alignments.sam")

        with tempfile.TemporaryDirectory() as tempdir:
            new_dir = os.path.join(tempdir, "new_dir")
            working_folder = FolderHandler()
            working_folder.open_folder(os.path.join(tempdir, "test_dir"))

            shutil.copytree(signal_file_reads, new_dir)

            args = create_signalAlignment_args(alignment_file=signal_file_guide_alignment, bwa_reference=ecoli_reference,
                                               forward_reference=ecoli_reference, in_templateHmm=template_model,
                                               path_to_bin=self.path_to_bin, destination=working_folder.path,
                                               embed=True)
            final_args = merge_dicts([args, dict(in_fast5=os.path.join(new_dir, "LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch6_read347_strand.fast5"))])
            handle = SignalAlignment(**final_args)
            handle.run()
            f5fh = Fast5(os.path.join(new_dir, "LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch6_read347_strand.fast5"))
            mea = f5fh.get_signalalign_events(mea=True)
            sam = f5fh.get_signalalign_events(sam=True)
            self.assertEqual(mea[0]["raw_start"], 153)
            self.assertEqual(sam[0], "9")
            self.assertEqual(len(os.listdir(working_folder.path)), 1)
            self.assertEqual(sorted(os.listdir(working_folder.path))[0], "9e4d14b1-8167-44ef-9fdb-5c29dd0763fd.sm.backward.tsv")

        # DNA WITH events
        signal_file_reads = os.path.join(self.HOME, "tests/minion_test_reads/1D")
        template_model = os.path.join(self.HOME, "models/testModelR9p4_5mer_acegt_template.model")
        ecoli_reference = os.path.join(self.HOME, "tests/test_sequences/E.coli_K12.fasta")
        signal_file_guide_alignment = os.path.join(self.HOME, "tests/minion_test_reads/oneD_alignments.sam")

        with tempfile.TemporaryDirectory() as tempdir:
            new_dir = os.path.join(tempdir, "new_dir")
            working_folder = FolderHandler()
            working_folder.open_folder(os.path.join(tempdir, "test_dir"))

            shutil.copytree(signal_file_reads, new_dir)

            args = create_signalAlignment_args(alignment_file=signal_file_guide_alignment, bwa_reference=ecoli_reference,
                                               forward_reference=ecoli_reference, in_templateHmm=template_model,
                                               path_to_bin=self.path_to_bin, destination=working_folder.path,
                                               embed=True)
            final_args = merge_dicts([args, dict(in_fast5=os.path.join(new_dir, "LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch6_read347_strand.fast5"))])
            handle = SignalAlignment(**final_args)
            handle.run()
            f5fh = Fast5(os.path.join(new_dir, "LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch6_read347_strand.fast5"))
            mea = f5fh.get_signalalign_events(mea=True)
            sam = f5fh.get_signalalign_events(sam=True)
            self.assertEqual(mea[0]["raw_start"], 153)
            self.assertEqual(sam[0], "9")
            self.assertEqual(len(os.listdir(working_folder.path)), 1)
            self.assertEqual(sorted(os.listdir(working_folder.path))[0], "9e4d14b1-8167-44ef-9fdb-5c29dd0763fd.sm.backward.tsv")

    def test_multithread_signal_alignment(self):
        with tempfile.TemporaryDirectory() as tempdir:
            working_folder = FolderHandler()
            working_folder.open_folder(os.path.join(tempdir, "test_dir"))
            # create signalalign args
            assert os.path.isfile(self.template_hmm)
            signal_align_arguments = create_signalAlignment_args(bwa_reference=self.ecoli_reference,
                                                                 in_templateHmm=self.template_hmm,
                                                                 destination=working_folder.path,
                                                                 forward_reference=self.ecoli_reference,
                                                                 path_to_bin=self.path_to_bin)

            fast5_files = self.fast5_paths[:1]
            with captured_output() as (out, err):
                output_files = multithread_signal_alignment(signal_align_arguments, fast5_files, 2,
                                                            forward_reference=self.ecoli_reference)
            self.assertEqual(len(output_files), len(fast5_files))
            # round 2
            working_folder = FolderHandler()
            working_folder.open_folder(os.path.join(tempdir, "test_dir2"))
            # create signalalign args
            assert os.path.isfile(self.template_hmm)
            signal_align_arguments = create_signalAlignment_args(bwa_reference=self.ecoli_reference,
                                                                 in_templateHmm=self.template_hmm,
                                                                 destination=working_folder.path,
                                                                 forward_reference=self.ecoli_reference,
                                                                 path_to_bin=self.path_to_bin)

            filter_reads_generator = filter_reads_to_string_wrapper(filter_reads(self.bam, self.readdb, [self.test_dir]))
            output_files = multithread_signal_alignment(signal_align_arguments, [], 2,
                                                        forward_reference=self.ecoli_reference,
                                                        filter_reads_to_string_wrapper=filter_reads_generator,
                                                        debug=True)
            self.assertEqual(len(output_files), 3)

    def test_multithread_signal_alignment_samples(self):
        with tempfile.TemporaryDirectory() as tempdir:
            working_folder = FolderHandler()
            test_fast5 = os.path.join(tempdir, "miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_strand.fast5")
            num_files = 1
            copyfile(self.fast5_paths[0], test_fast5)
            working_folder.open_folder(os.path.join(tempdir, "test_dir"))
            # create signalalign args
            signal_align_arguments = create_signalAlignment_args(in_templateHmm=self.template_hmm,
                                                                 destination=working_folder.path,
                                                                 path_to_bin=self.path_to_bin)
            # create samples
            samples = []
            options = create_sa_sample_args(fast5_dirs=[tempdir], name="some_name",
                                            fw_reference=self.ecoli_reference,
                                            bwa_reference=self.ecoli_reference,
                                            readdb=self.fast5_readdb,
                                            alignment_file=self.fast5_bam)
            samples.append(SignalAlignSample(working_folder=working_folder, **options))
            options["name"] = "some_name2"
            samples.append(SignalAlignSample(working_folder=working_folder, **options))
            # with captured_output() as (out, err):
            samples = multithread_signal_alignment_samples(samples, signal_align_arguments, 2)
            self.assertSetEqual(set([sample.name for sample in samples]), {'some_name', 'some_name2'})
            for sample in samples:
                if sample.name == "some_name":
                    self.assertEqual(len(sample.analysis_files), num_files)
                if sample.name == "some_name2":
                    self.assertEqual(len(sample.analysis_files), num_files)
                for file_path in sample.analysis_files:
                    self.assertTrue(os.path.isfile(file_path))
            options["name"] = "some_name"
            samples.append(SignalAlignSample(working_folder=working_folder, **options))
            self.assertRaises(AssertionError, multithread_signal_alignment_samples,samples, signal_align_arguments, 2)

    def test_SignalAlignSample(self):
        with tempfile.TemporaryDirectory() as tempdir:
            # create fast5 dir
            test_fast5 = os.path.join(tempdir, "test.fast5")
            copyfile(self.fast5_paths[0], test_fast5)
            # create fofn
            test_out = os.path.join(tempdir, "test.fofn")
            with open(test_out, 'w+') as fofn_file:
                print(test_fast5, file=fofn_file)

            test_args = create_sa_sample_args(fast5_dirs=[tempdir, tempdir],
                                              name="some_name",
                                              fofns=[test_out, test_out],
                                              fw_reference=self.ecoli_reference,
                                              bwa_reference=self.ecoli_reference)

            working_folder = FolderHandler()
            working_folder.open_folder(os.path.join(tempdir, "test_dir"))
            sample = SignalAlignSample(working_folder=working_folder, **test_args)

    def test_rna_reads(self):
        with tempfile.TemporaryDirectory() as tempdir:
            template_model = os.path.join(self.HOME, "models/testModelR9p4_5mer_acgt_RNA.model")
            args = create_signalAlignment_args(alignment_file=self.rna_bam, bwa_reference=self.rna_reference,
                                               forward_reference=self.rna_reference, in_templateHmm=template_model,
                                               path_to_bin=self.path_to_bin, destination=tempdir, embed=True,
                                               delete_tmp=False)

            in_rna_file = os.path.join(self.test_dir_rna, "DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_36_ch_218_strand.fast5")
            final_args = merge_dicts([args, dict(in_fast5=in_rna_file)])
            handle = SignalAlignment(**final_args)
            handle.run()
            fh = pysam.FastaFile(self.rna_reference)
            f5fh = Fast5(in_rna_file)
            sa_events = f5fh.get_signalalign_events()
            for i, event in enumerate(sa_events):
                kmer = fh.fetch(reference="rna_fake", start=event["reference_index"], end=event["reference_index"]+5)[::-1]
                self.assertEqual(event["path_kmer"].decode(), kmer)
                self.assertEqual(event["reference_kmer"].decode(), kmer)

            in_rna_file = os.path.join(self.test_dir_rna, "DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5")
            final_args = merge_dicts([args, dict(in_fast5=in_rna_file)])
            handle = SignalAlignment(**final_args)
            handle.run()
            rev_c = ReverseComplement()
            f5fh = Fast5(in_rna_file)
            sa_events = f5fh.get_signalalign_events()
            for i, event in enumerate(sa_events):
                kmer = fh.fetch(reference="rna_fake", start=event["reference_index"], end=event["reference_index"]+5)[::-1]
                rev_kmer = rev_c.reverse_complement(kmer)
                self.assertEqual(event["path_kmer"].decode(), rev_kmer)
                self.assertEqual(event["reference_kmer"].decode(), kmer)

    def test_read_in_signal_align_tsv(self):
        example = os.path.join(self.HOME, "tests/test_assignment_files/d6160b0b-a35e-43b5-947f-adaa1abade28.sm.assignments.tsv")
        data = SignalAlignment.read_in_signal_align_tsv(example, "assignments")
        self.assertEqual(NanoporeRead.bytes_to_string(data["k-mer"][0]), "GCCTTA")
        with tempfile.TemporaryDirectory() as tempdir:
            file1 = os.path.join(tempdir, "test.txt")
            with open(file1, "w") as f:
                subprocess.call(["head", "-1", example], stdout=f)
            data = SignalAlignment.read_in_signal_align_tsv(file1, "assignments")
            self.assertEqual(NanoporeRead.bytes_to_string(data["k-mer"][0]), "GCCTTA")

    def test_embed_with_both(self):
        signal_file_reads = os.path.join(self.HOME, "tests/minion_test_reads/pUC/")
        template_model = os.path.join(self.HOME, "models/testModelR9_5mer_acegt_template.model")
        complement_model = os.path.join(self.HOME, "models/testModelR9_5mer_acegt_complement.model")

        puc_reference = os.path.join(self.HOME, "tests/test_sequences/pUC19_SspI.fa")
        signal_file_guide_alignment = os.path.join(self.HOME, "tests/minion_test_reads/pUC/puc.bam")
        with tempfile.TemporaryDirectory() as tempdir:
            new_dir = os.path.join(tempdir, "new_dir")
            if os.path.exists(new_dir):
                shutil.rmtree(new_dir)
            working_folder = FolderHandler()
            working_folder.open_folder(os.path.join(tempdir, "test_dir"))

            shutil.copytree(signal_file_reads, new_dir)

            args = create_signalAlignment_args(alignment_file=signal_file_guide_alignment, bwa_reference=puc_reference,
                                               forward_reference=puc_reference, in_templateHmm=template_model,
                                               path_to_bin=self.path_to_bin, destination=working_folder.path,
                                               embed=True, output_format="both", filter_reads=0, twoD_chemistry=True,
                                               in_complementHmm=complement_model, delete_tmp=True)
            final_args = merge_dicts([args, dict(in_fast5=os.path.join(new_dir, "makeson_PC_20160807_FNFAD20242_MN17284_sequencing_run_MA_470_R9_pUC_g_PCR_BC_08_07_16_93165_ch1_read176_strand.fast5"))])
            handle = SignalAlignment(**final_args)
            handle.run()
            f5fh = Fast5(os.path.join(new_dir, "makeson_PC_20160807_FNFAD20242_MN17284_sequencing_run_MA_470_R9_pUC_g_PCR_BC_08_07_16_93165_ch1_read176_strand.fast5"))
            mea = f5fh.get_signalalign_events(mea=True)
            sam = f5fh.get_signalalign_events(sam=True)
            self.assertEqual(mea[0]["raw_start"], 2879)
            self.assertEqual(sam[0], "0")
            self.assertEqual(len(os.listdir(working_folder.path)), 2)

    def test_multithread_signal_alignment_samples(self):
        with tempfile.TemporaryDirectory() as tempdir:
            # tempdir = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/1testing"
            working_folder = FolderHandler()
            test_fast5 = os.path.join(tempdir, "miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_strand.fast5")
            num_files = 1
            copyfile(self.fast5_paths[0], test_fast5)
            working_folder.open_folder(os.path.join(tempdir, "test_dir"))
            # create signalalign args
            brdu_template_model = os.path.join(self.HOME, "models/BrdU_sa_model_threshold2.0.model")

            signal_align_arguments = create_signalAlignment_args(in_templateHmm=brdu_template_model,
                                                                 destination=working_folder.path,
                                                                 path_to_bin=self.path_to_bin,
                                                                 delete_tmp=False)
            # create samples
            samples = []
            options = create_sa_sample_args(fast5_dirs=[tempdir], name="some_name",
                                            motifs=[["T", "Z"]],
                                            # motifs=[["GATAAT", "GAXAAT"], ["GATAAT", "GATAAX"]],
                                            bwa_reference=self.ecoli_reference,
                                            readdb=self.fast5_readdb,
                                            alignment_file=self.fast5_bam)
            samples.append(SignalAlignSample(working_folder=working_folder, **options))
            # with captured_output() as (out, err):
            samples = multithread_signal_alignment_samples(samples, signal_align_arguments, 2)
            for sample in samples:
                if sample.name == "some_name":
                    self.assertEqual(len(sample.analysis_files), num_files)
                for file_path in sample.analysis_files:
                    self.assertTrue(os.path.isfile(file_path))

    def test_variant_calling_with_multiple_paths_rna(self):
        with tempfile.TemporaryDirectory() as tempdir:
            new_dir = os.path.join(tempdir, "new_dir")
            if os.path.exists(new_dir):
                shutil.rmtree(new_dir)
            working_folder = FolderHandler()
            working_folder.open_folder(os.path.join(tempdir, "test_dir"))

            shutil.copytree(self.test_dir_rna, new_dir)

            args = create_signalAlignment_args(alignment_file=self.rna_bam, bwa_reference=self.rna_reference,
                                               forward_reference=os.path.join(self.HOME, "tests/test_sequences/fake_rna_replace/forward.fake_rna_atg.fake_rna_ref.fa"),
                                               backward_reference=os.path.join(self.HOME, "tests/test_sequences/fake_rna_replace/backward.fake_rna_atg.fake_rna_ref.fa"),
                                               in_templateHmm=os.path.join(self.HOME, "models/fake_testModelR9p4_5mer_acfgt_RNA.model"),
                                               path_to_bin=self.path_to_bin, destination=working_folder.path,
                                               embed=False, output_format="full", filter_reads=0, twoD_chemistry=False,
                                               delete_tmp=True, check_for_temp_file_existance=False)

            multithread_signal_alignment(args, list_dir(new_dir, ext="fast5"), worker_count=8,
                                         forward_reference=None,
                                         debug=True, filter_reads_to_string_wrapper=None)
            self.assertEqual(len(os.listdir(working_folder.path)), 2)


if __name__ == '__main__':
    unittest.main()
