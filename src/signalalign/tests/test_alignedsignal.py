#!/usr/bin/env python
"""Tests for alignedsignal.py"""
########################################################################
# File: test_alignedsignal.py
#  executable: test_alignedsignal.py
#
# Author: Andrew Bailey
# History: Created 03/09/18
########################################################################

import unittest
import os
import numpy as np
import tempfile
import shutil
import pysam

from signalalign.alignedsignal import *
from signalalign.visualization.plot_labelled_read import PlotSignal
from signalalign.signalAlignment import SignalAlignment, create_signalAlignment_args
from py3helpers.utils import merge_dicts, binary_search
from py3helpers.seq_tools import ReverseComplement, sam_string_to_aligned_segment


class CreateLabelsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(CreateLabelsTest, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.fasta = os.path.join(cls.HOME,
                                 "tests/test_sequences/E.coli_K12.fasta")
        dna_file = os.path.join(cls.HOME,
                                "tests/minion_test_reads/1D/LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch112_read108_strand.fast5")
        rev_dna_file = os.path.join(cls.HOME,
                                    "tests/minion_test_reads/1D/LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch6_read347_strand.fast5")
        rev_rna_file = os.path.join(cls.HOME,
                                "tests/minion_test_reads/RNA_no_events/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5")
        forward_rna_file = os.path.join(cls.HOME,
                                "tests/minion_test_reads/RNA_no_events/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_36_ch_218_strand.fast5")

        rna_reference = os.path.join(cls.HOME, "tests/test_sequences/fake_rna_ref.fa")
        ecoli_dna_reference = os.path.join(cls.HOME, "tests/test_sequences/E.coli_K12.fasta")
        cls.dna_reference_handle = pysam.FastaFile(ecoli_dna_reference)
        cls.rna_reference_handle = pysam.FastaFile(rna_reference)
        cls.tmp_directory = tempfile.mkdtemp()

         # get file locations
        cls.tmp_dna_file = os.path.join(str(cls.tmp_directory), 'test_dna.fast5')
        cls.tmp_dna_file2 = os.path.join(str(cls.tmp_directory), 'test_dna2.fast5')

        cls.tmp_rna_file1 = os.path.join(str(cls.tmp_directory), 'test_rna.fast5')
        cls.tmp_rna_file2 = os.path.join(str(cls.tmp_directory), 'test_rna2.fast5')

        # run signalAlign on one file
        cls.rna_model_file = os.path.join(cls.HOME, "models/testModelR9p4_5mer_acgt_RNA.model")
        cls.dna_model_file_94 = os.path.join(cls.HOME, "models/testModelR9p4_5mer_acegt_template.model")
        cls.rna_sam = os.path.join(cls.HOME, "tests/minion_test_reads/RNA_edge_cases/rna_reads.bam")
        cls.dna_sam = os.path.join(cls.HOME, "tests/minion_test_reads/oneD.bam")
        cls.bin_path = os.path.join(cls.HOME, "bin")
        # kmer index
        cls.kmer_index = 2

        # copy file to tmp directory
        shutil.copy(dna_file, cls.tmp_dna_file)
        shutil.copy(rev_dna_file, cls.tmp_dna_file2)

        shutil.copy(forward_rna_file, cls.tmp_rna_file1)
        shutil.copy(rev_rna_file, cls.tmp_rna_file2)

        args = create_signalAlignment_args(destination=cls.tmp_directory,
                                           in_templateHmm=cls.rna_model_file,
                                           alignment_file=cls.rna_sam,
                                           forward_reference=rna_reference,
                                           embed=True,
                                           path_to_bin=cls.bin_path,
                                           diagonal_expansion=5,
                                           delete_tmp=False)
        sa_h = SignalAlignment(**merge_dicts([args, {'in_fast5': cls.tmp_rna_file1}]))
        sa_h.run()

        sa_h = SignalAlignment(**merge_dicts([args, {'in_fast5': cls.tmp_rna_file2}]))
        sa_h.run()

        args = create_signalAlignment_args(destination=cls.tmp_directory,
                                           in_templateHmm=cls.dna_model_file_94,
                                           alignment_file=cls.dna_sam,
                                           forward_reference=ecoli_dna_reference,
                                           embed=True,
                                           path_to_bin=cls.bin_path,
                                           diagonal_expansion=10,
                                           traceBackDiagonals=100,
                                           constraint_trim=3)
        sa_h = SignalAlignment(**merge_dicts([args, {'in_fast5': cls.tmp_dna_file}]))
        sa_h.run()

        sa_h = SignalAlignment(**merge_dicts([args, {'in_fast5': cls.tmp_dna_file2}]))
        sa_h.run()

        cls.dna_handle = CreateLabels(cls.tmp_dna_file, kmer_index=cls.kmer_index)
        cls.dna_handle2 = CreateLabels(cls.tmp_dna_file2, kmer_index=cls.kmer_index)

        cls.rna1_handle = CreateLabels(cls.tmp_rna_file1, kmer_index=cls.kmer_index)
        cls.rna2_handle = CreateLabels(cls.tmp_rna_file2, kmer_index=cls.kmer_index)
        cls.rev_comp = ReverseComplement()

        cls.tmp_dna_file3 = os.path.join(cls.HOME,
                                         "tests/minion_test_reads/embedded_files/miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch100_read2324_strand.fast5")
        cls.dna3_handle = CreateLabels(cls.tmp_dna_file3, kmer_index=cls.kmer_index)


    def test_initialize(self):
        self.assertEqual(self.dna_handle.aligned_signal.raw_signal[0], 1172)
        self.assertEqual(self.dna_handle2.aligned_signal.raw_signal[0], 964)

        self.assertEqual(self.rna1_handle.aligned_signal.raw_signal[0], 682)
        self.assertEqual(self.rna2_handle.aligned_signal.raw_signal[0], 647)

    def test_fix_sa_reference_indexes_dna_forward(self):
        mea_alignment1 = self.dna_handle.get_signalalign_events(mea=True)
        # test forward alignment of dna (ref index refers to start of 5'-3' kmer)
        self.dna_handle2.aligned_signal.minus_strand = False
        for row in mea_alignment1:
            kmer = row["kmer"].decode()
            base = kmer[self.kmer_index]

            ref_kmer = self.dna_reference_handle.fetch(reference="gi_ecoli", start=row["reference_index"],
                                                       end=row["reference_index"]+5)
            self.assertEqual(kmer, ref_kmer)
        # changes in place
        self.dna_handle.fix_sa_reference_indexes(mea_alignment1)
        for row in mea_alignment1:
            kmer = row["kmer"].decode()
            base = kmer[self.kmer_index]
            ref_base = self.dna_reference_handle.fetch(reference="gi_ecoli", start=row["reference_index"],
                                                       end=row["reference_index"]+1)
            self.assertEqual(base, ref_base)

    def test_fix_sa_reference_indexes_dna_reverse(self):
        mea_alignment2 = self.dna_handle2.get_signalalign_events(mea=True)
        # test reverse strand alignment of dna (ref index refers to start of 5'-3' kmer)
        self.dna_handle2.aligned_signal.minus_strand = True
        for row in mea_alignment2:
            kmer = row["kmer"].decode()
            ref_kmer = self.dna_reference_handle.fetch(reference="gi_ecoli", start=row["reference_index"],
                                                       end=row["reference_index"]+5)
            minus_strand_kmer = self.rev_comp.reverse_complement(ref_kmer)
            self.assertEqual(kmer, minus_strand_kmer)
        # changes in place
        self.dna_handle2.fix_sa_reference_indexes(mea_alignment2)
        for row in mea_alignment2:
            kmer = row["kmer"].decode()
            base = kmer[self.kmer_index]
            ref_base = self.dna_reference_handle.fetch(reference="gi_ecoli", start=row["reference_index"],
                                                       end=row["reference_index"]+1)
            minus_strand_base = self.rev_comp.reverse_complement(ref_base)

            self.assertEqual(base, minus_strand_base)

    def test_fix_sa_reference_indexes_rna_forward(self):
        mea_alignment1 = self.rna1_handle.get_signalalign_events(mea=True)
        # test forward strand alignment of rna (ref index refers to end of 3'-5' kmer)
        self.rna1_handle.aligned_signal.minus_strand = False
        for row in mea_alignment1:
            kmer = row["kmer"].decode()
            ref_kmer = self.rna_reference_handle.fetch(reference="rna_fake", start=row["reference_index"],
                                                       end=row["reference_index"]+5)
            # reverse kmer because it is 3'-5' for RNA
            self.assertEqual(kmer[::-1], ref_kmer)
        # changes in place
        self.rna1_handle.fix_sa_reference_indexes(mea_alignment1)
        for row in mea_alignment1:
            kmer = row["kmer"].decode()
            base = kmer[self.kmer_index]
            ref_base = self.rna_reference_handle.fetch(reference="rna_fake", start=row["reference_index"],
                                                       end=row["reference_index"]+1)

            self.assertEqual(base, ref_base)

    def test_fix_sa_reference_indexes_rna_reversed(self):
        mea_alignment2 = self.rna2_handle.get_signalalign_events(mea=True)
        # test reverse strand alignment of rna (ref index refers to end of 3'-5' kmer)
        self.rna2_handle.aligned_signal.minus_strand = True
        for row in mea_alignment2:
            kmer = row["kmer"].decode()
            ref_kmer = self.rna_reference_handle.fetch(reference="rna_fake", start=row["reference_index"],
                                                       end=row["reference_index"]+5)
            minus_strand_kmer = self.rev_comp.reverse_complement(ref_kmer)

            # reverse kmer because it is 3'-5' for RNA
            self.assertEqual(kmer[::-1], minus_strand_kmer)
        # changes in place
        self.rna2_handle.fix_sa_reference_indexes(mea_alignment2)
        for row in mea_alignment2:
            kmer = row["kmer"].decode()
            base = kmer[self.kmer_index]
            ref_base = self.rna_reference_handle.fetch(reference="rna_fake", start=row["reference_index"],
                                                       end=row["reference_index"]+1)
            minus_strand_base = self.rev_comp.reverse_complement(ref_base)
            self.assertEqual(base, minus_strand_base)

    def test_add_signal_align_predictions(self):
        self.dna_handle.add_signal_align_predictions(add_basecall=False)
        self.assertSequenceEqual(self.dna_handle.aligned_signal.prediction["full_signalalign"][0].tolist(),
                                 [774, 3, 3560628+self.kmer_index, 1., b'CGTTT'])

        self.dna_handle2.add_signal_align_predictions(add_basecall=False)
        self.assertSequenceEqual(self.dna_handle2.aligned_signal.prediction["full_signalalign"][0].tolist(),
                                 [153, 3, 1845108+self.kmer_index, 0.992924, b'CATTG'])

        self.rna1_handle.add_signal_align_predictions(add_basecall=False)
        self.assertSequenceEqual(self.rna1_handle.aligned_signal.prediction["full_signalalign"][0].tolist(),
                                 [0, 11, 1081+self.kmer_index, 0.999781, b'AACCT'])
        self.rna2_handle.add_signal_align_predictions(add_basecall=False)
        self.assertSequenceEqual(self.rna2_handle.aligned_signal.prediction["full_signalalign"][0].tolist(),
                                 [0, 7, 3+self.kmer_index, 1.0, b'AACCT'])

    def test_add_mea_labels(self):
        """Test add mea labels"""
        self.dna_handle.add_mea_labels()
        self.assertSequenceEqual(self.dna_handle.aligned_signal.label["mea_signalalign"][0].tolist(),
                                 [774, 3, 3560628+self.kmer_index, 1., b'CGTTT'])

        self.dna_handle2.add_mea_labels()
        self.assertSequenceEqual(self.dna_handle2.aligned_signal.label["mea_signalalign"][0].tolist(),
                                 [153, 3, 1845108+self.kmer_index, 0.992924, b'CATTG'])

        self.rna1_handle.add_mea_labels()
        self.assertSequenceEqual(self.rna1_handle.aligned_signal.label["mea_signalalign"][0].tolist(),
                                 [0, 11, 1081+self.kmer_index, 0.999781, b'AACCT'])
        self.rna2_handle.add_mea_labels()
        self.assertSequenceEqual(self.rna2_handle.aligned_signal.label["mea_signalalign"][0].tolist(),
                                 [0, 7, 3+self.kmer_index, 1.0, b'AACCT'])

    def test_add_basecall_alignment_prediction(self):
        self.rna1_handle.add_basecall_alignment_prediction()
        self.assertSequenceEqual(self.rna1_handle.aligned_signal.prediction["matches_guide_alignment"][0].tolist(),
                                 [0, 11, 1081+self.kmer_index, 0.025408448110250344, b'C'])
        self.rna2_handle.add_basecall_alignment_prediction()
        self.assertSequenceEqual(self.rna2_handle.aligned_signal.prediction["matches_guide_alignment"][0].tolist(),
                                 [0, 7, 3+self.kmer_index, 0.026885957628966044, b'C'])

        self.dna_handle.add_basecall_alignment_prediction()
        self.assertSequenceEqual(self.dna_handle.aligned_signal.prediction["matches_guide_alignment"][0].tolist(),
                                 [774, 3, 3560628+self.kmer_index, 0.8552098274230957, b'T'])

        self.dna_handle2.add_basecall_alignment_prediction()
        self.assertSequenceEqual(self.dna_handle2.aligned_signal.prediction["matches_guide_alignment"][0].tolist(),
                                 [153, 3, 1845108+self.kmer_index, 0.11882638931274414, b'T'])

    def test_match_cigar_with_basecall_guide_rna1(self):
        events = self.rna1_handle.get_basecall_data()
        sam = self.rna1_handle.get_signalalign_events(sam=True)
        rna = self.rna1_handle.rna
        matches, mismatches, raw_start = match_cigar_with_basecall_guide(events, sam, self.kmer_index,
                                                                         rna=rna, one_ref_indexing=False)
        for match in matches:
            ref_base = self.rna_reference_handle.fetch(reference="rna_fake", start=match["reference_index"],
                                                       end=match["reference_index"]+1)
            self.assertEqual(match["kmer"].decode(), ref_base)

    def test_match_cigar_with_basecall_guide_rna2(self):
        events = self.rna2_handle.get_basecall_data()
        sam = self.rna2_handle.get_signalalign_events(sam=True)
        rna = self.rna2_handle.rna
        matches, mismatches, raw_start = match_cigar_with_basecall_guide(events, sam, self.kmer_index,
                                                                         rna=rna, one_ref_indexing=False)
        for match in matches:
            ref_base = self.rna_reference_handle.fetch(reference="rna_fake", start=match["reference_index"],
                                                       end=match["reference_index"]+1)
            self.assertEqual(match["kmer"].decode(), self.rev_comp.complement(ref_base))

    def test_match_cigar_with_basecall_guide_dna1(self):
        events = self.dna_handle.get_basecall_data()
        sam = self.dna_handle.get_signalalign_events(sam=True)
        events = add_raw_start_and_raw_length_to_events(events, self.dna_handle.sample_rate,
                                                        self.dna_handle.raw_attributes["start_time"])

        rna = self.dna_handle.rna
        matches, mismatches, raw_start = match_cigar_with_basecall_guide(events, sam, self.kmer_index,
                                                                         rna=rna, one_ref_indexing=False)
        for match in matches:
            ref_base = self.dna_reference_handle.fetch(reference="gi_ecoli", start=match["reference_index"],
                                                       end=match["reference_index"]+1)
            self.assertEqual(match["kmer"].decode(), ref_base)

    def test_match_cigar_with_basecall_guide_dna2(self):
        events = self.dna_handle2.get_basecall_data()
        sam = self.dna_handle2.get_signalalign_events(sam=True)
        rna = self.dna_handle2.rna
        events = add_raw_start_and_raw_length_to_events(events, self.dna_handle2.sample_rate,
                                                        self.dna_handle2.raw_attributes["start_time"])

        matches, mismatches, raw_start = match_cigar_with_basecall_guide(events, sam, self.kmer_index,
                                                                         rna=rna, one_ref_indexing=False)
        for match in matches:
            ref_base = self.dna_reference_handle.fetch(reference="gi_ecoli", start=match["reference_index"],
                                                       end=match["reference_index"]+1)
            self.assertEqual(match["kmer"].decode(), self.rev_comp.complement(ref_base))

    def test_index_bases_from_events(self):
        events = self.rna1_handle.get_basecall_data()
        bases, base_raw_starts, base_raw_lengths, probs = index_bases_from_events(events, kmer_index=self.kmer_index)
        first_event = events[0]
        self.assertEqual(bases[self.kmer_index], first_event["model_state"].decode()[self.kmer_index])
        self.assertEqual(base_raw_starts[self.kmer_index], first_event["raw_start"])
        self.assertEqual(base_raw_lengths[self.kmer_index], first_event["raw_length"])
        self.assertEqual(probs[self.kmer_index], first_event["p_model_state"])

    def test_index_bases_from_events2(self):
        """Test index_bases_from_events"""
        # make sure each event is corresponding to correct nucleotide
        events = np.zeros(4, dtype=[('raw_start', int), ('raw_length', int), ('move', int),
                                    ('p_model_state', float), ('model_state', 'S5')])
        events["raw_start"] = [0, 1, 2, 3]
        events["raw_length"] = [1, 1, 1, 1]
        events["move"] = [1, 1, 1, 1]
        events["p_model_state"] = [1, 1, 1, 1]
        events["model_state"] = ["GATTA", "ATTAC", "TTACA", "TACAG"]

        bases, base_raw_starts, base_raw_lengths, probs = index_bases_from_events(events, kmer_index=2)
        self.assertSequenceEqual(bases, list("GATTACAG"))
        self.assertSequenceEqual(base_raw_lengths, [1, 1, 1, 1, 1, 1, 1, 1])
        self.assertSequenceEqual(probs, [1, 1, 1, 1, 1, 1, 1, 1])
        self.assertSequenceEqual(base_raw_starts, [0, 0, 0, 1, 2, 3, 3, 3])
        bases, base_raw_starts, base_raw_lengths, probs = index_bases_from_events(events, kmer_index=3)
        self.assertSequenceEqual(bases, list("GATTACAG"))
        self.assertSequenceEqual(base_raw_lengths, [1, 1, 1, 1, 1, 1, 1, 1])
        self.assertSequenceEqual(probs, [1, 1, 1, 1, 1, 1, 1, 1])
        self.assertSequenceEqual(base_raw_starts, [0, 0, 0, 0, 1, 2, 3, 3])
        bases, base_raw_starts, base_raw_lengths, probs = index_bases_from_events(events, kmer_index=4)
        self.assertSequenceEqual(bases, list("GATTACAG"))
        self.assertSequenceEqual(base_raw_lengths, [1, 1, 1, 1, 1, 1, 1, 1])
        self.assertSequenceEqual(probs, [1, 1, 1, 1, 1, 1, 1, 1])
        self.assertSequenceEqual(base_raw_starts, [0, 0, 0, 0, 0, 1, 2, 3])

    def test_get_distance_from_guide_alignment(self):
        data = self.dna3_handle.add_variant_data(number=0)
        basecall_data = self.dna3_handle.add_basecall_alignment_prediction(number=0)
        get_distance_from_guide_alignment(pd.DataFrame(data),pd.DataFrame(basecall_data[0]), reference_index_key="position", minus_strand=self.dna3_handle.aligned_signal.minus_strand)

    def test_add_variant_data(self):
        data = self.dna3_handle.add_variant_data(number=0)
        self.assertEqual(len(data), 8)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp_directory)


class AlignedSignalTest(unittest.TestCase):
    """Test the class AlignedSignal"""

    @classmethod
    def setUpClass(cls):
        super(AlignedSignalTest, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-2])
        cls.dna_file = os.path.join(cls.HOME,
                                    "tests/test_files/minion-reads/canonical/miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch100_read280_strand.fast5")
        cls.modified_file = os.path.join(cls.HOME,
                                         "tests/test_files/minion-reads/methylated/DEAMERNANOPORE_20160805_FNFAD19383_MN16450_sequencing_run_MA_821_R9_gEcoli_MG1655_08_05_16_89825_ch100_read5189_strand.fast5")
        cls.rna_file = os.path.join(cls.HOME,
                                    "tests/test_files/minion-reads/rna_reads/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5")
        cls.handle = AlignedSignal(scaled_signal=[1.1, 2.2, 1.1, 2.2, 1.1, 2.2])

    def test__add_label(self):
        """Test _add_label method"""
        label = np.zeros(4, dtype=[('raw_start', int), ('raw_length', int), ('reference_index', int),
                                   ('posterior_probability', float), ('kmer', 'S5')])
        label["raw_start"] = [0, 1, 2, 3]
        label["raw_length"] = [1, 1, 1, 1]
        label["reference_index"] = [0, 1, 2, 3]
        label["posterior_probability"] = [1, 1, 1, 1]
        label["kmer"] = ["AAT", "A", "B", "C"]

        self.handle.add_label(label, name="test", label_type='label')
        self.handle.add_label(label, name="test2", label_type='prediction')
        self.handle.add_label(label, name="test3", label_type='guide', guide_name="something")
        # catch wrong label type
        with self.assertRaises(AssertionError):
            self.handle.add_label(label, name="test3", label_type='fake')

        with self.assertRaises(KeyError):
            label = np.zeros(0, dtype=[('fake', int), ('raw_length', int), ('reference_index', int),
                                       ('posterior_probability', float), ('kmer', 'S5')])
            self.handle.add_label(label, name="test", label_type="label")

    def test_add_raw_signal(self):
        """Test add_raw_signal method"""
        self.handle.add_raw_signal(np.asanyarray([1, 2, 3, 4, 5, 6]))
        self.handle.add_raw_signal([1, 2, 3, 4, 5, 6])

        with self.assertRaises(AssertionError):
            self.handle.add_raw_signal([1.1, 2.2, 1.1, 4])
            self.handle.add_raw_signal([1.1, 2, 1, 2, 3, 6])

    def test__add_scaled_signal(self):
        """Test _add_scaled_signal method"""
        # add floats as scaled signal
        self.handle._add_scaled_signal(np.asanyarray([1.1, 2.2, 1.1, 2.2, 1.1, 2.2]))
        self.handle._add_scaled_signal([1.1, 2.2, 1.1, 2.2, 1.1, 2.2])
        # throw error if not (probably passed raw ADC counts)
        with self.assertRaises(AssertionError):
            self.handle._add_scaled_signal([1, 2.2, 1.1, 4])
            self.handle._add_scaled_signal([1, 2, 1, 2, 3, 6])

    def test_generate_label_mapping(self):
        label = np.zeros(4, dtype=[('raw_start', int), ('raw_length', int), ('reference_index', int),
                                   ('posterior_probability', float), ('kmer', 'S5')])
        label["raw_start"] = [0, 1, 2, 3]
        label["raw_length"] = [0, 0, 0, 1]
        label["reference_index"] = [0, 1, 2, 3]
        label["posterior_probability"] = [1, 1, 1, 1]
        label["kmer"] = ["AAT", "A", "B", "C"]
        handle = AlignedSignal(scaled_signal=[1.1, 2.2, 1.1, 2.2, 1.1, 2.2])
        # create labels
        handle.add_label(label, name="test", label_type='label')
        handle.add_label(label, name="test2", label_type='label')
        handle.add_label(label, name="test2", label_type='prediction')
        handle.add_label(label, name="test3", label_type='guide', guide_name="something2")
        # make sure we generate the correct mappings
        test = handle.generate_label_mapping(name='test')
        for i, return_tuple in enumerate(test):
            self.assertEqual(return_tuple[0], handle.scaled_signal[i:i + 1])
            self.assertEqual(return_tuple[1], label["kmer"][i])
            self.assertEqual(return_tuple[2], label["posterior_probability"][i])
            self.assertEqual(return_tuple[3], label["reference_index"][i])
        # make sure we generate the correct mappings for all labels added
        test = handle.generate_label_mapping(name='test2')
        for i, return_tuple in enumerate(test):
            self.assertEqual(return_tuple[0], handle.scaled_signal[i:i + 1])
            self.assertEqual(return_tuple[1], label["kmer"][i])
            self.assertEqual(return_tuple[2], label["posterior_probability"][i])
            self.assertEqual(return_tuple[3], label["reference_index"][i])
        # make sure the key exists and the raw data exists
        with self.assertRaises(AssertionError):
            handle.generate_label_mapping(name="test2", scaled=False).__next__()
            handle.generate_label_mapping(name="fake").__next__()

    def test_check_strand_mapping(self):
        label = np.zeros(4, dtype=[('raw_start', int), ('raw_length', int), ('reference_index', int),
                                   ('posterior_probability', float), ('kmer', 'S5')])
        label["raw_start"] = [0, 1, 2, 3]
        label["raw_length"] = [0, 0, 0, 1]
        label["reference_index"] = [0, 1, 2, 3]
        label["posterior_probability"] = [1, 1, 1, 1]
        label["kmer"] = ["AAT", "A", "B", "C"]
        handle = AlignedSignal(scaled_signal=[1.1, 2.2, 1.1, 2.2, 1.1, 2.2])
        handle.check_strand_mapping(label)
        self.assertEqual(handle.minus_strand, False)
        handle = AlignedSignal(scaled_signal=[1.1, 2.2, 1.1, 2.2, 1.1, 2.2], rna=True)
        handle.check_strand_mapping(label)
        self.assertEqual(handle.minus_strand, True)


if __name__ == "__main__":
    unittest.main()
