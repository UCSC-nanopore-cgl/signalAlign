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
from nanotensor.alignedsignal import *


class CreateLabelsTest(unittest.TestCase):
    """Test the class CreateLabels"""

    @classmethod
    def setUpClass(cls):
        super(CreateLabelsTest, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-2])
        cls.dna_file = os.path.join(cls.HOME,
                                    "tests/test_files/minion-reads/canonical/miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch100_read280_strand.fast5")
        cls.modified_file = os.path.join(cls.HOME,
                                         "tests/test_files/minion-reads/methylated/DEAMERNANOPORE_20160805_FNFAD19383_MN16450_sequencing_run_MA_821_R9_gEcoli_MG1655_08_05_16_89825_ch100_read5189_strand.fast5")
        cls.rna_file = os.path.join(cls.HOME,
                                    "tests/test_files/minion-reads/rna_reads/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5")
        cls.fasta = os.path.join(cls.HOME,
                                 "tests/test_files/ecoli_k12_mg1655.fa")

    # TODO
    def test_initialize(self):
        """Test initialize function in CreateLabels"""
        # handle = CreateLabels(self.dna_file)

    # TODO
    def test_add_signal_align_predictions(self):
        """Test add_signal_align_predictions function"""
        pass

    # TODO
    def test_add_mea_labels(self):
        """Test add mea labels"""
        pass

    # TODO
    def test_add_guide_alignment(self):
        """Test add guide alignemnt """
        pass
        # handle = CreateLabels(self.dna_file)
        # handle.add_guide_alignment()

    def test_create_labels_from_guide_alignment(self):
        """Test create_labels_from_guide_alignment"""
        # make sure alignments track correct reference indexes
        test_sam = "@SQ	SN:ref	LN:45\n@SQ	SN:ref2	LN:40\nr001	163	ref	7	30	8M	=	37	39	GATTACTG	*	XX:B:S,12561,2,20,112	MD:Z:6A1"
        events = np.zeros(4, dtype=[('raw_start', int), ('raw_length', int), ('move', int),
                                    ('p_model_state', float), ('model_state', 'S5')])
        events["raw_start"] = [0, 1, 2, 3]
        events["raw_length"] = [1, 1, 1, 1]
        events["move"] = [1, 1, 1, 1]
        events["p_model_state"] = [1, 1, 1, 1]
        events["model_state"] = ["GATTA", "ATTAC", "TTACT", "TACTG"]
        cigar_labels = create_labels_from_guide_alignment(events=events, sam_string=test_sam, kmer_index=2,
                                                          one_ref_indexing=True)[0]
        self.assertEqual("GATTACAG", ''.join([bytes.decode(x) for x in cigar_labels['kmer']]))
        self.assertSequenceEqual([0, 0, 0, 1, 2, 3, 3, 3], cigar_labels['raw_start'].tolist())
        self.assertSequenceEqual([7, 8, 9, 10, 11, 12, 13, 14], cigar_labels['reference_index'].tolist())

        # zero reference indexing and kmer index of 1
        cigar_labels = create_labels_from_guide_alignment(events=events, sam_string=test_sam, kmer_index=1,
                                                          one_ref_indexing=False)[0]
        self.assertEqual("GATTACAG", ''.join([bytes.decode(x) for x in cigar_labels['kmer']]))
        self.assertSequenceEqual([0, 0, 1, 2, 3, 3, 3, 3], cigar_labels['raw_start'].tolist())
        self.assertSequenceEqual([6, 7, 8, 9, 10, 11, 12, 13], cigar_labels['reference_index'].tolist())

        test_header = "@SQ	SN:Chromosome	LN:4641652 \n@PG	ID:bwa	PN:bwa	VN:0.7.15-r1142-dirty	CL:bwa mem -x ont2d /Users/andrewbailey/CLionProjects/nanopore-RNN/signalAlign/bin/test_output/tempFiles_alignment/temp_bwaIndex /Users/andrewbailey/CLionProjects/nanopore-RNN/signalAlign/bin/test_output/tempFiles_alignment/tempFiles_miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_strand/temp_seq_5048dffc-a463-4d84-bd3b-90ca183f488a.fa\n"

        no_mdz = "r001\t163\tChromosome\t1\t30\t7M\t=\t37\t39\tAGCTTTC\t*\tXX:B:S,12561,2,20,112"  # \tMD:Z:6T"
        events = np.zeros(4, dtype=[('raw_start', int), ('raw_length', int), ('move', int),
                                    ('p_model_state', float), ('model_state', 'S5')])
        events["raw_start"] = [0, 1, 2, 3]
        events["raw_length"] = [1, 1, 1, 1]
        events["move"] = [1, 1, 0, 1]
        events["p_model_state"] = [1, 1, 1, 1]
        events["model_state"] = ["AGCTT", "GCTTT", "GCTTT", "CTTTC"]

        with self.assertRaises(AssertionError):
            create_labels_from_guide_alignment(events=events, sam_string=test_header + no_mdz)

        cigar_labels = create_labels_from_guide_alignment(events=events, sam_string=test_header + no_mdz, kmer_index=1,
                                                          one_ref_indexing=False, reference_path=self.fasta)[0]
        self.assertEqual("AGCTTTT", ''.join([bytes.decode(x) for x in cigar_labels['kmer']]))
        self.assertSequenceEqual([0, 0, 1, 3, 3, 3, 3], cigar_labels['raw_start'].tolist())
        self.assertSequenceEqual([0, 1, 2, 3, 4, 5, 6], cigar_labels['reference_index'].tolist())

    def test_index_bases_from_events(self):
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
        self.handle.add_label(label, name="test3", label_type='guide')
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
        """Test generate_label_mapping method"""
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
        handle.add_label(label, name="test3", label_type='guide')
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


if __name__ == "__main__":
    unittest.main()
    raise SystemExit
