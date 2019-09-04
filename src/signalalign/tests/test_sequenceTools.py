#!/usr/bin/env python3
"""Test sequenceTools.py"""
########################################################################
# File: test_sequenceTools.py
#  executable: test_sequenceTools.py
#
# Author: Andrew Bailey
# History: 5/21/18 Created
########################################################################


import sys
import os
import numpy as np
import pandas as pd
import unittest
import tempfile
import timeit
import filecmp
from itertools import product
from py3helpers.utils import get_random_string, find_substring_indices
from py3helpers.seq_tools import ReverseComplement, ReferenceHandler

from signalalign.utils.sequenceTools import *
from signalalign.utils.fileHandlers import FolderHandler


class TestMakePositionsFiles(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestMakePositionsFiles, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.reference = os.path.join(cls.HOME, "tests/test_sequences/pUC19_SspI_Zymo.fa")
        cls.ecoli_reference = os.path.join(cls.HOME, "tests/test_sequences/E.coli_K12.fasta")
        cls.ambiguity_positions_file = os.path.join(cls.HOME, "tests/test_position_files/test_positions_file3.tsv")

    def test_make_positions_file(self):
        with tempfile.TemporaryDirectory() as tempdir:
            out_path = os.path.join(tempdir, "test.txt")
            file_path = make_positions_file(self.reference, out_path, [["ATTATTGAAG", "ABTATTGAAG"]])
            self.assertTrue(filecmp.cmp(file_path, self.ambiguity_positions_file))
            file_path = make_positions_file(self.ecoli_reference, out_path, [["CCAGG", "CEAGG"]])
            data = CustomAmbiguityPositions.parseAmbiguityFile(file_path)
            rh = ReferenceHandler(self.ecoli_reference)
            forward_strand_data = data[data["strand"] == '+']
            rev_strand_data = data[data["strand"] == '-']

            for i, line in forward_strand_data.iterrows():
                self.assertEqual("CCAGG", rh.get_sequence(line["contig"], line["position"] - 1, line["position"] + 4))

            for i, line in rev_strand_data.iterrows():
                self.assertEqual("CCAGG", reverse_complement(rh.get_sequence(line["contig"], line["position"] - 3, line["position"] + 2)))

    def test_find_gatc_motifs(self):
        indices = find_gatc_motifs("gatcgatc")
        self.assertEqual([x for x in indices], [1, 5])
        indices = find_gatc_motifs(
            "ATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACCATTATTATCATGACATTAACCTATAAAAATAGGCGTATCACGAGGCCCTTTCGTCTCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAAT")
        self.assertEqual([x for x in indices],
                         [460, 601, 1556, 1631, 1642, 1650, 1728, 1740, 1845, 2186, 2204, 2250, 2508, 2525, 2561])
        indices = find_gatc_motifs("ATGC")
        self.assertEqual([x for x in indices], [])

    def test_find_different_char_index(self):
        self.assertRaises(AssertionError, find_different_char_index, "asf", "asdf")
        self.assertRaises(AssertionError, find_different_char_index, "asfd", "asdf")
        index = find_different_char_index("asf", "asd")
        self.assertEqual(index, 2)
        index = find_different_char_index("CCGTT", "CEGTT")
        self.assertEqual(index, 1)

    def test_find_modification_index_and_character(self):
        self.assertRaises(AssertionError, find_modification_index_and_character, "asf", "asdf")
        self.assertRaises(AssertionError, find_modification_index_and_character, "asfd", "asdf")
        self.assertRaises(AssertionError, find_modification_index_and_character, "actg", "actg")

        index, t, f = find_modification_index_and_character("act", "acf")
        self.assertEqual(index, 2)
        self.assertEqual(t, 'T')
        self.assertEqual(f, 'F')
        index, c, e = find_modification_index_and_character("CCGTT", "CEGTT")
        self.assertEqual(index, 1)
        self.assertEqual(c, 'C')
        self.assertEqual(e, 'E')

    def test_find_motifs_sequence_positions(self):
        output = [[i, o, s] for i, o, s in find_motifs_sequence_positions("ATGCATGC", [["ATGC", "ACGC"]])]

        self.assertEqual(len(output), 2)
        self.assertSequenceEqual(output, [[1, "T", "C"], [5, "T", "C"]])

        output = [[i, o, s] for i, o, s in find_motifs_sequence_positions("AAAAAAA", [["AA", "AC"]])]

        self.assertEqual(len(output), 3)
        self.assertSequenceEqual(output, [[1, "A", "C"], [3, "A", "C"], [5, "A", "C"]])

        output = [[i, o, s] for i, o, s in find_motifs_sequence_positions("AAAAAAA", [["AA", "AC"]], overlap=True)]

        self.assertEqual(len(output), 6)
        self.assertSequenceEqual(output,
                                 [[1, "A", "C"], [2, "A", "C"],
                                  [3, "A", "C"], [4, "A", "C"], [5, "A", "C"], [6, "A", "C"]])

    def test_replace_motifs_sequence_positions(self):
        new_seq = replace_motifs_sequence_positions("ATGCATGC", [["ATGC", "ACGC"]])
        self.assertEqual(new_seq, "ACGCACGC")
        new_seq = replace_motifs_sequence_positions("AAAAAAA", [["AA", "AC"]])
        self.assertEqual(new_seq, "ACACACA")
        new_seq = replace_motifs_sequence_positions("AAAAAAA", [["AA", "AC"]], overlap=True)
        self.assertEqual(new_seq, "ACCCCCC")

    def test_replace_periodic_sequence_positions(self):
        new_seq = replace_periodic_sequence_positions("ATGCATGC", 3, 1, "F")
        self.assertEqual(new_seq, "AFGCFTGF")
        new_seq = replace_periodic_sequence_positions("ATGCATGC", 3, 0, "F")
        self.assertEqual(new_seq, "FTGFATFC")

    def test_replace_periodic_reference_positions(self):
        with tempfile.TemporaryDirectory() as tempdir:
            new_path = os.path.join(tempdir, "test.fasta")
            new_fasta = replace_periodic_reference_positions(self.reference, new_path, 10, 0, substitution_char='X')
            for header, comment, sequence in read_fasta(new_fasta):
                for x in sequence[::10]:
                    self.assertEqual(x, "X")

    def test_replace_motifs_reference_positions(self):
        with tempfile.TemporaryDirectory() as tempdir:
            old_path = os.path.join(tempdir, "test.fasta")
            with open(old_path, 'w') as input_fasta:
                print("> header", file=input_fasta)
                print("AAATTTGGGCCCATGATG", file=input_fasta)
            # test normal examples
            new_fasta_path = os.path.join(tempdir, "test2.fasta")
            new_fasta = replace_motif_reference_positions(old_path, new_fasta_path, [["AAA", "AAT"], ["ATG", "ACG"]])
            for header, comment, sequence in read_fasta(new_fasta):
                self.assertEqual(sequence, "AATTTTGGGCCCACGACG")
            # test overlap error
            new_fasta_path = os.path.join(tempdir, "test3.fasta")
            self.assertRaises(AssertionError, replace_motif_reference_positions, old_path, new_fasta_path,
                              [["AAA", "AAT"], ["ATT", "TTT"]])
            new_fasta_path = os.path.join(tempdir, "test4.fasta")
            self.assertRaises(AssertionError, replace_motif_reference_positions, old_path, new_fasta_path,
                              ["ATT", "TTT"])

    def test_samtools_faidx_fasta(self):
        path = samtools_faidx_fasta(self.reference)
        os.remove(path)
        self.assertRaises(AssertionError, samtools_faidx_fasta, "fake.fasta")

    def test_count_all_sequence_kmers(self):
        fake_dict = dict(ATGC=3, TGCA=2, GCAT=3, CATG=2)
        counted_kmers = count_all_sequence_kmers("ATGCATGC", 4, rev_comp=True)
        self.assertEqual(fake_dict, counted_kmers)
        fake_dict = dict(ATGC=2, TGCA=1, GCAT=1, CATG=1)
        counted_kmers = count_all_sequence_kmers("ATGCATGC", 4, rev_comp=False)
        self.assertEqual(fake_dict, counted_kmers)

    def test_get_sequence_kmers(self):
        fake_set = {"ATGC", "TGCA", "GCAT", "CATG"}
        kmer_set = get_sequence_kmers("ATGCATGC", 4, rev_comp=True)
        self.assertEqual(fake_set, kmer_set)
        fake_set = {"AAAA", "TTTT"}
        kmer_set = get_sequence_kmers("AAAAAA", 4, rev_comp=True)
        self.assertEqual(fake_set, kmer_set)
        fake_set = {"AAAA"}
        kmer_set = get_sequence_kmers("AAAAAA", 4, rev_comp=False)
        self.assertEqual(fake_set, kmer_set)

    def test_get_front_back_kmer_overlap(self):
        # CCAGG, CEAGG
        front, back = get_front_back_kmer_overlap(5, 5, 1)
        self.assertEqual(front, 3)
        self.assertEqual(back, 1)
        # CCAGG, CEAGG
        front, back = get_front_back_kmer_overlap(6, 5, 1)
        self.assertEqual(front, 4)
        self.assertEqual(back, 2)
        # CCAGG, CEAGG
        front, back = get_front_back_kmer_overlap(3, 5, 1)
        self.assertEqual(front, 1)
        self.assertEqual(back, -1)
        # CG, EG
        front, back = get_front_back_kmer_overlap(3, 2, 0)
        self.assertEqual(front, 2)
        self.assertEqual(back, 1)
        # CG, EG
        front, back = get_front_back_kmer_overlap(4, 2, 0)
        self.assertEqual(front, 3)
        self.assertEqual(back, 2)
        front, back = get_front_back_kmer_overlap(1, 2, 0)
        self.assertEqual(front, 0)
        self.assertEqual(back, -1)
        front, back = get_front_back_kmer_overlap(2, 2, 0)
        self.assertEqual(front, 1)
        self.assertEqual(back, 0)

        self.assertRaises(AssertionError, get_front_back_kmer_overlap, 0, 2, 0)

    def test_get_motif_kmers(self):
        kmers = get_motif_kmers(["C", "E"], k=2, alphabet="ATGC")
        self.assertSetEqual(set(kmers), {'AE', 'TE', 'GE', 'CE', 'EA', 'ET', 'EG', 'EC'})
        kmers = get_motif_kmers(["ACG", "AEG"], k=2, alphabet="ATGC")
        self.assertSetEqual(set(kmers), {'AE', 'EG'})
        kmers = get_motif_kmers(["CG", "EG"], k=2, alphabet="ATGC")
        self.assertSetEqual(set(kmers), {'AE', 'TE', 'GE', 'CE', 'EG'})
        kmers = get_motif_kmers(["C", "E"], k=3, alphabet="AT")
        self.assertSetEqual(set(kmers), {'AAE', 'ATE', 'TAE', 'TTE', 'AEA',
                                         'AET', 'TEA', 'TET', 'EAA', 'EAT',
                                         'ETA', 'ETT'})
        kmers = get_motif_kmers(["C", "E"], k=1, alphabet="AT")
        self.assertSetEqual(set(kmers), {'E'})
        kmers = self.ccwgg_kmers([x for x in all_string_permutations("ATGC", 5)], 5)
        kmers2 = get_motif_kmers(["CCAGG", "CEAGG"], k=5, alphabet="ATGC")
        kmers3 = get_motif_kmers(["CCTGG", "CETGG"], k=5, alphabet="ATGC")
        kmers2 |= kmers3
        self.assertSetEqual(kmers, kmers2)

    def test_motif_creation_speed(self):
        def ccwgg_motif():
            kmers2 = get_motif_kmers(["CCAGG", "CEAGG"], k=5, alphabet="ATGC")
            kmers3 = get_motif_kmers(["CCTGG", "CETGG"], k=5, alphabet="ATGC")
            kmers2 |= kmers3
            return kmers2

        def ccwgg_motif_old():
            kmers = self.ccwgg_kmers([x for x in all_string_permutations("ATGC", 5)], 5)
            return kmers

        old_time = timeit.timeit(ccwgg_motif_old, number=100)
        new_time = timeit.timeit(ccwgg_motif, number=100)
        self.assertLess(new_time, old_time)

    @staticmethod
    def ccwgg_kmers(sequence_kmers, kmer_length):
        def check_and_add(methyl_kmer):
            normal_kmer = str.translate(methyl_kmer, demethylate)
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(methyl_kmer)

        labeled_kmers = []

        methyl_core1 = "CEAGG"
        methyl_core2 = "CETGG"
        demethylate = str.maketrans("E", "C")

        nucleotides = "ACGT"
        fourmers = [''.join(x) for x in product(nucleotides, repeat=4)]
        threemers = [''.join(x) for x in product(nucleotides, repeat=3)]
        twomers = [''.join(x) for x in product(nucleotides, repeat=2)]
        # NNNNCC*WGGNN

        # NNNNCC*
        if kmer_length == 6:
            for fourmer in fourmers:
                labeled_kmer1 = (fourmer + methyl_core1)[:kmer_length]
                labeled_kmer2 = (fourmer + methyl_core2)[:kmer_length]
                check_and_add(labeled_kmer1)
                check_and_add(labeled_kmer2)

        # NNNCC*W and NNNCC*
        for threemer in threemers:
            labeled_kmer1 = (threemer + methyl_core1)[:kmer_length]
            labeled_kmer2 = (threemer + methyl_core2)[:kmer_length]
            check_and_add(labeled_kmer1)
            check_and_add(labeled_kmer2)

        # NNCC*WG and NNCC*W
        for twomer in twomers:
            labeled_kmer1 = (twomer + methyl_core1)[:kmer_length]
            labeled_kmer2 = (twomer + methyl_core2)[:kmer_length]
            check_and_add(labeled_kmer1)
            check_and_add(labeled_kmer2)
            # C*WGGNN
            if kmer_length == 6:
                labeled_kmer1 = (methyl_core1 + twomer)[1:]
                labeled_kmer2 = (methyl_core2 + twomer)[1:]
                check_and_add(labeled_kmer1)
                check_and_add(labeled_kmer2)

        for onemer in nucleotides:
            # CC*WGGN and C*WGGN
            labeled_kmer1 = methyl_core1 + onemer
            labeled_kmer2 = methyl_core2 + onemer
            if kmer_length == 6:
                check_and_add(labeled_kmer1)
                check_and_add(labeled_kmer2)
            if kmer_length == 5:
                check_and_add(labeled_kmer1[1:])
                check_and_add(labeled_kmer2[1:])
            labeled_kmer1 = (onemer + methyl_core1)[:kmer_length]
            labeled_kmer2 = (onemer + methyl_core2)[:kmer_length]
            check_and_add(labeled_kmer1)
            check_and_add(labeled_kmer2)

        if kmer_length == 5:
            check_and_add(methyl_core1)
            check_and_add(methyl_core2)

        return set(labeled_kmers)


class SignalAlignUtilsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(SignalAlignUtilsTest, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.reference = os.path.join(cls.HOME, "tests/test_sequences/pUC19_SspI_Zymo.fa")
        cls.fast5_dir = os.path.join(cls.HOME, "tests/minion_test_reads/canonical_ecoli_R9")
        cls.files = [
            "miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_strand.fast5",
            "miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch101_read456_strand.fast5",
            "miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch101_read544_strand1.fast5",
            "miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch103_read333_strand1.fast5"]
        cls.fast5_paths = [os.path.join(cls.fast5_dir, f) for f in os.listdir(cls.fast5_dir)
                           if os.path.isfile(os.path.join(cls.fast5_dir, f))]
        cls.ambiguity_positions_file = os.path.join(cls.HOME, "tests/test_position_files/test_positions_file.tsv")
        cls.alignment_file = os.path.join(cls.HOME,
                                          "tests/test_alignments/ecoli1D_test_alignments_sm3/5cc86bac-79fd-4897-8631-8f1c55954a45.sm.backward.tsv")

    def test_processReferenceFasta_positions(self):
        with tempfile.TemporaryDirectory() as tempdir:
            work_folder = FolderHandler()
            work_folder.open_folder(os.path.join(tempdir, "test_outdir"))
            forward_ref, backward_ref = processReferenceFasta(self.reference, work_folder, motifs=None,
                                                              positions_file=self.ambiguity_positions_file,
                                                              name="")
            title, comment, seq = read_fasta(forward_ref).__next__()
            self.assertEqual(seq,
                             "ABTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACCATTATTATCATGACATTAACCTATAAAAATAGGCGTATCACGAGGCCCTTTCGTCTCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAAT")
            title, comment, seq = read_fasta(backward_ref).__next__()
            self.assertEqual(seq,
                             "TBATAACTTCGTAAATAGTCCCAATAACAGAGTACTCGCCTATGTATAAACTTACATAAATCTTTTTATTTGTTTATCCCCAAGGCGCGTGTAAAGGGGCTTTTCACGGTGGACTGCAGATTCTTTGGTAATAATAGTACTGTAATTGGATATTTTTATCCGCATAGTGCTCCGGGAAAGCAGAGCGCGCAAAGCCACTACTGCCACTTTTGGAGACTGTGTACGTCGAGGGCCTCTGCCAGTGTCGAACAGACATTCGCCTACGGCCCTCGTCTGTTCGGGCAGTCCCGCGCAGTCGCCCACAACCGCCCACAGCCCCGACCGAATTGATACGCCGTAGTCTCGTCTAACATGACTCTCACGTGGTATACGCCACACTTTATGGCGTGTCTACGCATTCCTCTTTTATGGCGTAGTCCGCGGTAAGCGGTAAGTCCGACGCGTTGACAACCCTTCCCGCTAGCCACGCCCGGAGAAGCGATAATGCGGTCGACCGCTTTCCCCCTACACGACGTTCCGCTAATTCAACCCATTGCGGTCCCAAAAGGGTCAGTGCTGCAACATTTTGCTGCCGGTCACTTAAGCTCGAGCCATGGGCCCCTAGGAGATCTCAGCTGGACGTCCGTACGTTCGAACCGCATTAGTACCAGTATCGACAAAGGACACACTTTAACAATAGGCGAGTGTTAAGGTGTGTTGTATGCTCGGCCTTCGTATTTCACATTTCGGACCCCACGGATTACTCACTCGATTGAGTGTAATTAACGCAACGCGAGTGACGGGCGAAAGGTCAGCCCTTTGGACAGCACGGTCGACGTAATTACTTAGCCGGTTGCGCGCCCCTCTCCGCCAAACGCATAACCCGCGAGAAGGCGAAGGAGCGAGTGACTGAGCGACGCGAGCCAGCAAGCCGACGCCGCTCGCCATAGTCGAGTGAGTTTCCGCCATTATGCCAATAGGTGTCTTAGTCCCCTATTGCGTCCTTTCTTGTACACTCGTTTTCCGGTCGTTTTCCGGTCCTTGGCATTTTTCCGGCGCAACGACCGCAAAAAGGTATCCGAGGCGGGGGGACTGCTCGTAGTGTTTTTAGCTGCGAGTTCAGTCTCCACCGCTTTGGGCTGTCCTGATATTTCTATGGTCCGCAAAGGGGGACCTTCGAGGGAGCACGCGAGAGGACAAGGCTGGGACGGCGAATGGCCTATGGACAGGCGGAAAGAGGGAAGCCCTTCGCACCGCGAAAGAGTATCGAGTGCGACATCCATAGAGTCAAGCCACATCCAGCAAGCGAGGTTCGACCCGACACACGTGCTTGGGGGGCAAGTCGGGCTGGCGACGCGGAATAGGCCATTGATAGCAGAACTCAGGTTGGGCCATTCTGTGCTGAATAGCGGTGACCGTCGTCGGTGACCATTGTCCTAATCGTCTCGCTCCATACATCCGCCACGATGTCTCAAGAACTTCACCACCGGATTGATGCCGATGTGATCTTCTTGTCATAAACCATAGACGCGAGACGACTTCGGTCAATGGAAGCCTTTTTCTCAACCATCGAGAACTAGGCCGTTTGTTTGGTGGCGACCATCGCCACCAAAAAAACAAACGTTCGTCGTCTAATGCGCGTCTTTTTTTCCTAGAGTTCTTCTAGGAAACTAGAAAAGATGCCCCAGACTGCGAGTCACCTTGCTTTTGAGTGCAATTCCCTAAAACCAGTACTCTAATAGTTTTTCCTAGAAGTGGATCTAGGAAAATTTAATTTTTACTTCAAAATTTAGTTAGATTTCATATATACTCATTTGAACCAGACTGTCAATGGTTACGAATTAGTCACTCCGTGGATAGAGTCGCTAGACAGATAAAGCAAGTAGGTATCAACGGACTGAGGGGCAGCACATCTATTGATGCTATGCCCTCCCGAATGGTAGACCGGGGTCACGACGTTACTATGGCGCTCTGGGTGCGAGTGGCCGAGGTCTAAATAGTCGTTATTTGGTCGGTCGGCCTTCCCGGCTCGCGTCTTCACCAGGACGTTGAAATAGGCGGAGGTAGGTCAGATAATTAACAACGGCCCTTCGATCTCATTCATCAAGCGGTCAATTATCAAACGCGTTGCAACAACGGTAACGATGTCCGTAGCACCACAGTGCGAGCAGCAAACCATACCGAAGTAAGTCGAGGCCAAGGGTTGCTAGTTCCGCTCAATGTACTAGGGGGTACAACACGTTTTTTCGCCAATCGAGGAAGCCAGGAGGCTAGCAACAGTCTTCATTCAACCGGCGTCACAATAGTGAGTACCAATACCGTCGTGACGTATTAAGAGAATGACAGTACGGTAGGCATTCTACGAAAAGACACTGACCACTCATGAGTTGGTTCAGTAAGACTCTTATCACATACGCCGCTGGCTCAACGAGAACGGGCCGCAGTTATGCCCTATTATGGCGCGGTGTATCGTCTTGAAATTTTCACGAGTAGTAACCTTTTGCAAGAAGCCCCGCTTTTGAGAGTTCCTAGAATGGCGACAACTCTAGGTCAAGCTACATTGGGTGAGCACGTGGGTTGACTAGAAGTCGTAGAAAATGAAAGTGGTCGCAAAGACCCACTCGTTTTTGTCCTTCCGTTTTACGGCGTTTTTTCCCTTATTCCCGCTGTGCCTTTACAACTTATGAGTATGAGAAGGAAAAAGTTA")
            forward_ref, backward_ref = processReferenceFasta(self.reference, work_folder, motifs=None,
                                                              positions_file=None, name="")
            title, comment, seq = read_fasta(forward_ref).__next__()
            self.assertEqual(seq,
                             "ATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACCATTATTATCATGACATTAACCTATAAAAATAGGCGTATCACGAGGCCCTTTCGTCTCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAAT")
            self.assertIsNone(backward_ref)
            self.assertEqual(forward_ref, self.reference)

            self.assertRaises(RuntimeError, processReferenceFasta, self.reference, work_folder, motifs="something",
                              positions_file=self.ambiguity_positions_file, name="")

    def test_processReferenceFasta_positions(self):
        with tempfile.TemporaryDirectory() as tempdir:
            work_folder = FolderHandler()
            work_folder.open_folder(os.path.join(tempdir, "test_outdir"))
            forward_ref, backward_ref = processReferenceFasta(self.reference, work_folder, motifs=[["AC", "EC"]],
                                                              positions_file=None, name="")
            title, comment, seq = read_fasta(forward_ref).__next__()
            self.assertEqual(seq.find("AC"), -1)
            self.assertEqual(seq.find("EC"), 42)
            title, comment, seq = read_fasta(backward_ref).__next__()
            self.assertEqual(seq.find("AC"), -1)
            self.assertEqual(seq.find("EC"), 5)

    def test_kmer_iterator(self):
        error_test = kmer_iterator('', 1)
        self.assertRaises(AssertionError, error_test.__next__)
        error_test = kmer_iterator('asdf', 0)
        self.assertRaises(AssertionError, error_test.__next__)

        for x in range(10):
            rand_len = np.random.randint(1, 1000)
            random_string = get_random_string(rand_len)
            rand_k = np.random.randint(1, 10)
            kmer_generator = kmer_iterator(random_string, rand_k)
            for i, kmer in enumerate(kmer_generator):
                self.assertEqual(rand_k, len(kmer))
                self.assertEqual(random_string[i:i + rand_k], kmer)

    def test_reverse_complement(self):
        rev_comp = ReverseComplement(find="ACGTMKRYBVDHNacgtmkrybvdhn", replace="TGCAKMYRVBHDNtgcakmyrvbhdn")
        for x in range(10):
            rand_len = np.random.randint(0, 1000)
            random_dna = get_random_string(rand_len, chars=list(set("ACGTMKRYBVDHN")))

            self.assertEqual(reverse_complement(random_dna, reverse=True, complement=True),
                             rev_comp.reverse_complement(random_dna))
            self.assertEqual(reverse_complement(random_dna, reverse=False, complement=True),
                             rev_comp.complement(random_dna))
            self.assertEqual(reverse_complement(random_dna, reverse=True, complement=False),
                             rev_comp.reverse(random_dna))
            self.assertEqual(reverse_complement(random_dna, reverse=False, complement=False),
                             random_dna)

    def test_count_kmers(self):
        self.assertRaises(AssertionError, count_kmers, '', 1)
        self.assertRaises(AssertionError, count_kmers, 'asdf', 0)
        for x in range(10):
            rand_len = np.random.randint(1, 1000)
            random_string = get_random_string(rand_len)
            rand_k = np.random.randint(1, 10)
            kmer_counts = count_kmers(random_string, rand_k)
            for kmer in kmer_counts.keys():
                self.assertEqual(kmer_counts[kmer], len([x for x in find_substring_indices(random_string, kmer)]),
                                 "random_string: {} \n kmer: {}".format(random_string, kmer))

    def test_parse_full_alignment_file(self):
        data = parse_full_alignment_file(self.alignment_file)
        self.assertEqual(len(data), 16852)
        self.assertRaises(ValueError, parse_full_alignment_file, self.ambiguity_positions_file)


class CustomAmbiguityPositionsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(CustomAmbiguityPositionsTest, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.reference = os.path.join(cls.HOME, "tests/test_sequences/pUC19_SspI_Zymo.fa")
        cls.fast5_dir = os.path.join(cls.HOME, "tests/minion_test_reads/canonical_ecoli_R9")
        cls.files = [
            "miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_strand.fast5",
            "miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch101_read456_strand.fast5",
            "miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch101_read544_strand1.fast5",
            "miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch103_read333_strand1.fast5"]
        cls.fast5_paths = [os.path.join(cls.fast5_dir, f) for f in os.listdir(cls.fast5_dir)
                           if os.path.isfile(os.path.join(cls.fast5_dir, f))]

        cls.ambiguity_positions_file = os.path.join(cls.HOME, "tests/test_position_files/test_positions_file.tsv")
        cls.bad_ambiguity_positions_file = os.path.join(cls.HOME, "tests/test_position_files/test_positions_file2.tsv")

    def test_parseAmbiguityFile(self):
        handle = CustomAmbiguityPositions.parseAmbiguityFile(self.ambiguity_positions_file)
        self.assertRaises(ValueError, CustomAmbiguityPositions.parseAmbiguityFile, self.reference)
        self.assertSequenceEqual(list(handle["position"]), [1, 1, 1, 1])

    def test__get_contig_positions(self):
        handle = CustomAmbiguityPositions(self.ambiguity_positions_file)
        answer = handle.ambig_df.loc[
            (handle.ambig_df["contig"] == 'pUC19') & (handle.ambig_df["strand"] == '+')].drop_duplicates()
        contig_positions = handle._get_contig_positions('pUC19', '+')
        self.assertTrue(answer.equals(contig_positions))

    def test__get_substituted_sequence(self):
        handle = CustomAmbiguityPositions(self.ambiguity_positions_file)
        raw_sequence = handle._get_substituted_sequence("pUC19", "ATGA", "+")
        self.assertEqual("ABGA", raw_sequence)
        raw_sequence = handle._get_substituted_sequence("pUC19", "AAGA", "-")
        self.assertEqual("ABGA", raw_sequence)
        handle = CustomAmbiguityPositions(self.bad_ambiguity_positions_file)
        self.assertRaises(AssertionError, handle._get_substituted_sequence, "pUC19", "AAGA", "+")

    def test_getBackwardSequence(self):
        handle = CustomAmbiguityPositions(self.ambiguity_positions_file)
        raw_sequence = handle.getBackwardSequence("pUC19", "ATGA")
        self.assertEqual("TBCT", raw_sequence)

    def test_getForwardSequence(self):
        handle = CustomAmbiguityPositions(self.ambiguity_positions_file)
        raw_sequence = handle.getForwardSequence("pUC19", "ATGA")
        self.assertEqual("ABGA", raw_sequence)


if __name__ == '__main__':
    unittest.main()
