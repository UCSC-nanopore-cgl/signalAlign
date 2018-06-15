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
import unittest
import tempfile
import timeit
from itertools import product
from collections import defaultdict
from scipy import sparse
from signalalign.utils.sequenceTools import *
from signalalign.utils import CustomAmbiguityPositions
from signalalign.utils.fileHandlers import FolderHandler

class TestMakePositionsFiles(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestMakePositionsFiles, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.reference = os.path.join(cls.HOME, "tests/test_sequences/pUC19_SspI_Zymo.fa")
        # cls.reference = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/tests/test_sequences/pUC19_SspI_Zymo.fa"

    def test_find_gatc_motifs(self):
        indices = find_gatc_motifs("gatcgatc")
        self.assertEqual([x for x in indices], [1, 5])
        indices = find_gatc_motifs("ATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACCATTATTATCATGACATTAACCTATAAAAATAGGCGTATCACGAGGCCCTTTCGTCTCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAAT")
        self.assertEqual([x for x in indices], [460, 601, 1556, 1631, 1642, 1650, 1728, 1740, 1845, 2186, 2204, 2250, 2508, 2525, 2561])
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

    def replace_motifs_sequence_positions(self):
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
            self.assertRaises(AssertionError, replace_motif_reference_positions, old_path, new_fasta_path, [["AAA", "AAT"], ["ATT", "TTT"]])
            new_fasta_path = os.path.join(tempdir, "test4.fasta")
            self.assertRaises(AssertionError, replace_motif_reference_positions, old_path, new_fasta_path, ["ATT", "TTT"])

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


if __name__ == '__main__':
    unittest.main()
