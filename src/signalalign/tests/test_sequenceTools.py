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

    def test_find_subsequence_indices(self):
        # no substrings
        indices = find_substring_indices("AAAAA", 'G', offset=0, overlap=True)
        self.assertEqual([x for x in indices], [])
        # yes substrings
        indices = find_substring_indices("AAAAA", 'A', offset=0, overlap=True)
        self.assertEqual([x for x in indices], [0, 1, 2, 3, 4])
        indices = find_substring_indices("AAAAA", 'AA', offset=0, overlap=True)
        self.assertEqual([x for x in indices], [0, 1, 2, 3])
        # test overlap
        indices = find_substring_indices("AAAAA", 'AA', offset=0, overlap=False)
        self.assertEqual([x for x in indices], [0, 2])
        # test offset
        indices = find_substring_indices("ATGCATGC", 'ATGCATGC', offset=1, overlap=True)
        self.assertEqual([x for x in indices], [1])
        indices = find_substring_indices("ATGCATGC", 'ATGCATGC', offset=1, overlap=True)
        self.assertEqual([x for x in indices], [1])
        # compare with gatc modtifs
        indices2 = find_substring_indices("gatcgatc", 'gatc', offset=1, overlap=True)
        self.assertEqual([x for x in indices2], [1, 5])

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


if __name__ == '__main__':
    unittest.main()
