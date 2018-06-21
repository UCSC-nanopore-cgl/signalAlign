#!/usr/bin/env python3
"""
    Place unit tests for banded_alignment.py
"""
########################################################################
# File: banded_alignment.py
#  executable: banded_alignment.py
# Purpose: banded alignment functions
#
# Author: Andrew Bailey
# History: 5/31/2018 Created
########################################################################
import unittest
import os
import numpy as np
import threading
import time
from signalalign.banded_alignment import *
from signalalign.fast5 import Fast5
from py3helpers.utils import all_string_permutations, get_random_string


class BandedAlignmentTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(BandedAlignmentTests, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.rna_model_file = os.path.join(cls.HOME, "models/testModelR9p4_5mer_acgt_RNA.model")
        cls.dna_template_model_file = os.path.join(cls.HOME, "models/testModelR9p4_5mer_acegt_template.model")
        cls.rna_model = HmmModel(model_file=cls.rna_model_file)
        cls.dna_model = HmmModel(model_file=cls.dna_template_model_file)

        cls.dna_fast5_path = os.path.join(cls.HOME,
                                          "tests/minion_test_reads/1D/LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch112_read108_strand.fast5")
        cls.dna_handle = Fast5(cls.dna_fast5_path)

    def test_simple_banded_alignment1(self):
        # event_table = self.dna_handle.get_basecall_data()
        # fastq = self.dna_handle.get_fastq(analysis="Basecall_1D", section='template')
        # print(fastq)
        # nuc_sequence = fastq.split('\n')[1]
        # print(nuc_sequence)
        event_table = np.empty(3, dtype=[('start', float), ('length', float),
                                         ('mean', float), ('stdv', float),
                                         ('model_state', 'S5'), ('move', '<i4'),
                                         ('p_model_state', float)])

        event_table["mean"] = [86.9, 75.3, 80.2]
        nuc_sequence = "AAAAATA"

        events, sum_emission = simple_banded_event_align(event_table, self.dna_model, nuc_sequence)
        self.assertSequenceEqual(events['mean'].tolist(), [86.9, 75.3, 80.2])
        self.assertSequenceEqual([bytes.decode(x) for x in events['model_state']], ["AAAAA", "AAAAT", "AAATA"])
        self.assertSequenceEqual(events['move'].tolist(), [0, 1, 1])

    def test_simple_banded_alignment1(self):
        # event_table = self.dna_handle.get_basecall_data()
        # fastq = self.dna_handle.get_fastq(analysis="Basecall_1D", section='template')
        # print(fastq)
        # nuc_sequence = fastq.split('\n')[1]
        # print(nuc_sequence)
        event_table = np.empty(3, dtype=[('start', float), ('length', float),
                                         ('mean', float), ('stdv', float),
                                         ('model_state', 'S5'), ('move', '<i4'),
                                         ('p_model_state', float)])

        event_table["mean"] = [86.9, 75.3, 80.2]
        nuc_sequence = "AAAAATA"

        events, sum_emission = simple_banded_event_align(event_table, self.dna_model, nuc_sequence)
        self.assertSequenceEqual(events['mean'].tolist(), [86.9, 75.3, 80.2])
        self.assertSequenceEqual([bytes.decode(x) for x in events['model_state']], ["AAAAA", "AAAAT", "AAATA"])
        self.assertSequenceEqual(events['move'].tolist(), [0, 1, 1])

    def test_adaptive_banded_alignment1(self):
        # event_table = self.dna_handle.get_basecall_data()
        # fastq = self.dna_handle.get_fastq(analysis="Basecall_1D", section='template')
        # print(fastq)
        # nuc_sequence = fastq.split('\n')[1]
        # print(nuc_sequence)
        event_table = np.empty(3, dtype=[('start', float), ('length', float),
                                         ('mean', float), ('stdv', float),
                                         ('model_state', 'S5'), ('move', '<i4'),
                                         ('p_model_state', float)])

        event_table["mean"] = [86.9, 75.3, 80.2]
        nuc_sequence = "AAAAATA"
        # self.assertEqual(1, 2)
        events, sum_emission = adaptive_banded_simple_event_align(event_table, self.dna_model, nuc_sequence)
        self.assertSequenceEqual(events['mean'].tolist(), [86.9, 75.3, 80.2])
        self.assertSequenceEqual([bytes.decode(x) for x in events['model_state']], ["AAAAA", "AAAAT", "AAATA"])
        self.assertSequenceEqual(events['move'].tolist(), [0, 1, 1])


if __name__ == '__main__':
    unittest.main()
