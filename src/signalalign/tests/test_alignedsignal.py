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
from signalalign.alignedsignal import *
from signalalign.visualization.plot_labelled_read import PlotSignal
from signalalign.signalAlignment import SignalAlignment, create_signalAlignment_args
from py3helpers.utils import merge_dicts


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


if __name__ == "__main__":
    unittest.main()
    raise SystemExit
