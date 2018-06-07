#!/usr/bin/env python3
"""
    Place unit tests for fast5.py
"""
########################################################################
# File: fast5_test.py
#  executable: fast5_test.py
# Purpose: fast5 test functions
#
# Author: Andrew Bailey
# History: 12/20/2017 Created
########################################################################
import unittest
import os
import numpy as np
import threading
import time
from signalalign.fast5 import Fast5


class Fast5Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(Fast5Test, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        fast5_file = os.path.join(cls.HOME, "tests/minion_test_reads/canonical_ecoli_R9/miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_strand.fast5")
        fast5handle = Fast5(fast5_file, 'r+')
        cls.fast5handle = fast5handle.create_copy("test.fast5")

    def test_get_fastq(self):
        # """Test get_fastq method of Fast5 class"""
        self.fast5handle.get_fastq()

    def test_get_basecall_data(self):
        # """Test get_basecall_data method of Fast5 class"""
        self.fast5handle.get_basecall_data()

    def test_get_read_stats(self):
        # """Test get_read_stats from Fast5 class"""
        self.fast5handle.get_read_stats()

    def test_get_read(self):
        # """Test get_read from Fast5 class"""
        self.fast5handle.get_read(raw=True, scale=False)

    def test_get_corrected_events(self):
        # """Test get_corrected_events from Fast5 class"""
        self.fast5handle.get_corrected_events()

    def test_delete(self):
        # """Test delete fast5 file"""
        channel_id = self.fast5handle["UniqueGlobalKey/channel_id"].attrs
        self.fast5handle.delete("UniqueGlobalKey/channel_id")
        self.fast5handle.delete("UniqueGlobalKey/fakedata", ignore=True)
        with self.assertRaises(KeyError):
            error = self.fast5handle["UniqueGlobalKey/channel_id"]
            self.fast5handle.delete("UniqueGlobalKey/channel_id")
            self.fast5handle.delete("UniqueGlobalKey/fakedata")
        self.fast5handle._add_attrs(channel_id, "UniqueGlobalKey/channel_id", convert=None)

    def test_test_event_table(self):
        # """Test the method test_event_table"""
        pass

    @classmethod
    def tearDownClass(cls):
        """Remove test fast5 file"""
        os.remove("test.fast5")


if __name__ == '__main__':
    unittest.main()
