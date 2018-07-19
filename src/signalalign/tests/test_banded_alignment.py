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
import tempfile
import shutil
import os
import numpy as np
import math
import threading
import time
from signalalign.banded_alignment import *
from signalalign.fast5 import Fast5
from py3helpers.utils import all_string_permutations, get_random_string
import sys
from contextlib import closing
from signalalign import nanoporeRead
import kmeralign

class BandedAlignmentTests(unittest.TestCase):
    UNIT_TEST_NAME = "PyUnittest"

    @classmethod
    def setUpClass(cls):
        super(BandedAlignmentTests, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.rna_model_file = os.path.join(cls.HOME, "models/testModelR9p4_5mer_acgt_RNA.model")
        cls.dna_template_model_file = os.path.join(cls.HOME, "models/testModelR9p4_5mer_acegt_template.model")
        cls.rna_model = HmmModel(model_file=cls.rna_model_file)
        cls.dna_model = HmmModel(model_file=cls.dna_template_model_file)

        cls.tmp_directory = tempfile.mkdtemp()

        dna_name = "LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch112_read108_strand.fast5"
        dna_src = os.path.join(cls.HOME, "tests/minion_test_reads/1D", dna_name)
        dna_dst = os.path.join(cls.tmp_directory, dna_name)
        shutil.copy(dna_src, dna_dst)
        cls.dna_fast5_path =  dna_dst

        rna_name = "DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5"
        rna_src = os.path.join(cls.HOME, "tests/minion_test_reads/RNA_edge_cases", rna_name)
        rna_dst = os.path.join(cls.tmp_directory, rna_name)
        shutil.copy(rna_src, rna_dst)
        cls.rna_fast5_path =  rna_dst

        # clear line for output
        print("")
        print("", file=sys.stderr)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp_directory)
        pass

    @staticmethod
    def run_kmeralign(rna_fast5_path, nuc_sequence, rna_model_file, dest):

        # because the file handles stay open when I use the cython integration :/
        try:
            subprocess.check_call(['python', '-c', 'import kmeralign ; kmeralign.load_from_raw(fast5_path="{}", nuc_sequence="{}", template_model_file="{}", path_in_fast5="{}")'.format(rna_fast5_path, nuc_sequence, rna_model_file, dest)])
            status = 0
        except Exception as e:
            print("Exception in run_kmeralign: {}".format(e))
            status = -1

        # todo this is how we want to invoke this
        # status = kmeralign.load_from_raw(fast5_path=self.rna_fast5_path, template_model_file=self.rna_model_file,
        #                                  nuc_sequence=nuc_sequence, path_in_fast5=dest)

        return status

    def compare_event_alignment(self, old_events, new_events):

        old_event_map = nanoporeRead.NanoporeRead.make_event_map(old_events, len(old_events[0]['model_state']))
        new_event_map = nanoporeRead.NanoporeRead.make_event_map(new_events, len(new_events[0]['model_state']))
        self.assertEquals(len(old_event_map), len(new_event_map))

        starts = []
        means = []
        stds = []
        err_cnt = 0.0
        for i in range(len(old_event_map)):
            old_event_idx = old_event_map[i]
            new_event_idx = new_event_map[i]

            old_event = old_events[old_event_idx]
            new_event = new_events[new_event_idx]
            if (old_event['model_state'] != new_event['model_state']):
                print("ERROR  i:{}  oei:{}  nei:{}  oems:{}  nems:{}".format(i, old_event_idx, new_event_idx, old_event['model_state'], new_event['model_state']))
                err_cnt += 1
            else:
                starts.append(abs(old_event['start'] - new_event['start']))
                means.append(abs(old_event['mean'] - new_event['mean']))
                stds.append(abs(old_event['stdv'] - new_event['stdv']))

        start_diff_avg = np.mean(starts)
        mean_diff_avg = np.mean(means)
        std_diff_avg = np.mean(stds)
        print("total differently-aligned model states: %d/%d (%2.5f%%)\n" % (err_cnt, len(old_event_map), 100.0 * err_cnt / len(old_event_map)))
        print("DIFF: start:   avg:%.8f  max:%.8f" % (start_diff_avg, max(starts)))
        print("DIFF: mean:    avg:%.8f  max:%.8f" % (mean_diff_avg, max(means)))
        print("DIFF: std:     avg:%.8f  max:%.8f" % (std_diff_avg, max(stds)))

        return start_diff_avg, mean_diff_avg, std_diff_avg

    def test_kmeralign_rna(self):
        with closing(Fast5(self.rna_fast5_path, read='r+')) as fast5_handle:
            nuc_sequence = fast5_handle.get_fastq(analysis="Basecall_1D", section="template").split()[2]
            dest = fast5_handle.get_analysis_events_path_new(self.UNIT_TEST_NAME)
            fast5_handle.ensure_path(dest)

        status = self.run_kmeralign(self.rna_fast5_path, nuc_sequence, self.rna_model_file, dest)
        self.assertEquals(status, 0, "error aligning DNA file")

        with closing(Fast5(self.rna_fast5_path, read='r')) as fast5_handle:
            new_events = fast5_handle.get_custom_analysis_events(self.UNIT_TEST_NAME)
            old_events = fast5_handle.get_custom_analysis_events(fast5_handle.__default_basecall_1d_analysis__)
            self.compare_event_alignment(old_events, new_events)

            start_diff_avg, mean_diff_avg, std_diff_avg = self.compare_event_alignment(old_events, new_events)
            self.assertTrue(mean_diff_avg < 1.0, "Aligned event means are too varied")
            self.assertTrue(std_diff_avg < 1.0, "Aligned event stds are too varied")
            self.assertTrue(start_diff_avg < 0.001, "Aligned event start times are too varied")

    def test_kmeralign_dna(self):
        with closing(Fast5(self.dna_fast5_path, read='r+')) as fast5_handle:
            nuc_sequence = fast5_handle.get_fastq(analysis="Basecall_1D", section="template").split()[2]
            dest = fast5_handle.get_analysis_events_path_new(self.UNIT_TEST_NAME)
            fast5_handle.ensure_path(dest)

        status = self.run_kmeralign(self.dna_fast5_path, nuc_sequence, self.dna_template_model_file, dest)
        self.assertEquals(status, 0, "error aligning DNA file")

        with closing(Fast5(self.dna_fast5_path, read='r')) as fast5_handle:
            new_events = fast5_handle.get_custom_analysis_events(self.UNIT_TEST_NAME)
            old_events = fast5_handle.get_custom_analysis_events(fast5_handle.__default_basecall_1d_analysis__)

            start_diff_avg, mean_diff_avg, std_diff_avg = self.compare_event_alignment(old_events, new_events)
            self.assertTrue(mean_diff_avg < 1.0, "Aligned event means are too varied")
            self.assertTrue(std_diff_avg < 1.0, "Aligned event stds are too varied")
            self.assertTrue(start_diff_avg < 0.001, "Aligned event start times are too varied")




if __name__ == '__main__':
    unittest.main()
