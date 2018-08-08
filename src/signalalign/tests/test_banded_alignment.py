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
import sys
import numpy as np
import subprocess
from py3helpers.seq_tools import pairwise_alignment_accuracy
from signalalign.fast5 import Fast5
from contextlib import closing
from signalalign import nanoporeRead
from signalalign.hiddenMarkovModel import HmmModel

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

        # clear line for output
        print("")
        print("", file=sys.stderr)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp_directory)
        pass

    def run_kmeralign_exe(self, rna_fast5_path, nuc_sequence, rna_model_file, dest):
        executable = os.path.join(self.HOME, "bin/kmerEventAlign")
        try:
            subprocess.check_call([executable, '-f', rna_fast5_path, '-m', rna_model_file, '-N', nuc_sequence, '-p', dest])
            status = 0
        except Exception as e:
            print("Exception in run_kmeralign: {}".format(e))
            status = -1

        return status

    def compare_event_alignment(self, old_events, new_events):

        old_event_map = nanoporeRead.NanoporeRead.make_event_map(old_events, len(old_events[0]['model_state']))
        new_event_map = nanoporeRead.NanoporeRead.make_event_map(new_events, len(new_events[0]['model_state']))
        new_nuc_seq = nanoporeRead.NanoporeRead.sequence_from_events(new_events)
        old_nuc_seq = nanoporeRead.NanoporeRead.sequence_from_events(old_events)

        self.assertEqual(len(old_event_map), len(old_nuc_seq))
        self.assertEqual(len(new_event_map), len(new_nuc_seq))
        self.assertGreater(pairwise_alignment_accuracy(old_nuc_seq, new_nuc_seq), 0.99)

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
                # print("ERROR  i:{}  oei:{}  nei:{}  oems:{}  nems:{}".format(i, old_event_idx, new_event_idx, old_event['model_state'], new_event['model_state']))
                err_cnt += 1
            else:
                starts.append(abs(old_event['start'] - new_event['start']))
                means.append(abs(old_event['mean'] - new_event['mean']))
                stds.append(abs(old_event['stdv'] - new_event['stdv']))
            # print(err_cnt/(i+1))

        start_diff_avg = np.mean(starts)
        mean_diff_avg = np.mean(means)
        std_diff_avg = np.mean(stds)
        print("total differently-aligned model states: %d/%d (%2.5f%%)" % (err_cnt, len(old_event_map), 100.0 * err_cnt / len(old_event_map)))
        print("DIFF: start:   avg:%.8f  max:%.8f" % (start_diff_avg, max(starts)))
        print("DIFF: mean:    avg:%.8f  max:%.8f" % (mean_diff_avg, max(means)))
        print("DIFF: std:     avg:%.8f  max:%.8f" % (std_diff_avg, max(stds)))

        return start_diff_avg, mean_diff_avg, std_diff_avg

    def run_alignment_comparison(self, src_file_path, dna=True):
        # copy file to tmp directory
        file_name = os.path.basename(src_file_path)
        file_path = os.path.join(self.tmp_directory, file_name)
        shutil.copy(src_file_path, file_path)

        # pre-alignment work
        with closing(Fast5(file_path, read='r+')) as fast5_handle:
            nuc_sequence = fast5_handle.get_fastq(analysis="Basecall_1D", section="template").split()[2]
            dest = fast5_handle.get_analysis_events_path_new(self.UNIT_TEST_NAME)
            #todo verify we don't need this
            # fast5_handle.ensure_path(dest)

        # run kmeralign
        model_file = self.dna_template_model_file if dna else self.rna_model_file
        status = self.run_kmeralign_exe(file_path, nuc_sequence, model_file, dest)
        self.assertEqual(status, 0, "error aligning file {}".format(file_name))

        with closing(Fast5(file_path, read='r')) as fast5_handle:
            new_events = fast5_handle.get_custom_analysis_events(self.UNIT_TEST_NAME)
            if dna:
                old_events = fast5_handle.get_custom_analysis_events(fast5_handle.__default_basecall_1d_analysis__)
            else:
                old_events = fast5_handle.get_resegment_basecall()
                # og_events = fast5_handle.get_custom_analysis_events(fast5_handle.__default_basecall_1d_analysis__)

            print("\nComparing {} events from {}".format("DNA" if dna else "RNA", file_name))

            start_diff_avg, mean_diff_avg, std_diff_avg = self.compare_event_alignment(old_events, new_events)
            self.assertLess(mean_diff_avg, 1.0, "Aligned event means are too varied")
            self.assertLess(std_diff_avg, 1.0, "Aligned event stds are too varied")
            self.assertLess(start_diff_avg, 0.001, "Aligned event start times are too varied")


    # def test_kmeralign_rna1(self):
    #     rna_name = "DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5"
    #     rna_src = os.path.join(self.HOME, "tests/minion_test_reads/RNA_edge_cases", rna_name)
    #     self.run_alignment_comparison(rna_src, dna=False)
    #
    # def test_kmeralign_rna2(self):
    #     rna_name = "DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_36_ch_218_strand.fast5"
    #     rna_src = os.path.join(self.HOME, "tests/minion_test_reads/RNA_edge_cases", rna_name)
    #     self.run_alignment_comparison(rna_src, dna=False)

    def test_kmeralign_dna1(self):
        dna_name = "LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch112_read108_strand.fast5"
        dna_src = os.path.join(self.HOME, "tests/minion_test_reads/1D", dna_name)
        self.run_alignment_comparison(dna_src, dna=True)

    def test_kmeralign_dna2(self):
        dna_name = "LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch6_read347_strand.fast5"
        dna_src = os.path.join(self.HOME, "tests/minion_test_reads/1D", dna_name)
        self.run_alignment_comparison(dna_src, dna=True)

    def test_kmeralign_dna3(self):
        dna_name = "LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch92_read1108_strand.fast5"
        dna_src = os.path.join(self.HOME, "tests/minion_test_reads/1D", dna_name)
        self.run_alignment_comparison(dna_src, dna=True)

if __name__ == '__main__':
    unittest.main()
