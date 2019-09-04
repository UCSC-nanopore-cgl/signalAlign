#!/usr/bin/env python
"""
    Place unit tests for event_detection.py
"""
########################################################################
# File: event_detection_test.py
#  executable: event_detection_test.py
# Purpose: event_detection test functions
#
# Author: Andrew Bailey
# History: 12/21/2017 Created
########################################################################
from signalalign.event_detection import *
import unittest
from signalalign.nanoporeRead import NanoporeRead
import tempfile
import shutil
from numpy.lib.recfunctions import rename_fields
import math


class EventDetectTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(EventDetectTests, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.rna_model_file = os.path.join(cls.HOME, "models/testModelR9p4_5mer_acgt_RNA.model")
        cls.dna_template_model_file = os.path.join(cls.HOME, "models/testModelR9p4_5mer_acegt_template.model")

        cls.tmp_directory = tempfile.mkdtemp()

        dna_file = os.path.join(cls.HOME,
                                "tests/minion_test_reads/canonical_ecoli_R9/miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_strand.fast5")
        rna_file = os.path.join(cls.HOME,
                                "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5")
        rna_file2 = os.path.join(cls.HOME,
                                 "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_36_ch_218_strand.fast5")
        # get file locations
        cls.tmp_dna_file = os.path.join(str(cls.tmp_directory), 'test_dna.fast5')
        cls.tmp_rna_file1 = os.path.join(str(cls.tmp_directory), 'test_rna1.fast5')
        cls.tmp_rna_file2 = os.path.join(str(cls.tmp_directory), 'test_rna2.fast5')
        cls.tmp_rna_file3 = os.path.join(str(cls.tmp_directory), 'test_rna3.fast5')

        # copy file to tmp directory
        shutil.copy(dna_file, cls.tmp_dna_file)
        shutil.copy(rna_file, cls.tmp_rna_file1)
        shutil.copy(rna_file, cls.tmp_rna_file2)
        shutil.copy(rna_file2, cls.tmp_rna_file3)

        # create handles
        cls.dna_handle = Fast5(cls.tmp_dna_file, 'r+')
        cls.rna_handle2 = Fast5(cls.tmp_rna_file2, 'r+')

        # clear line for output
        print("")
        print("", file=sys.stderr)

    def test_sequence_from_events(self):
        # """Test sequence from events method"""
        events = np.empty(9, dtype=[('model_state', 'S5'), ('move', '<i4')])
        events["model_state"] = ["AAAAA", "AAAAA", "AAAAT", "AAATA", "AATAA", "ATAAA", "TAAAA", "TAAAA", "AAAAA"]
        events["move"] = [0, 1, 1, 1, 1, 1, 1, 0, 1]
        seq = sequence_from_events(events)
        self.assertEqual(seq, "AAAAAATAAAAA")
        events = np.empty(9, dtype=[('model_state', 'S5'), ('move', '<i4')])
        events["model_state"] = ["AAAAA", "AAAAA", "AAAAT", "AAATA", "AATAA", "ATAAA", "TAAAA", "TAAAA", "AAAAA"]
        events["move"] = [0, 6, 1, 1, 1, 1, 1, 0, 1]
        seq = sequence_from_events(events)
        self.assertEqual(seq, "AAAAAAAAAATAAAAA")
        with self.assertRaises(TypeError):
            sequence_from_events("string")
        with self.assertRaises(KeyError):
            events = np.empty(9, dtype=[('model_state', 'S5')])
            sequence_from_events(events)

    def test_add_start_and_length_modifications(self):
        """test add_raw_start_and_raw_length_to_events and add_start_and_length_to_events methods"""

        # get input data
        sampling_freq = self.dna_handle.sample_rate
        start_time = self.dna_handle.raw_attributes['start_time']
        event_table = self.dna_handle.get_basecall_data()
        self.assertTrue(check_event_table_time(event_table), "Invalid initial start times")

        # add raw fields
        event_table = add_raw_start_and_raw_length_to_events(event_table, start_time, sampling_freq)
        raw_starts = event_table['raw_start'].tolist()
        raw_lengths = event_table['raw_length'].tolist()

        # save old fields
        event_table = rename_fields(event_table, {'start': 'original_start', 'length': 'original_length'})

        # add non-raw fields
        event_table = add_start_and_length_to_events(event_table, start_time, sampling_freq)
        self.assertTrue(check_event_table_time(event_table, min_difference=1.0/sampling_freq), "Invalid modified start times")

        # get fields
        original_starts = event_table['original_start'].tolist()
        original_lengths = event_table['original_length'].tolist()
        starts = event_table['start'].tolist()
        lengths = event_table['length'].tolist()

        # compare elementwise
        places = int(math.log10(sampling_freq)) + 1
        for original_start, start in zip(original_starts, starts):
            self.assertAlmostEqual(original_start, start, places=places)
        for original_length, length in zip(original_lengths, lengths):
            self.assertAlmostEqual(original_length, length, places=places)

    def test_check_event_table_time(self):
        # """test check_event_table_time"""
        events = np.empty(3, dtype=[('start', float), ('length', float)])
        events["start"] = [0, 1, 2]
        events["length"] = [1, 1, 1]
        self.assertTrue(check_event_table_time(events))
        events = np.empty(3, dtype=[('start', float), ('length', float)])
        events["start"] = [0, 2, 2]
        events["length"] = [1, 1, 1]
        self.assertFalse(check_event_table_time(events))

    def test_load_from_raw(self):
        path_to_bin = os.path.join(self.HOME, "bin")
        np_handle = NanoporeRead(os.path.abspath(self.tmp_rna_file3))
        np_handle._initialize_metadata()
        alignment_file = os.path.join(self.HOME, "tests/minion_test_reads/RNA_edge_cases/rna_reads.sam")
        saved_location = load_from_raw(np_handle, alignment_file, self.rna_model_file, path_to_bin,
                                       write_failed_alignments=True, rna=True)
        # close and reopen
        np_handle.close()
        np_handle = NanoporeRead(os.path.abspath(self.tmp_rna_file3))
        # get events and validate
        events = np.array(np_handle.fastFive["/Analyses/Basecall_1D_001/BaseCalled_template/Events"])
        self.assertEqual(events[0]["raw_length"], 11)
        self.assertTrue("/Analyses/Basecall_1D_001/BaseCalled_template/Fastq" in np_handle.fastFive)
        self.assertEqual(saved_location, "/Analyses/Basecall_1D_001")

    def test_run_kmeralign_exe(self):
        path_to_bin = os.path.join(self.HOME, "bin")
        rna_fast5_path = os.path.abspath(self.tmp_rna_file2)
        nuc_sequence = "CAUCCUGCCCUGUGUUAUCCAGUUAUGAGAUAAAAAAUGAAUAUAAGAGUGCUUGUCAUUAUAAAAGUUUUCCUUUUUAUUACCAUCCAAGCCACCAGCUGCCAGCCACCAGCAGCCAGCUGCCAGCACUAGCUUUUUUUUUUUAGCACUUAGUAUUUAGCAGCAUUUAUUAACAGGUACUUUAAGAAUGAUGAAGCAUUGUUUUAAUCUCACUGACUAUGAAGGUUUUAGUUUCUGCUUUUGCAAUUGUGUUUGUGAAAUUUGAAUACUUGCAGGCUUUGUAUGUGAAUAAUUUUAGCGGCUGGUUGGAGAUAAUCCUACGGGAAUUACUUAAAACUGUGCUUUAACUAAAAUGAAUGAGCUUUAAAAUCCCUCCUCCUACUCCAUCAUCAUCCCACUAUUCAUCUUAUCUCAUUAUCAUCAACCUAUCCCACAUCCCUAUCACCACAGCAAUCCAA"
        rna_model_file = self.rna_model_file
        np_handle = NanoporeRead(os.path.abspath(self.tmp_rna_file3))
        np_handle._initialize_metadata()

        dest = "/Analyses/SignalAlign_Basecall_1D_001/BaseCalled_template"
        self.rna_handle2.close()

        status = run_kmeralign_exe(rna_fast5_path, nuc_sequence, rna_model_file, dest, path_to_bin,
                                   write_failed_alignments=True, rna=True)
        rna_handle = Fast5(self.tmp_rna_file2, 'r+')

        events = np.array(rna_handle[dest])

        self.assertEqual(events[0]["raw_length"], 7)
        self.assertTrue(status)

    @classmethod
    def tearDownClass(cls):
        """Remove test fast5 file"""
        shutil.rmtree(cls.tmp_directory)


if __name__ == '__main__':
    unittest.main()
