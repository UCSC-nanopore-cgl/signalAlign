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
import unittest
import os
import numpy as np
import threading
import time
from signalalign.fast5 import Fast5
from signalalign.event_detection import *
import unittest


class EventDetectTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super(EventDetectTests, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.dna_file = os.path.join(cls.HOME, "tests/minion_test_reads/canonical_ecoli_R9/miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_strand.fast5")
        cls.rna_file = os.path.join(cls.HOME,
                                    "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5")
        dna_handle = Fast5(cls.dna_file, 'r+')
        rna_handle = Fast5(cls.rna_file, 'r+')
        cls.dna_handle = dna_handle.create_copy("test_dna.fast5")
        cls.rna_handle = rna_handle.create_copy("test_rna.fast5")

    def test_create_speedy_event_table(self):
        # """Test create_speedy_event_table"""
        for fast5handle in [self.dna_handle, self.rna_handle]:
            sampling_freq = fast5handle.sample_rate
            signal = fast5handle.get_read(raw=True, scale=True)
            start_time = fast5handle.raw_attributes['start_time']
            events = create_speedy_event_table(signal=signal, sampling_freq=sampling_freq, start_time=start_time,
                                               min_width=5, max_width=80, min_gain_per_sample=0.008,
                                               window_width=800)
            events_to_check = np.random.randint(0, len(events), 10)
            # print("SPEEDY")
            # # print(events[:10])
            # print(np.mean(signal[0:5]), events[0]['mean'])
            for x in events_to_check:
                event = events[x]
                signal_mean = np.mean(signal[event["raw_start"]:event["raw_start"] + event["raw_length"]])
                signal_std = np.std(signal[event["raw_start"]:event["raw_start"] + event["raw_length"]])

                self.assertAlmostEqual(event["mean"], signal_mean)
                self.assertAlmostEqual(event["stdv"], signal_std)
                self.assertAlmostEqual(event['raw_start'],
                                       (event['start'] - (start_time / sampling_freq)) * sampling_freq)
                self.assertAlmostEqual(event['raw_length'], event['length'] * sampling_freq)

            with self.assertRaises(TypeError):
                create_speedy_event_table(signal=1, sampling_freq=sampling_freq, start_time=start_time,
                                          min_width=5, max_width=80, min_gain_per_sample=0.008,
                                          window_width=800)
            with self.assertRaises(AssertionError):
                create_speedy_event_table(signal=signal, sampling_freq=sampling_freq, start_time=-120,
                                          min_width=5, max_width=80, min_gain_per_sample=0.008,
                                          window_width=800)
                signal = fast5handle.get_read(raw=True, scale=False)
                create_speedy_event_table(signal=signal, sampling_freq=sampling_freq, start_time=start_time,
                                          min_width=5, max_width=80, min_gain_per_sample=0.008,
                                          window_width=800)

    def test_create_minknow_event_table(self):
        # """Test create_minknow_event_table"""
        for fast5handle in [self.rna_handle, self.dna_handle]:
            sampling_freq = fast5handle.sample_rate
            signal = fast5handle.get_read(raw=True, scale=True)
            start_time = fast5handle.raw_attributes['start_time']
            events = create_minknow_event_table(signal=signal, sampling_freq=sampling_freq, start_time=start_time,
                                                window_lengths=(16, 40), thresholds=(8.0, 4.0), peak_height=1)
            events_to_check = np.random.randint(0, len(events), 10)
            for x in events_to_check:
                event = events[x]
                signal_mean = np.mean(signal[event["raw_start"]:(event["raw_start"] + event["raw_length"])])
                signal_std = np.std(signal[event["raw_start"]:(event["raw_start"] + event["raw_length"])])
                self.assertAlmostEqual(event["mean"], signal_mean)
                self.assertAlmostEqual(event["stdv"], signal_std)
                self.assertAlmostEqual(event['raw_start'],
                                       (event['start'] - (start_time / sampling_freq)) * sampling_freq)
                self.assertAlmostEqual(event['raw_length'], event['length'] * sampling_freq)
            with self.assertRaises(TypeError):
                create_minknow_event_table(signal=1, sampling_freq=sampling_freq, start_time=start_time,
                                           window_lengths=(16, 40), thresholds=(8.0, 4.0), peak_height=1)
            with self.assertRaises(AssertionError):
                create_minknow_event_table(signal=signal, sampling_freq=sampling_freq, start_time=-10,
                                           window_lengths=(16, 40), thresholds=(8.0, 4.0), peak_height=1)
                signal = fast5handle.get_read(raw=True, scale=False)
                create_minknow_event_table(signal=signal, sampling_freq=sampling_freq, start_time=start_time,
                                           window_lengths=(16, 40), thresholds=(8.0, 4.0), peak_height=1)

    def test_create_anchor_kmers1(self):
        # """Test create anchor kmers method with direct matching"""
        new = np.empty(3, dtype=[('start', float), ('length', float),
                                 ('mean', float), ('stdv', float),
                                 ('model_state', 'S5'), ('move', '<i4'),
                                 ('raw_start', int), ('raw_length', int),
                                 ('p_model_state', float)])
        old = np.empty(3, dtype=[('start', float), ('length', float),
                                 ('mean', float), ('stdv', float),
                                 ('model_state', 'S5'), ('move', '<i4'),
                                 ('raw_start', int), ('raw_length', int),
                                 ('p_model_state', float)])
        old["start"] = [0, 1, 2]
        old["length"] = [1, 1, 1]
        old["model_state"] = ["AAAAA", "AAAAT", "AAATT"]
        old["move"] = [1, 1, 1]
        old["p_model_state"] = [0.1, 0.2, 0.2]

        new["start"] = [0, 1, 2]
        new["length"] = [1, 1, 1]

        realignment = create_anchor_kmers(new_events=new, old_events=old)

        self.assertSequenceEqual(realignment['start'].tolist(), new['start'].tolist())
        self.assertSequenceEqual(realignment['length'].tolist(), new['length'].tolist())
        self.assertSequenceEqual(realignment['move'].tolist(), old['move'].tolist())
        self.assertSequenceEqual(realignment['p_model_state'].tolist(), old['p_model_state'].tolist())
        self.assertSequenceEqual(realignment['model_state'].tolist(), old['model_state'].tolist())

    def test_create_anchor_kmers2(self):
        # """Test create anchor kmers method with trimming events if old is longer than new"""

        new = np.empty(3, dtype=[('start', float), ('length', float),
                                 ('mean', float), ('stdv', float),
                                 ('model_state', 'S5'), ('move', '<i4'),
                                 ('raw_start', int), ('raw_length', int),
                                 ('p_model_state', float)])
        old = np.empty(4, dtype=[('start', float), ('length', float),
                                 ('mean', float), ('stdv', float),
                                 ('model_state', 'S5'), ('move', '<i4'),
                                 ('raw_start', int), ('raw_length', int),
                                 ('p_model_state', float)])

        old["start"] = [0, 0.9, 2, 4]
        old["length"] = [0.9, 1.1, 2, 1]
        old["model_state"] = ["AAAAA", "AAAAT", "AAATA", "AAATA"]
        old["move"] = [1, 1, 1, 0]
        old["p_model_state"] = [0.1, 0.2, 0.3, 0.2]

        new["start"] = [1, 2, 3]
        new["length"] = [1, 1, 1]

        realignment = create_anchor_kmers(new_events=new, old_events=old)
        self.assertSequenceEqual(realignment['p_model_state'].tolist(), [0.2, 0.3, 0.3])
        self.assertSequenceEqual(realignment['model_state'].tolist(), [b"AAAAT", b"AAATA", b"AAATA"])
        self.assertSequenceEqual(realignment['move'].tolist(), [2, 1, 0])
        self.assertSequenceEqual(realignment['start'].tolist(), new['start'].tolist())
        self.assertSequenceEqual(realignment['length'].tolist(), new['length'].tolist())

    def test_create_anchor_kmers3(self):
        # """Test create anchor kmers method with trimming events if new is longer than old"""

        new = np.empty(5, dtype=[('start', float), ('length', float),
                                 ('mean', float), ('stdv', float),
                                 ('model_state', 'S5'), ('move', '<i4'),
                                 ('raw_start', int), ('raw_length', int),
                                 ('p_model_state', float)])
        old = np.empty(3, dtype=[('start', float), ('length', float),
                                 ('mean', float), ('stdv', float),
                                 ('model_state', 'S5'), ('move', '<i4'),
                                 ('raw_start', int), ('raw_length', int),
                                 ('p_model_state', float)])

        old["start"] = [2.0, 2.9, 4.0]
        old["length"] = [0.9, 1.1, 1]
        old["model_state"] = ["AAAAA", "AAAAT", "AAATA"]
        old["move"] = [1, 1, 1]
        old["p_model_state"] = [0.1, 0.2, 0.2]

        new["start"] = [1, 2, 3, 4, 5]
        new["length"] = [1, 1, 1, 1, 1]

        realignment = create_anchor_kmers(new_events=new, old_events=old)
        self.assertSequenceEqual(realignment['p_model_state'].tolist(), [0.1, 0.2, 0.2])
        self.assertSequenceEqual(realignment['model_state'].tolist(), [b"AAAAA", b"AAAAT", b"AAATA"])
        self.assertSequenceEqual(realignment['move'].tolist(), [1, 1, 1])
        self.assertSequenceEqual(realignment['start'].tolist(), new['start'].tolist()[1:4])
        self.assertSequenceEqual(realignment['length'].tolist(), new['length'].tolist()[1:4])

    def test_create_anchor_kmers4(self):
        # """Test create anchor kmers method to test multiple moves of same kmer in new event"""

        new = np.empty(4, dtype=[('start', float), ('length', float),
                                 ('mean', float), ('stdv', float),
                                 ('model_state', 'S5'), ('move', '<i4'),
                                 ('raw_start', int), ('raw_length', int),
                                 ('p_model_state', float)])
        old = np.empty(3, dtype=[('start', float), ('length', float),
                                 ('mean', float), ('stdv', float),
                                 ('model_state', 'S5'), ('move', '<i4'),
                                 ('raw_start', int), ('raw_length', int),
                                 ('p_model_state', float)])

        old["start"] = [2.0, 2.9, 4.0]
        old["length"] = [0.9, 1.1, 1]
        old["model_state"] = ["AAAAA", "AAAAA", "AAAAT"]
        old["move"] = [1, 1, 1]
        old["p_model_state"] = [0.1, 0.2, 0.2]

        new["start"] = [1, 2, 4, 5]
        new["length"] = [1, 2, 1, 1]

        realignment = create_anchor_kmers(new_events=new, old_events=old)
        self.assertSequenceEqual(realignment['p_model_state'].tolist(), [0.2, 0.2])
        self.assertSequenceEqual(realignment['model_state'].tolist(), [b"AAAAA", b"AAAAT"])
        self.assertSequenceEqual(realignment['move'].tolist(), [2, 1])
        self.assertSequenceEqual(realignment['start'].tolist(), new['start'].tolist()[1:3])
        self.assertSequenceEqual(realignment['length'].tolist(), new['length'].tolist()[1:3])

    def test_create_anchor_kmers5(self):
        # """Test create anchor kmers method for multiple old events in one new event and multiple new events in one old
        # event"""

        new = np.empty(5, dtype=[('start', float), ('length', float),
                                 ('mean', float), ('stdv', float),
                                 ('model_state', 'S5'), ('move', '<i4'),
                                 ('raw_start', int), ('raw_length', int),
                                 ('p_model_state', float)])
        old = np.empty(5, dtype=[('start', float), ('length', float),
                                 ('mean', float), ('stdv', float),
                                 ('model_state', 'S5'), ('move', '<i4'),
                                 ('raw_start', int), ('raw_length', int),
                                 ('p_model_state', float)])

        old["start"] = [0, 0.2, 0.41, 0.5, 1]
        old["length"] = [0.2, 0.21, 0.09, .5, 1]
        old["model_state"] = ["AAAAA", "AAAAT", "AAATA", "AATAA", "AATAA"]
        old["move"] = [1, 1, 1, 1, 0]
        old["p_model_state"] = [0.1, 0.2, 0.3, 0.4, 0.5]

        new["start"] = [0, 1, 1.1, 1.2, 1.4]
        new["length"] = [1, 0.1, 0.1, 0.2, 1]

        realignment = create_anchor_kmers(new_events=new, old_events=old)
        self.assertSequenceEqual(realignment['p_model_state'].tolist(), [0.4, 0.5, 0.5, 0.5, 0.5])
        self.assertSequenceEqual(realignment['model_state'].tolist(),
                                 [b"AATAA", b"AATAA", b"AATAA", b"AATAA", b"AATAA"])
        self.assertSequenceEqual(realignment['move'].tolist(), [4, 0, 0, 0, 0])
        self.assertSequenceEqual(realignment['start'].tolist(), new['start'].tolist())
        self.assertSequenceEqual(realignment['length'].tolist(), new['length'].tolist())

    def test_create_anchor_kmers6(self):
        # """Test create anchor kmers method for more than 5 skips of old index and if new extends longer than old"""
        new = np.empty(5, dtype=[('start', float), ('length', float),
                                 ('mean', float), ('stdv', float),
                                 ('model_state', 'S5'), ('move', '<i4'),
                                 ('raw_start', int), ('raw_length', int),
                                 ('p_model_state', float)])
        old = np.empty(9, dtype=[('start', float), ('length', float),
                                 ('mean', float), ('stdv', float),
                                 ('model_state', 'S5'), ('move', '<i4'),
                                 ('raw_start', int), ('raw_length', int),
                                 ('p_model_state', float)])

        old["start"] = [0, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 2]
        old["length"] = [1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.4, 0.1]
        old["model_state"] = ["AAAAA", "AAAAA", "AAAAT", "AAATA", "AATAA", "ATAAA", "TAAAA", "TAAAA", "AAAAA"]
        old["move"] = [0, 1, 1, 1, 1, 1, 1, 0, 1]
        old["p_model_state"] = [0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.5, 0.5]

        new["start"] = [0, 1, 2, 2.1, 2.2]
        new["length"] = [1, 1, 0.1, 1, 1]

        realignment = create_anchor_kmers(new_events=new, old_events=old)
        self.assertSequenceEqual(realignment['p_model_state'].tolist(), [0.1, 0.5, 0.5])
        self.assertSequenceEqual(realignment['model_state'].tolist(), [b"AAAAA", b"TAAAA", b"AAAAA"])
        self.assertSequenceEqual(realignment['move'].tolist(), [0, 5, 1])
        self.assertSequenceEqual(realignment['start'].tolist(), new['start'].tolist()[:3])
        self.assertSequenceEqual(realignment['length'].tolist(), new['length'].tolist()[:3])

    def test_create_anchor_kmers7(self):
        # """Test create anchor kmers method for more than 5 skips of old index, select middle kmer and then add moves
        # to next assignment"""

        new = np.empty(5, dtype=[('start', float), ('length', float),
                                 ('mean', float), ('stdv', float),
                                 ('model_state', 'S5'), ('move', '<i4'),
                                 ('raw_start', int), ('raw_length', int),
                                 ('p_model_state', float)])
        old = np.empty(12, dtype=[('start', float), ('length', float),
                                  ('mean', float), ('stdv', float),
                                  ('model_state', 'S5'), ('move', '<i4'),
                                  ('raw_start', int), ('raw_length', int),
                                  ('p_model_state', float)])

        old["start"] = [0, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2]
        old["length"] = [1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
        old["model_state"] = ["AAAAA", "AAAAA", "AAAAT", "AAATA", "AATAA", "ATAAA", "TAAAA", "AAAAA", "AAAAT",
                              "AAATA", "AATAA", "ATAAA"]
        old["move"] = [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        old["p_model_state"] = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

        new["start"] = [0, 1, 1.99, 2.0, 3]
        new["length"] = [1, .99, 0.01, 1, 1]

        realignment = create_anchor_kmers(new_events=new, old_events=old)
        # print(realignment)
        self.assertSequenceEqual(realignment['p_model_state'].tolist(), [0.1, 0.1, 0.1, 0.1])
        self.assertSequenceEqual(realignment['model_state'].tolist(), [b"AAAAA", b"AAAAA", b"AATAA", b"ATAAA"])
        self.assertSequenceEqual(realignment['move'].tolist(), [0, 1, 5, 1])
        self.assertSequenceEqual(realignment['start'].tolist(), new['start'].tolist()[:4])
        self.assertSequenceEqual(realignment['length'].tolist(), new['length'].tolist()[:4])

        with self.assertRaises(TypeError):
            create_anchor_kmers(new_events=1, old_events=old)
        with self.assertRaises(KeyError):
            mean_table = new["mean"]
            create_anchor_kmers(new_events=np.array(range(10)), old_events=old)
            create_anchor_kmers(new_events=mean_table, old_events=old)

    def test_create_anchor_kmers8(self):
        # """Test create anchor kmers method for old event spanning across a new event"""

        new = np.empty(3, dtype=[('start', float), ('length', float),
                                 ('mean', float), ('stdv', float),
                                 ('model_state', 'S5'), ('move', '<i4'),
                                 ('raw_start', int), ('raw_length', int),
                                 ('p_model_state', float)])
        old = np.empty(3, dtype=[('start', float), ('length', float),
                                 ('mean', float), ('stdv', float),
                                 ('model_state', 'S5'), ('move', '<i4'),
                                 ('raw_start', int), ('raw_length', int),
                                 ('p_model_state', float)])

        old["start"] = [0, 1, 2]
        old["length"] = [1, 1, 1]
        old["model_state"] = ["AAAAA", "AAAAA", "AAAAT"]
        old["move"] = [0, 1, 1]
        old["p_model_state"] = [0.1, 0.1, 0.1]

        new["start"] = [0, 1.1, 1.9]
        new["length"] = [1.1, 0.8, 1.1]

        realignment = create_anchor_kmers(new_events=new, old_events=old)
        # print(realignment)
        self.assertSequenceEqual(realignment['p_model_state'].tolist(), [0.1, 0.1, 0.1])
        self.assertSequenceEqual(realignment['model_state'].tolist(), [b"AAAAA", b"AAAAA", b"AAAAT"])
        self.assertSequenceEqual(realignment['move'].tolist(), [1, 0, 1])
        self.assertSequenceEqual(realignment['start'].tolist(), new['start'].tolist())
        self.assertSequenceEqual(realignment['length'].tolist(), new['length'].tolist())

        with self.assertRaises(TypeError):
            create_anchor_kmers(new_events=1, old_events=old)
        with self.assertRaises(KeyError):
            mean_table = new["mean"]
            create_anchor_kmers(new_events=np.array(range(10)), old_events=old)
            create_anchor_kmers(new_events=mean_table, old_events=old)

    def test_create_anchor_kmers9(self):
        # """Test create anchor kmers method for homopolymers"""

        new = np.empty(3, dtype=[('start', float), ('length', float),
                                 ('mean', float), ('stdv', float),
                                 ('model_state', 'S5'), ('move', '<i4'),
                                 ('raw_start', int), ('raw_length', int),
                                 ('p_model_state', float)])
        old = np.empty(5, dtype=[('start', float), ('length', float),
                                 ('mean', float), ('stdv', float),
                                 ('model_state', 'S5'), ('move', '<i4'),
                                 ('raw_start', int), ('raw_length', int),
                                 ('p_model_state', float)])

        old["start"] = [0, 1, 2, 3, 4]
        old["length"] = [1, 1, 1, 1, 1]
        old["model_state"] = ["AAAAA", "AAAAA", "AAAAA", "AAAAA", "AAAAA"]
        old["move"] = [0, 1, 1, 0, 1]
        old["p_model_state"] = [0.1, 0.1, 0.1, 0.1, 0.1]

        new["start"] = [0, 1.1, 4.5]
        new["length"] = [1.1, 3.4, 1]

        realignment = create_anchor_kmers(new_events=new, old_events=old)
        # print(realignment)
        self.assertSequenceEqual(realignment['p_model_state'].tolist(), [0.1, 0.1, 0.1])
        self.assertSequenceEqual(realignment['model_state'].tolist(), [b"AAAAA", b"AAAAA", b"AAAAA"])
        self.assertSequenceEqual(realignment['move'].tolist(), [1, 2, 0])
        self.assertSequenceEqual(realignment['start'].tolist(), new['start'].tolist())
        self.assertSequenceEqual(realignment['length'].tolist(), new['length'].tolist())

        with self.assertRaises(TypeError):
            create_anchor_kmers(new_events=1, old_events=old)
        with self.assertRaises(KeyError):
            mean_table = new["mean"]
            create_anchor_kmers(new_events=np.array(range(10)), old_events=old)
            create_anchor_kmers(new_events=mean_table, old_events=old)

    def test_resegment_reads(self):
        # """Test resegment_reads method"""
        minknow_params = dict(window_lengths=(5, 10), thresholds=(2.0, 1.1), peak_height=1.2)
        speedy_params = dict(min_width=5, max_width=80, min_gain_per_sample=0.008, window_width=800)
        with self.assertRaises(AssertionError):
            resegment_reads("fakepath/path", speedy_params, speedy=True, overwrite=True)
        with self.assertRaises(TypeError):
            resegment_reads(self.dna_file, speedy_params, speedy=False, overwrite=True)

        for fast5_file in [self.dna_file, self.rna_file]:
            resegment_reads(fast5_file, minknow_params, speedy=False, overwrite=True)
            fasthandle = resegment_reads(fast5_file, speedy_params, speedy=True, overwrite=True)
            # TODO make sure this is working with test files
            fasthandle.delete("Analyses/ReSegmentBasecall_000")

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

    def test_index_to_time(self):
        # """test index_to_time method"""
        # get input data
        sampling_freq = self.rna_handle.sample_rate
        start_time = self.rna_handle.raw_attributes['start_time']
        event_table = self.rna_handle.get_basecall_data()
        # run method
        new_table = index_to_time(event_table, sampling_freq=sampling_freq, start_time=start_time)
        start = event_table["start"] / sampling_freq + (start_time / sampling_freq)
        length = event_table["length"] / float(sampling_freq)
        # compare elementwise
        self.assertSequenceEqual(new_table["start"].tolist(), start.tolist())
        self.assertSequenceEqual(new_table["length"].tolist(), length.tolist())

        with self.assertRaises(AssertionError):
            index_to_time(event_table, start_time=start_time)
        with self.assertRaises(AssertionError):
            index_to_time(event_table, sampling_freq=sampling_freq)
        with self.assertRaises(KeyError):
            index_to_time(np.array([1, 2, 3]), sampling_freq=sampling_freq, start_time=start_time)
        with self.assertRaises(TypeError):
            index_to_time("Numpy", sampling_freq=sampling_freq, start_time=start_time)
        with self.assertRaises(AssertionError):
            event_table = self.dna_handle.get_basecall_data()
            index_to_time(event_table, sampling_freq=sampling_freq, start_time=start_time)

    def test_time_to_index(self):
        # """test time_to_index method"""
        # get input data

        sampling_freq = self.dna_handle.sample_rate
        start_time = self.dna_handle.raw_attributes['start_time']
        event_table = self.dna_handle.get_basecall_data()
        # run method
        start = np.round((event_table["start"] - (start_time / float(sampling_freq))) * sampling_freq)
        length = np.round(event_table["length"] * sampling_freq)

        new_table = time_to_index(event_table, sampling_freq=sampling_freq, start_time=start_time)
        # compare elementwise
        self.assertSequenceEqual(new_table["start"][0:100].tolist(), start[0:100].tolist())
        self.assertSequenceEqual(new_table["length"][0:100].tolist(), length[0:100].tolist())

        with self.assertRaises(AssertionError):
            time_to_index(event_table, start_time=start_time)
        with self.assertRaises(AssertionError):
            time_to_index(event_table, sampling_freq=sampling_freq)
        with self.assertRaises(KeyError):
            time_to_index(np.array([1, 2, 3]), sampling_freq=sampling_freq, start_time=start_time)
        with self.assertRaises(TypeError):
            time_to_index("Numpy", sampling_freq=sampling_freq, start_time=start_time)
        with self.assertRaises(AssertionError):
            event_table = self.rna_handle.get_basecall_data()
            time_to_index(event_table, sampling_freq=sampling_freq, start_time=start_time)


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

    def test_get_resegment_accuracy(self):
        # """Test get_resegment_accuracy"""
        minknow_params = dict(window_lengths=(5, 10), thresholds=(2.0, 1.1), peak_height=1.2)
        f5fh = resegment_reads("test_rna.fast5", minknow_params, speedy=False, overwrite=True)
        self.assertAlmostEqual(get_resegment_accuracy(f5fh), 1.0)

    def test_create_minknow_events_from_fast5(self):
        for path in ["test_rna.fast5", "test_dna.fast5"]:
            events, f5handle = create_minknow_events_from_fast5(path)
            passing = check_numpy_table(events, req_fields=('start', 'length', 'mean', 'stdv', 'model_state', 'move', 'p_model_state'))
            self.assertTrue(passing)
            self.assertIsInstance(f5handle, Fast5)

    @classmethod
    def tearDownClass(cls):
        """Remove test fast5 file"""
        os.remove("test_rna.fast5")
        os.remove("test_dna.fast5")


if __name__ == '__main__':
    unittest.main()
