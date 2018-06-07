#!/usr/bin/env python3
"""Test maximum expected accuracy alignment algorithm"""
########################################################################
# File: mea_algorithm_test.py
#  executable: mea_algorithm_test.py
#
# Author: Andrew Bailey
# History: 1/25/18 Created
########################################################################


import sys
import numpy as np
import unittest
from collections import defaultdict
from scipy import sparse
from signalalign.mea_algorithm import *
from py3helpers.utils import time_it


class MeaTest(unittest.TestCase):
    """Test the functions in mea_algorithm.py"""

    def test_maximum_expected_accuracy_alignment(self):
        # """Test maximum_expected_accuracy_alignment function"""
        # ref x event
        posterior_matrix = [[0.2, 0.2, 0.3, 0.0, 0.0],
                            [0.3, 0.5, 0.1, 0.0, 0.0],
                            [0.2, 0.3, 0.0, 0.0, 0.0],
                            [0.2, 0.0, 0.3, 0.4, 0.2],
                            [0.1, 0.0, 0.3, 0.1, 0.5]]

        posterior_matrix = np.asanyarray(posterior_matrix).T
        # correct input
        shortest_ref_per_event = [0, 0, 0, 3, 3]
        forward_edges = maximum_expected_accuracy_alignment(posterior_matrix, shortest_ref_per_event, return_all=True)
        # trim unnecessary edges
        self.assertEqual(3, len(forward_edges))
        # 0.2->0.5->0.1 = 0.7 dont count horizontal move through 0.1
        self.assertAlmostEqual(sum_forward_edge_accuracy(forward_edges[0]) - 0.1, forward_edges[0][3])
        self.assertAlmostEqual(0.7, forward_edges[0][3])
        # 0.2->0.5->0.1->0.4->0.2 = 1.1 (Don't count last horizontal move and move from 0.4 to 0.2)
        self.assertAlmostEqual(sum_forward_edge_accuracy(forward_edges[1]) - 0.2 - 0.1, forward_edges[1][3])
        self.assertAlmostEqual(1.1, forward_edges[1][3])
        # 0.2->0.5->0.1->0.4->0.5 = 1.6 (Don't count horizontal move from 0.5 to 0.1)
        self.assertAlmostEqual(sum_forward_edge_accuracy(forward_edges[2]) - 0.1, forward_edges[2][3])
        self.assertAlmostEqual(1.6, forward_edges[2][3])

        # test passing through a sparse matrix
        forward_edges = maximum_expected_accuracy_alignment(sparse.coo_matrix(posterior_matrix), shortest_ref_per_event,
                                                            sparse_posterior_matrix=True,
                                                            return_all=True)
        self.assertEqual(3, len(forward_edges))
        # 0.2->0.5->0.1 = 0.7 dont count horizontal move through 0.1
        self.assertAlmostEqual(sum_forward_edge_accuracy(forward_edges[0]) - 0.1, forward_edges[0][3])
        self.assertAlmostEqual(0.7, forward_edges[0][3])
        # 0.2->0.5->0.1->0.4->0.2 = 1.1 (Don't count last horizontal move and move from 0.4 to 0.2)
        self.assertAlmostEqual(sum_forward_edge_accuracy(forward_edges[1]) - 0.2 - 0.1, forward_edges[1][3])
        self.assertAlmostEqual(1.1, forward_edges[1][3])
        # 0.2->0.5->0.1->0.4->0.5 = 1.6 (Don't count horizontal move from 0.5 to 0.1)
        self.assertAlmostEqual(sum_forward_edge_accuracy(forward_edges[2]) - 0.1, forward_edges[2][3])
        self.assertAlmostEqual(1.6, forward_edges[2][3])

        # incorrect min ref lengths
        shortest_ref_per_event = [0, 1, 1, 1, 1]
        forward_edges = maximum_expected_accuracy_alignment(posterior_matrix, shortest_ref_per_event, return_all=True)
        self.assertEqual(4, len(forward_edges))

    def test_binary_search_for_edge(self):
        # """Test binary_search_for_edge"""
        forward_edges = [[0, 1, .1, .1], [1, 1, .1, .1], [2, 1, .1, .1], [3, 1, .1, .1], [4, 1, .1, .1],
                         [5, 1, .1, .1], [6, 1, .1, .1]]

        binary_edge = binary_search_for_edge(forward_edges, 1.1, 1, 0.1)
        self.assertEqual([1.1, 1, .1, .2, [1, 1, .1, .1]], binary_edge)
        random_indexes = np.random.uniform(0, 7, 10)
        for index in random_indexes:
            edge = slow_search_for_edge(forward_edges, index, 1, 0.1)
            binary_edge = binary_search_for_edge(forward_edges, index, 1, 0.1)
            self.assertEqual(edge, binary_edge)

        with self.assertRaises(AssertionError):
            binary_search_for_edge(forward_edges, -1, 1, 0.1)

    def test_mae_random_matrix(self):
        # """Create random alignment matrix for mae alignment"""
        max_size = 40
        for j in range(20):
            row, col = np.random.randint(int(max_size/2), max_size, 2)
            # create random inputs
            posterior_matrix, shortest_ref_per_event = create_random_prob_matrix(row=row, col=col)
            # test 3 implementations
            most_probable_edge = maximum_expected_accuracy_alignment(posterior_matrix, shortest_ref_per_event)
            another_probable_edge = mea_slower(posterior_matrix, shortest_ref_per_event)
            yet_another_implementation = mea_slow(posterior_matrix, shortest_ref_per_event)

            self.assertAlmostEqual(most_probable_edge[3], another_probable_edge[3])
            self.assertAlmostEqual(yet_another_implementation[3], most_probable_edge[3])

            # check total probability and traceback
            ref_indx = most_probable_edge[0]
            event_indx = most_probable_edge[1]
            prob = most_probable_edge[2]
            sum_prob = 0
            total_prob = most_probable_edge[3]
            prev_event = most_probable_edge[4]
            while prev_event:
                # step or stay with reference
                self.assertGreaterEqual(ref_indx, prev_event[0])
                # must step for each event
                self.assertGreater(event_indx, prev_event[1])
                # gather correct probabilities
                if ref_indx != prev_event[0]:
                    sum_prob += prob
                ref_indx = prev_event[0]
                event_indx = prev_event[1]
                prob = prev_event[2]
                prev_event = prev_event[4]
            # include first probability
            sum_prob += prob
            self.assertAlmostEqual(sum_prob, total_prob)

    def test_get_indexes_from_best_path(self):
        # """test get_get_indexes_from_best_path"""
        fake_mea = [4, 4, 1, 0.1, [3, 3, 0.9, 0.1, [2, 2, 0.8, 0.1, [1, 1, 0.7, 0.1, [0, 0, 0.6, 0.6, None]]]]]
        alignment = get_indexes_from_best_path(fake_mea)
        self.assertEqual(alignment, [[0, 0], [1, 1], [2, 2], [3, 3], [4, 4]])
        fake_mea = [0.1, [3, 3, 0.9, 0.1, [2, 2, 0.8, 0.1, [1, 1, 0.7, 0.1, [0, 0, 0.6, 0.6, None]]]]]
        with self.assertRaises(IndexError):
            get_indexes_from_best_path(fake_mea)

    def test_get_events_from_path(self):
        # """Test get_events_from_path"""
        path = [[0, 0], [1, 1], [2, 2], [3, 3]]
        event1 = np.zeros(4, dtype=[('contig', 'S10'), ('reference_index', '<i8'), ('reference_kmer', 'S5'),
                                    ('strand', 'S1'),
                                    ('event_index', '<i8'), ('event_mean', '<f8'), ('event_noise', '<f8'),
                                    ('event_duration', '<f8'), ('aligned_kmer', 'S5'),
                                    ('scaled_mean_current', '<f8'), ('scaled_noise', '<f8'),
                                    ('posterior_probability', '<f8'), ('descaled_event_mean', '<f8'),
                                    ('ont_model_mean', '<f8'), ('path_kmer', 'S5')])
        event_matrix = [[event1[0], 0, 0, 0],
                        [0, event1[0], 0, 0],
                        [0, 0, event1[0], 0],
                        [0, 0, 0, event1[0]]]
        events = get_events_from_path(event_matrix, path)
        self.assertSequenceEqual(events.tolist(), event1.tolist())
        event_matrix = [[0, 0, 0, 0],
                        [0, 0, 0, 0],
                        [0, 0, 0, 0]]

        self.assertRaises(IndexError, get_events_from_path, event_matrix, path)

    def test_get_mea_params_from_events(self):
        # """Test get_mea_params_from_events"""
        for _ in range(20):
            max_size = 20
            row, col = np.random.randint(int(max_size/2), max_size, 2)
            # create random probability matrix
            true_posterior_matrix, true_shortest_ref_per_event = create_random_prob_matrix(row=row, col=col,
                                                                                                gaps=False)
            # print(true_posterior_matrix.T)
            # print(true_shortest_ref_per_event)
            # generate events and event matrix to match
            events, true_event_matrix = generate_events_from_probability_matrix(true_posterior_matrix)
            # get events from random signal align output
            # _, time = time_it(get_mea_params_from_events, events)
            # print(time)
            posterior_matrix, shortest_ref, event_matrix = get_mea_params_from_events(events)
            self.assertSequenceEqual(posterior_matrix.tolist(), true_posterior_matrix.tolist())
            self.assertSequenceEqual(shortest_ref, true_shortest_ref_per_event.tolist())
            self.assertSequenceEqual(event_matrix, true_event_matrix)
        with self.assertRaises(KeyError):
            events = np.zeros(4, dtype=[('reference_index', '<i8'), ('reference_kmer', 'S5'),
                                        ('strand', 'S1'),
                                        ('event_index', '<i8'), ('event_mean', '<f8'), ('event_noise', '<f8'),
                                        ('event_duration', '<f8'), ('aligned_kmer', 'S5'),
                                        ('scaled_mean_current', '<f8'), ('scaled_noise', '<f8'),
                                        ('posterior_probability', '<f8'), ('descaled_event_mean', '<f8'),
                                        ('ont_model_mean', '<f8'), ('path_kmer', 'S5')])
            get_mea_params_from_events(events)

    def test_match_events_with_signalalign(self):
        # """Test match_events_with_signalalign"""
        # RNA is sequenced 3'-5'
        # reversed for fasta/q sequence
        # if mapped to reverse strand
        # reverse reverse complement = complement
        # DNA is sequenced 5'-3'
        # if mapped to reverse strand
        # reverse complement

        alignment = np.zeros(4, dtype=[('reference_index', int), ('event_index', int), ('posterior_probability', float),
                                       ('reference_kmer', 'S5')])

        alignment["reference_index"] = [0, 1, 2, 3]
        alignment["event_index"] = [0, 1, 2, 3]
        alignment["posterior_probability"] = [1, 1, 1, 1]
        alignment["reference_kmer"] = ["AAGG", "AGGC", "AGGC", "GGCT"]

        event_detects = np.zeros(4, dtype=[('raw_start', int), ('raw_length', int)])
        event_detects["raw_start"] = [10, 11, 12, 13]
        event_detects["raw_length"] = [1, 1, 1, 1]

        labels = match_events_with_signalalign(sa_events=alignment,
                                               event_detections=event_detects,
                                               minus=False,
                                               rna=False)
        # DNA, plus strand
        self.assertSequenceEqual([bytes.decode(x) for x in labels["kmer"]], ["AAGG", "AGGC", "AGGC", "GGCT"])
        self.assertSequenceEqual(labels["raw_start"].tolist(), [10, 11, 12, 13])

        # DNA, minus strand
        labels = match_events_with_signalalign(sa_events=alignment,
                                               event_detections=event_detects,
                                               minus=True,
                                               rna=False)

        self.assertSequenceEqual([bytes.decode(x) for x in labels["kmer"]], ["CCTT", "GCCT", "GCCT", "AGCC"])
        self.assertSequenceEqual(labels["raw_start"].tolist(), [10, 11, 12, 13])

        labels = match_events_with_signalalign(sa_events=alignment,
                                               event_detections=event_detects,
                                               minus=False,
                                               rna=True)
        # RNA, plus strand
        self.assertSequenceEqual([bytes.decode(x) for x in labels["kmer"]], ["GGAA", "CGGA", "CGGA", "TCGG"])
        self.assertSequenceEqual(labels["raw_start"].tolist(), [10, 11, 12, 13])

        labels = match_events_with_signalalign(sa_events=alignment,
                                               event_detections=event_detects,
                                               minus=True,
                                               rna=True)
        # RNA, minus strand # just complement
        self.assertSequenceEqual([bytes.decode(x) for x in labels["kmer"]], ["TTCC", "TCCG", "TCCG", "CCGA"])
        self.assertSequenceEqual(labels["raw_start"].tolist(), [10, 11, 12, 13])
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
