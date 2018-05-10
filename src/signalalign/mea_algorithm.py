#!/usr/bin/env python
"""Create maximum expected accuracy alignment algorithm for signal align best alignment path"""
########################################################################
# File: mea_algorithm.py
#  executable: mea_algorithm.py
#
# Author: Andrew Bailey
# History: 1/17/18 Created
########################################################################


import sys
import os
import numpy as np
from scipy import sparse
from timeit import default_timer as timer
from signalalign.fast5 import Fast5
from py3helpers.utils import list_dir, check_numpy_table
from py3helpers.seq_tools import ReverseComplement
from collections import defaultdict
import traceback


def maximum_expected_accuracy_alignment(posterior_matrix, shortest_ref_per_event, return_all=False,
                                        sparse_posterior_matrix=None):
    """Computes the maximum expected accuracy alignment along a reference with given events and probabilities

    :param posterior_matrix: matrix of posterior probabilities with reference as columns and
                            events as rows.
    :param shortest_ref_per_event: list of the highest possible reference position for all future events at a
                                    given index
    :param return_all: option to return all the paths through the matrix
    :param sparse_posterior_matrix: bool sparse matrix option

    :return best_path: a nested list of lists-> last_event = [ref_index, event_index, prob, sum_prob, [prev_event]]
    """
    # optional convert to sparse matrix
    if sparse_posterior_matrix:
        sparse_posterior_matrix = posterior_matrix
    else:
        sparse_posterior_matrix = sparse.coo_matrix(posterior_matrix)
    forward_edges = list()
    # get the index of the largest probability for the first event
    smallest_event = min(sparse_posterior_matrix.row)
    first_events = sparse_posterior_matrix.row == smallest_event
    num_first_event = sum(first_events)
    max_prob = 0
    largest_start_index = np.argmax(sparse_posterior_matrix.data[first_events])
    # gather leading edges for all references above max
    for x in range(int(largest_start_index) + 1):
        event_index = sparse_posterior_matrix.row[x]
        posterior = sparse_posterior_matrix.data[x]
        ref_index = sparse_posterior_matrix.col[x]
        event_data = [ref_index, event_index,
                      posterior, posterior, None]
        if posterior >= max_prob:
            forward_edges.append(event_data)
            max_prob = posterior
    # skip the inner loop for finding new events by setting prev_event to the next event
    prev_event = sparse_posterior_matrix.row[num_first_event]
    new_edges = list()
    first_pass = True
    max_prob = 0
    i = 0
    # go through rest of events
    num_events = len(sparse_posterior_matrix.row)
    # print("NEW CALL")
    for j in range(num_first_event, num_events):
        event_index = sparse_posterior_matrix.row[j]
        posterior = sparse_posterior_matrix.data[j]
        ref_index = sparse_posterior_matrix.col[j]
        # print("NEW REF", ref_index, event_index, posterior)
        # print("new edges ----", new_edges)
        # update forward edges if new event
        if prev_event != event_index:
            # print("NEW EVENT - OLD EDGES", forward_edges)
            prev_event = event_index
            # capture edges that are further along than the previous last event
            # [forward_edges.pop(x) for x in pop]
            while i < len(forward_edges):
                # print(i)
                if forward_edges[i][3] > max_prob:
                    new_edges.append(forward_edges[i])
                    max_prob = forward_edges[i][3]
                i += 1
                # print(i)

            # forward_edges = forward_edges[:i]

            # reset edges
            first_pass = True
            forward_edges = new_edges
            # print("NEW FORWARD EDGES", forward_edges)

        if first_pass:
            first_pass = False
            max_i = -1
            new_edges = list()
            i = 0
            max_prob = 0
            # keep edge to capture shortest ref position after current event
            edge_found = False
            # print("FIRST pass", i, forward_edges)
            while forward_edges[i][0] < shortest_ref_per_event[event_index]:
                i += 1
                edge_found = True
                if i == len(forward_edges):
                    break
            # make sure we don't double dip on assigning events for shortest event refs
            if edge_found:
                # print("EDGE FOUND")
                test_edge = forward_edges[i - 1]
                new_edges.append(test_edge)
                max_prob = test_edge[3]

                # forward_edges = forward_edges[i-1:]
                # max_prob = forward_edges[0][3]
                assert max_prob == test_edge[3]
            i = 0
        # else:
        assigning = True
        while assigning:
            # print("Another_loop", i)
            if i < len(forward_edges):
                if forward_edges[i][0] < ref_index:
                    if i > max_i:
                        if max_prob < forward_edges[i][3]:
                            new_edges.append(forward_edges[i])
                            max_prob = forward_edges[i][3]
                            max_i = i
                    i += 1
                elif forward_edges[i][0] == ref_index:
                    if i == 0:
                        if forward_edges[i][3] > max_prob:
                            new_edges.append([ref_index, event_index, posterior, forward_edges[i][3], forward_edges[i]])
                            max_prob = forward_edges[i][3]
                    # compare with prev_edge
                    elif forward_edges[i][3] > forward_edges[i - 1][3] + posterior:
                        if forward_edges[i][3] > max_prob:
                            new_edges.append([ref_index, event_index, posterior, forward_edges[i][3], forward_edges[i]])
                            max_prob = forward_edges[i][3]
                    # assign move
                    elif forward_edges[i - 1][3] + posterior > max_prob:
                        new_edges.append([ref_index, event_index, posterior, forward_edges[i - 1][3] + posterior,
                                          forward_edges[i - 1]])
                        max_prob = forward_edges[i - 1][3] + posterior
                    assigning = False
                    # print("Assigned", new_edges)
                    max_i = i
                else:
                    if i == 0:
                        if posterior > max_prob:
                            new_edges.append([ref_index, event_index, posterior, posterior, None])
                            max_prob = posterior
                    elif forward_edges[i - 1][3] + posterior > max_prob:
                        new_edges.append([ref_index, event_index, posterior,
                                          forward_edges[i - 1][3] + posterior, forward_edges[i - 1]])
                        max_prob = forward_edges[i - 1][3] + posterior
                    assigning = False
                    # print("Assigned", new_edges)

            # ref past all edges
            elif forward_edges[i - 1][3] + posterior > max_prob:
                new_edges.append(
                    [ref_index, event_index, posterior, forward_edges[i - 1][3] + posterior, forward_edges[i - 1]])
                max_prob = forward_edges[i - 1][3] + posterior
                assigning = False
                # print("Assigned", new_edges)
            else:
                # print("Didnt Assign")
                assigning = False
        # print("?????_new_edges", new_edges)

    # catch trailing edges
    while i < len(forward_edges):
        # print(i)
        if forward_edges[i][3] > max_prob:
            new_edges.append(forward_edges[i])
        i += 1

    forward_edges = new_edges
    # print('\n')
    # print("FINAL EDGES", forward_edges)
    # print('\n')

    # grab and return the highest probability edge
    if return_all:
        return forward_edges
    else:
        highest_prob = 0
        best_forward_edge = 0
        for x in forward_edges:
            if x[3] > highest_prob:
                highest_prob = x[3]
                best_forward_edge = x
        return best_forward_edge


def binary_search_for_edge(forward_edges, ref_index, event_index, posterior):
    """Search the forward edges list for best ref index comparison
    :param forward_edges: list of forward edges to search
    :param ref_index: index to match with forward edges list
    :param event_index: information to be passed into new forward edge link
    :param posterior: posterior probability of event for a kmer at ref_index
    :return: new forward edge
    """
    assert forward_edges[0][0] <= ref_index, "Ref index cannot be smaller than smallest forward edge"
    searching = True
    last_index = len(forward_edges) - 1
    r = last_index
    l = 0
    i = ((r + l) // 2)

    while searching:
        # check events above ref index
        edge = forward_edges[i]
        edge_ref = edge[0]
        # print(i, edge, edge_ref, ref_index, l)
        if edge_ref <= ref_index:
            if edge_ref == ref_index:
                # if first edge
                if i == 0:
                    return [ref_index, event_index, posterior, edge[3], edge]
                else:
                    earlier_edge = forward_edges[i - 1]
                    if edge[3] > earlier_edge[3] + posterior:
                        return [ref_index, event_index, posterior, edge[3], edge]
                    else:
                        return [ref_index, event_index, posterior, earlier_edge[3] + posterior, earlier_edge]
            # if before last index
            else:
                if i == last_index:
                    # last event is only one above ref index
                    return [ref_index, event_index, posterior, edge[3] + posterior, edge]
                # check if next edge is after ref index
                elif forward_edges[i + 1][0] > ref_index:
                    return [ref_index, event_index, posterior, edge[3] + posterior, edge]
                # next event is either before ref or equal to it
                else:
                    l = i + 1
                    i = (l + r) // 2
        else:
            r = i - 1
            i = (l + r) // 2


def get_indexes_from_best_path(best_path):
    """Grab the reference and event index of the best path from the maximum_expected_accuracy_alignment function.
    :param best_path: output from maximum_expected_accuracy_alignment
    :return: list of events [[ref_pos, event_pos]...]
    """
    path = []
    while best_path[4]:
        ref_pos = best_path[0]
        event_pos = best_path[1]
        path.append([ref_pos, event_pos])
        best_path = best_path[4]
    # gather last event
    ref_pos = best_path[0]
    event_pos = best_path[1]
    path.append([ref_pos, event_pos])
    # flip ordering of path
    return path[::-1]


def get_mea_params_from_events(events):
    """Get the posterior matrix, shortest_ref_per_event and event matrix from events table

    :param events: events table with required fields"""
    check_numpy_table(events, req_fields=('contig', 'reference_index', 'reference_kmer', 'strand', 'event_index',
                                         'event_mean', 'event_noise', 'event_duration', 'aligned_kmer',
                                         'scaled_mean_current', 'scaled_noise', 'posterior_probability',
                                         'descaled_event_mean', 'ont_model_mean', 'path_kmer'))
    # get min/max args
    ref_start = min(events["reference_index"])
    ref_end = max(events["reference_index"])

    # sort events to collect the min ref position per event
    events = np.sort(events, order=['event_index'], kind='mergesort')
    event_start = events["event_index"][0]
    event_end = events["event_index"][-1]

    # check strand of the read
    minus_strand = False
    if events[0]["reference_index"] > events[-1]["reference_index"]:
        minus_strand = True
    # print("minus_strand", minus_strand)

    ref_length = int(ref_end - ref_start + 1)
    event_length = int(event_end - event_start + 1)

    # initialize data structures
    event_matrix = [[0 for _ in range(ref_length)] for _ in range(event_length)]
    posterior_matrix = np.zeros([event_length, ref_length])
    shortest_ref_per_event = [np.inf for _ in range(event_length)]

    min_shortest_ref = np.inf
    # print(ref_start, ref_end)
    # go through events backward to make sure the shortest ref per event is calculated at the same time
    for i in range(1, len(events) + 1):
        event = events[-i]
        event_indx = event["event_index"] - event_start
        if minus_strand:
            ref_indx = event["reference_index"] - ref_end
            ref_indx *= -1
        else:
            ref_indx = event["reference_index"] - ref_start

        # using full event so we can use the same matrix to assemble the training data later
        posterior_matrix[event_indx][ref_indx] = event['posterior_probability']
        event_matrix[event_indx][ref_indx] = event
        # edit shortest ref per event list
        if shortest_ref_per_event[event_indx] > ref_indx:
            if min_shortest_ref > ref_indx:
                min_shortest_ref = ref_indx
            # print(event_indx, ref_indx, min_shortest_ref)
            shortest_ref_per_event[event_indx] = min_shortest_ref

    return posterior_matrix, shortest_ref_per_event, event_matrix


def mea_alignment_from_signal_align(fast5_path, events=None):
    """Get the maximum expected alignment from a nanopore read fast5 file which has signalalign data

    :param fast5_path: path to fast5 file
    :param events: directly pass events in via a numpy array
    """
    if events is None:
        assert os.path.isfile(fast5_path)
        fileh = Fast5(fast5_path)
        events = fileh.get_signalalign_events()

    posterior_matrix, shortest_ref_per_event, event_matrix = get_mea_params_from_events(events)
    # get mea alignment
    mea_alignments = maximum_expected_accuracy_alignment(posterior_matrix, shortest_ref_per_event)
    # get raw index values from alignment data structure
    best_path = get_indexes_from_best_path(mea_alignments)
    # corrected_path = fix_path_indexes(best_path)
    final_event_table = get_events_from_path(event_matrix, best_path)
    return final_event_table


def get_events_from_path(event_matrix, path):
    """Return an event table from a list of index pairs generated from the mea alignment

    :param event_matrix: matrix [ref x event] with event info at positions in matrix
    :param path: [[ref_pos, event_pos]...] to gather events
    """
    events = np.zeros(0, dtype=[('contig', 'S10'), ('reference_index', '<i8'), ('reference_kmer', 'S5'),
                                ('strand', 'S1'),
                                ('event_index', '<i8'), ('event_mean', '<f8'), ('event_noise', '<f8'),
                                ('event_duration', '<f8'), ('aligned_kmer', 'S5'),
                                ('scaled_mean_current', '<f8'), ('scaled_noise', '<f8'),
                                ('posterior_probability', '<f8'), ('descaled_event_mean', '<f8'),
                                ('ont_model_mean', '<f8'), ('path_kmer', 'S5')])
    events_dtype = events.dtype
    # for each pair, access event info from matrix
    for index_pair in path:
        ref_pos = index_pair[0]
        event_pos = index_pair[1]
        try:
            events = np.append(events, np.array(event_matrix[event_pos][ref_pos], dtype=events_dtype))
        except IndexError:
            traceback.print_exc(file=sys.stderr)
            raise IndexError("Selected non event location in event matrix. Check path for correct indexes")
    return events


def match_events_with_signalalign(sa_events=None, event_detections=None, minus=False, rna=False):
    """Match event index with event detection data to label segments of signal for each kmer

    # RNA is sequenced 3'-5'
    # reversed for fasta/q sequence
    # if mapped to reverse strand
    # reverse reverse complement = complement

    # DNA is sequenced 5'-3'
    # if mapped to reverse strand
    # reverse complement

    :param sa_events: events table reference_index', 'event_index', 'aligned_kmer', 'posterior_probability
    :param event_detections: event detection event table
    :param minus: boolean option to for minus strand mapping
    :param rna: boolean for RNA read
    """
    assert sa_events is not None, "Must pass signal alignment events"
    assert event_detections is not None, "Must pass event_detections events"

    check_numpy_table(sa_events, req_fields=('reference_index', 'event_index',
                                            'reference_kmer', 'posterior_probability'))

    check_numpy_table(event_detections, req_fields=('raw_start', 'raw_length'))

    label = np.zeros(len(sa_events), dtype=[('raw_start', int), ('raw_length', int), ('reference_index', int),
                                            ('posterior_probability', float), ('kmer', 'S5')])

    label['raw_start'] = [event_detections[x]["raw_start"] for x in sa_events["event_index"]]
    label['raw_length'] = [event_detections[x]["raw_length"] for x in sa_events["event_index"]]
    label['reference_index'] = sa_events["reference_index"]

    def convert_to_str(string):
        """Helper function to catch bytes as strings"""
        if type(string) is str:
            return string
        else:
            return bytes.decode(string)

    flip = ReverseComplement()
    if minus:
        if rna:
            kmers = [flip.complement(convert_to_str(x)) for x in sa_events["reference_kmer"]]
        else:
            kmers = [flip.reverse_complement(convert_to_str(x)) for x in sa_events["reference_kmer"]]
    else:
        if rna:
            kmers = [flip.reverse(convert_to_str(x)) for x in sa_events["reference_kmer"]]
        else:
            kmers = sa_events["reference_kmer"]
    label['kmer'] = kmers
    label['posterior_probability'] = sa_events["posterior_probability"]
    np.sort(label, order='raw_start', kind='mergesort')

    return label


def slow_search_for_edge(forward_edges, ref_index, event_index, posterior):
    """Search the forward edges list for best ref index comparison
    :param forward_edges: list of forward edges to search
    :param ref_index: index to match with forward edges list
    :param event_index: information to be passed into new forward edge link
    :param posterior: posterior probability of event for a kmer at ref_index
    :return: new forward edge
    """
    assert forward_edges[0][0] <= ref_index, "Ref index cannot be smaller than smallest forward edge"

    inxs = []
    probs = []
    for j, forward_edge in enumerate(forward_edges):
        if forward_edge[0] < ref_index:
            # track which probabilities with prev edge
            inxs.append(j)
            probs.append(posterior + forward_edge[3])
        elif forward_edge[0] == ref_index:
            # stay at reference position
            # add probability of event if we want to promote sideways movement
            inxs.append(j)
            probs.append(forward_edge[3])
    # add most probable connecting edge if better than creating an new edge
    # deal with multiple edges with equal probabilty
    probs = probs[::-1]
    inxs = inxs[::-1]
    connecting_edge = forward_edges[inxs[int(np.argmax(probs))]]
    return [ref_index, event_index, posterior, max(probs), connecting_edge]


def sum_forward_edge_accuracy(forward_edge):
    """Add up all probabilities from a forward edge"""
    sum_prob = 0

    def get_probability_info(forward_edge):
        """Get event information from maximum_expected_accuracy_alignment"""
        nonlocal sum_prob
        sum_prob = sum_prob + forward_edge[2]
        if forward_edge[4] is not None:
            get_probability_info(forward_edge[4])
        else:
            pass

    get_probability_info(forward_edge)

    return sum_prob


def create_random_prob_matrix(row=None, col=None, gaps=True):
    """Create a matrix of random probability distributions along each row

    :param row: number of rows
    :param col: number of columns
    :param gaps: if have start gap, middle gaps or end gap
    """
    assert row is not None, "Must set row option"
    assert col is not None, "Must set col option"
    prob_matrix = np.zeros([row, col])
    shortest_future_col_per_row = np.zeros(row)
    shortest_col = col
    start = np.random.randint(0, 3)
    skip = 0
    if not gaps:
        start = 0
        skip = 1
    col_indexes = [x for x in range(col)]
    for row_i in range(row - 1, start - 1, -1):
        a = np.random.random(np.random.randint(skip, col))
        a /= a.sum()
        # go through events backward to make sure the shortest ref per event is calculated at the same time
        a = np.sort(a)
        np.random.shuffle(col_indexes)
        # make sure we dont have gaps at ends of columns
        if not gaps and row_i == row - 1:
            col_indexes.remove(col - 1)
            col_indexes.insert(0, col - 1)
        if not gaps and row_i == 0:
            col_indexes.remove(0)
            col_indexes.insert(0, 0)

        for prob, col_i in zip(a, col_indexes):
            # using full event so we can use the same matrix to assemble the training data later
            prob_matrix[row_i][col_i] = prob
            if col_i <= shortest_col:
                shortest_col = col_i
        shortest_future_col_per_row[row_i] = shortest_col
    # check with naive NxM algorithm
    prob_matrix = np.asanyarray(prob_matrix)
    assert matrix_event_length_pairs_test(prob_matrix, shortest_future_col_per_row), \
        "Did not create accurate prob matrix and shortest_future_col_per_row"
    return prob_matrix, shortest_future_col_per_row


def matrix_event_length_pairs_test(posterior_matrix, shortest_events):
    """Test if the shortest events list matches what is in the posterior matrix"""
    # posterior matrix is events x ref
    current_min = np.inf
    shortest_events = shortest_events[::-1]
    for i, row in enumerate(posterior_matrix[::-1]):
        if sum(row) > 0:
            min_event_in_row = min(np.nonzero(row)[0])
            if min_event_in_row < current_min:
                current_min = min_event_in_row
            if shortest_events[i] != current_min:
                return False
    return True


def generate_events_from_probability_matrix(matrix):
    """Create events from probability matrix for testing get_mea_params_from_events"""
    event_indexs = []
    probs = []
    ref_indexs = []

    for row_index, row in enumerate(matrix):
        for col_index, prob in enumerate(row):
            if prob != 0:
                probs.append(prob)
                event_indexs.append(row_index)
                ref_indexs.append(col_index)

    ref_indexs = np.asanyarray(ref_indexs)
    event_indexs = np.asanyarray(event_indexs)
    # create events table
    n_events = len(probs)
    events = np.zeros(n_events, dtype=[('contig', 'S10'), ('reference_index', '<i8'), ('reference_kmer', 'S5'),
                                       ('strand', 'S1'),
                                       ('event_index', '<i8'), ('event_mean', '<f8'), ('event_noise', '<f8'),
                                       ('event_duration', '<f8'), ('aligned_kmer', 'S5'),
                                       ('scaled_mean_current', '<f8'), ('scaled_noise', '<f8'),
                                       ('posterior_probability', '<f8'), ('descaled_event_mean', '<f8'),
                                       ('ont_model_mean', '<f8'), ('path_kmer', 'S5')])

    #  add to ref and event starts
    ref_start = np.random.randint(0, 10000)
    event_start = np.random.randint(0, 10000)

    events["reference_index"] = ref_indexs + ref_start
    events["posterior_probability"] = probs
    events["event_index"] = event_indexs + event_start

    # create comparison event matrix
    event_matrix = np.zeros(matrix.shape).tolist()

    for i in range(n_events):
        event_matrix[event_indexs[i]][ref_indexs[i]] = events[i]
        assert events[i]["reference_index"] == ref_indexs[i] + ref_start
        assert events[i]["event_index"] == event_indexs[i] + event_start

    return events, event_matrix


def mea_slower(posterior_matrix, shortest_ref_per_event, return_all=False,
               sparse_posterior_matrix=None):
    """Computes the maximum expected accuracy alignment along a reference with given events and probabilities

    NOTE: Slower than other version

    :param posterior_matrix: matrix of posterior probabilities with reference along x axis (col) and
                            events along y axis (row).
    :param shortest_ref_per_event: list of the highest possible reference position for all future events at a
                                    given index
    :param return_all: option to return all the paths through the matrix
    :param sparse_posterior_matrix: bool sparse matrix option

    :return best_path: a nested list of lists-> last_event = [ref_index, event_index, prob, sum_prob, [prev_event]]
    """
    # optional convert to sparse matrix
    if sparse_posterior_matrix:
        sparse_posterior_matrix = posterior_matrix
    else:
        sparse_posterior_matrix = sparse.coo_matrix(posterior_matrix)

    forward_edges = list()
    # get the index of the largest probability for the first event
    smallest_event = min(sparse_posterior_matrix.row)
    first_events = sparse_posterior_matrix.row == smallest_event
    num_first_event = sum(first_events)

    largest_start_prob = np.argmax(sparse_posterior_matrix.data[first_events])
    # gather leading edges for all references above max
    for x in range(int(largest_start_prob) + 1):
        event_data = [sparse_posterior_matrix.col[x], sparse_posterior_matrix.row[x],
                      sparse_posterior_matrix.data[x], sparse_posterior_matrix.data[x], None]
        forward_edges.append(event_data)
    # number of values for first event
    prev_event = sparse_posterior_matrix.row[num_first_event]
    new_edges = list()
    first_pass = True
    prev_ref_pos = 0
    fill_gap = False
    # go through rest of events
    num_events = len(sparse_posterior_matrix.row)
    for i in range(num_first_event, num_events):
        event_index = sparse_posterior_matrix.row[i]
        posterior = sparse_posterior_matrix.data[i]
        ref_index = sparse_posterior_matrix.col[i]
        # update forward edges if new event
        if prev_event != event_index:
            prev_event = event_index
            # capture edges that are further along than the current event
            for forward_edge in forward_edges:
                if forward_edge[0] > prev_ref_pos:
                    new_edges.append(forward_edge)
            forward_edges = new_edges
            new_edges = list()
            first_pass = True
        # check if there is a gap between reference
        if prev_ref_pos + 1 != ref_index and not first_pass:
            fill_gap = True
            gap_indicies = [x for x in range(prev_ref_pos + 1, ref_index)]
        # keep track of probabilities to select best connecting edge to new node
        inxs = []
        probs = []
        # event_data = [ref_index, event_index, prob, sum_prob, None]
        for j, forward_edge in enumerate(forward_edges):
            if forward_edge[0] < ref_index:
                # track which probabilities with prev edge
                inxs.append(j)
                probs.append(posterior + forward_edge[3])
                # if needed, keep edges aligned to ref positions previous than the current ref position
                if first_pass and shortest_ref_per_event[event_index] < forward_edge[0] + 2:
                    new_edges.append(forward_edge)
            elif forward_edge[0] == ref_index:
                # stay at reference position
                # add probability of event if we want to promote sideways movement
                inxs.append(j)
                probs.append(forward_edge[3])
            if fill_gap:
                if forward_edge[0] in gap_indicies:
                    # add edges that pass through gaps in the called events
                    new_edges.append(forward_edge)
        # add most probable connecting edge if better than creating an new edge
        if probs:
            connecting_edge = forward_edges[inxs[int(np.argmax(probs))]]
            new_edges.append([ref_index, event_index, posterior, max(probs), connecting_edge])
        else:
            # no possible connecting edges or connecting edges decrease probability, create a new one
            new_edges.append([ref_index, event_index, posterior, posterior, None])

        # reset trackers
        first_pass = False
        prev_ref_pos = ref_index
        fill_gap = False

    # add back last edges which may not have been connected
    for forward_edge in forward_edges:
        if forward_edge[0] > prev_ref_pos:
            new_edges.append(forward_edge)
    forward_edges = new_edges
    # grab and return the highest probability edge
    if return_all:
        return forward_edges
    else:
        highest_prob = 0
        best_forward_edge = 0
        for x in forward_edges:
            if x[3] > highest_prob:
                highest_prob = x[3]
                best_forward_edge = x
        return best_forward_edge


def mea_slow(posterior_matrix, shortest_ref_per_event, return_all=False):
    """Computes the maximum expected accuracy alignment along a reference with given events and probabilities.

    Computes a very slow but thorough search through the matrix

    :param posterior_matrix: matrix of posterior probabilities with reference along x axis and events along y
    :param shortest_ref_per_event: shortest ref position per event
    :param return_all: return all forward edges
    """
    ref_len = len(posterior_matrix[0])
    events_len = len(posterior_matrix)
    initialize = True
    forward_edges = list()
    new_edges = list()
    # step through all events
    for event_index in range(events_len):
        max_prob = 0
        if initialize:
            ref_index = 0
            while ref_index < ref_len:
                # intitialize forward edges with first event alignments
                # if type(posterior_matrix[ref_index][event_index]) is not int:
                posterior = posterior_matrix[event_index][ref_index]
                event_data = [ref_index, event_index, posterior, posterior, None]
                if 0 < posterior >= max_prob:
                    # print("True", posterior, max_prob)
                    new_edges.append(event_data)
                    max_prob = posterior
                ref_index += 1
            # print("INITIALIZE", new_edges, max_prob)
            if len(new_edges) != 0:
                forward_edges = new_edges
                new_edges = list()
                initialize = False
        else:
            # print(forward_edges)
            ref_index = 0
            top_edge = []
            while ref_index < ref_len:
                posterior = posterior_matrix[event_index][ref_index]
                if posterior >= max_prob:
                    # no possible connecting edges and is needed for other other events create a new one
                    if ref_index < shortest_ref_per_event[event_index]:
                        top_edge.append([ref_index, event_index, posterior, posterior, None])
                        max_prob = posterior
                ref_index += 1
            # add top edge if needed
            if top_edge:
                new_edges.append(top_edge[-1])
            ref_index = 0
            while ref_index < ref_len:
                inxs = []
                probs = []
                posterior = posterior_matrix[event_index][ref_index]
                for j, forward_edge in enumerate(forward_edges):
                    if forward_edge[0] < ref_index:
                        # track which probabilities with prev edge
                        inxs.append(j)
                        probs.append(posterior + forward_edge[3])
                        # if needed, keep edges aligned to ref positions previous than the current ref position
                    elif forward_edge[0] == ref_index:
                        # stay at reference position
                        # add probability of event if we want to promote sideways movement
                        inxs.append(j)
                        probs.append(forward_edge[3])

                # add new edge
                inxs = inxs[::-1]
                probs = probs[::-1]
                if len(probs) != 0:
                    if max(probs) > max_prob:
                        connecting_edge = forward_edges[inxs[int(np.argmax(probs))]]
                        new_edges.append([ref_index, event_index, posterior, max(probs), connecting_edge])
                        max_prob = max(probs)
                else:
                    if forward_edges[0][0] > ref_index and posterior > max_prob:
                        new_edges.append([ref_index, event_index, posterior, posterior, None])
                        max_prob = posterior
                ref_index += 1
            # print("END_NEW_EDGES", new_edges)
            forward_edges = new_edges
            new_edges = list()
    # grab and return the highest probability edge
    if return_all:
        return forward_edges
    else:
        highest_prob = 0
        best_forward_edge = 0
        for x in forward_edges:
            if x[3] > highest_prob:
                highest_prob = x[3]
                best_forward_edge = x
        return best_forward_edge


def main():
    """Main docstring"""
    start1 = timer()

    # m, l = create_random_prob_matrix(row=3, col=3)
    # print(m.shape)

    # posterior_matrix = [[0.2, 0.3, 0.2, 0.2, 0.1],
    #                     [0.2, 0.5, 0.3, 0.0, 0.0],
    #                     [0.3, 0.1, 0.0, 0.3, 0.3],
    #                     [0.0, 0.0, 0.0, 0.4, 0.1],
    #                     [0.0, 0.0, 0.0, 0.2, 0.5],
    #                     ]
    #
    # # correct input
    # shortest_ref_per_event = [0, 0, 0, 3, 3]
    # # best_edge = binary_search_for_edge([[0, 1, .1, .1], [1, 1, .1, .1], [2, 1, .1, .1], [3, 1, .1, .1], [4, 1, .1, .1], [5, 1, .1, .1], [6, 1, .1, .1]], 6.1, 2, 0.1)
    # # print(best_edge)
    #
    # forward_edges = maximum_expected_accuracy_alignment_edits(posterior_matrix, shortest_ref_per_event, return_all=True)
    # # print(forward_edges)
    # forward_edges2 = maximum_expected_accuracy_alignment(posterior_matrix, shortest_ref_per_event, return_all=True)
    # print("COMPARE")
    #
    # for x, y in zip(forward_edges, forward_edges2):
    #     print(x)
    #     print(y)
    #     print("NEW PAIR")
    # fast5_path = "/Users/andrewbailey/CLionProjects/nanopore-RNN/test_files/minion-reads/canonical"
    # files = list_dir(fast5_path, ext='fast5')
    # for file1 in files:
    #     start = timer()
    #     print(file1)
    #     f5fh = Fast5(file1)
    #     mae_events = f5fh.get_signalalign_events(mae=True)
    #     event_detection = f5fh.get_resegment_basecall()
    #     # events = mea_alignment_from_signal_align(file1)
    #     label = match_events_with_signalalign(mae_events, event_detection)
    #     print(label)
    #     event1 = 0
    #     for event in label:
    #         if type(event1) is int:
    #             event1 = event
    #         else:
    #             if event1["raw_start"]+event1["raw_length"] == event["raw_start"]:
    #                 pass
    #             else:
    #                 print(event1["raw_start"]+event1["raw_length"])
    #                 print(event["raw_start"])
    #                 print("error")
    #             event1 = event
    #     # break
    #     stop = timer()
    #     print("Total Time = {} seconds".format(stop - start), file=sys.stderr)

    stop = timer()
    print("Running Time = {} seconds".format(stop - start1), file=sys.stderr)


if __name__ == "__main__":
    main()
    raise SystemExit
