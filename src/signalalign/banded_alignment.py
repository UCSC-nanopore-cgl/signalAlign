#!/usr/bin/env python3
"""Segment raw signal and create banded alignment to fasta sequence"""
########################################################################
# File: banded_alignment.py
#  executable: banded_alignment.py
#
# Author: Andrew Bailey
# History: 5/30/18 Created
########################################################################

import os
import pandas as pd
import numpy as np
from collections import namedtuple

from py3helpers.utils import check_numpy_table, TimeStamp, merge_dicts
from signalalign.fast5 import Fast5
from signalalign.event_detection import create_minknow_events_from_fast5
from signalalign.train.trainModels import get_model
from signalalign.hiddenMarkovModel import SignalHmm


def simple_banded_event_align(events, model, nucleotide_seq):
    """Generate a banded alignment between events and a nucleotide sequence

    :param events: event table with required fields: ('start', 'length', 'mean', 'stdv', 'model_state', 'move', 'p_model_state')
    :param model: SignalHmm model
    :param nucleotide_seq: nucleotide sequence to match up
    """
    check_numpy_table(events, req_fields=('start', 'length', 'mean', 'stdv', 'model_state', 'move', 'p_model_state'))
    assert isinstance(model, SignalHmm), "Input model needs to be SignalHmm"
    k = model.kmer_length

    FROM_D = 0
    FROM_U = 1
    FROM_L = 2

    # # // qc
    min_average_log_emission = -5.0

    n_kmers = len(nucleotide_seq) - k + 1
    # // banding
    bandwidth = min(1000, n_kmers+1)
    half_band = bandwidth / 2

    # // transitions
    lp_skip = np.log(0.001)
    lp_stay = np.log(0.5)
    lp_step = np.log(1.0 - np.exp(lp_skip) - np.exp(lp_stay))
    lp_trim = np.log(0.1)

    n_events = len(events)
    events_per_kmer = float(n_kmers) / n_events
    min_event_idx_by_kmer = []
    # // Calculate the minimum event index that is within the band for each read kmer
    # // We determine this using the expected number of events observed per kmer
    for ki in range(n_kmers):
        expected_event_idx = ki * events_per_kmer
        min_event_idx_by_kmer.append(max(int(expected_event_idx - half_band), 0))

    n_rows = bandwidth
    n_cols = n_kmers + 1
    print(min_event_idx_by_kmer)
    print("n_rows", n_rows)
    print("n_cols", n_cols)
    print("seq_len", len(nucleotide_seq))

    viterbi_matrix = np.zeros((n_rows, n_cols), dtype=np.float64)
    backtrack_matrix = np.zeros((n_rows, n_cols), dtype=np.int32)

    for i in range(n_cols):
        viterbi_matrix[0][i] = -np.infty
        # backtrack_matrix[0][i] = 0

    for i in range(bandwidth):
        viterbi_matrix[i][0] = i * lp_trim
        # backtrack_matrix[i][0] = 0

    fills = 0
    for col in range(1, n_cols):
        # print("New_col {}/{}".format(col, n_cols))
        kmer_idx = col - 1
        min_event_idx = min_event_idx_by_kmer[kmer_idx]
        if kmer_idx > 0:
            min_event_idx_prev_col = min_event_idx_by_kmer[kmer_idx - 1]
        else:
            min_event_idx_prev_col = 0
        # kmer_rank = alphabet->kmer_rank(sequence.substr(kmer_idx, k).c_str(), k);

        for row in range(n_rows):
            event_idx = min_event_idx + row
            if event_idx >= n_events:
                viterbi_matrix[row][col] = -np.infty
                continue
            # // dp update
            # // here we are calculating whether the event for each neighboring cell is within the band
            # // and calculating its position within the column
            row_up = event_idx - min_event_idx - 1
            row_diag = event_idx - min_event_idx_prev_col - 1
            row_left = event_idx - min_event_idx_prev_col

            # try:
            if 0 <= row_up < n_rows:
                up = viterbi_matrix[row_up][col]
            else:
                up = -np.infty
            if 0 <= row_diag < n_rows:
                diag = viterbi_matrix[row_diag][col-1]
            else:
                diag = -np.infty
            if 0 <= row_left < n_rows:
                left = viterbi_matrix[row_left][col-1]
            else:
                left = -np.infty
            # except IndexError:
            #     print(row_up, row_diag, row_left, col)
            #     SystemExit
            # lp_emission = log_probability_match_r9(read, pore_model, kmer_rank, event_idx, strand_idx);
            event_mean = events['mean'][event_idx]
            # print(event_idx, event_mean)
            kmer = nucleotide_seq[kmer_idx:kmer_idx+k]
            lp_emission = model.log_event_mean_gaussian_probability_match(event_mean, kmer)

            score_d = diag + lp_step + lp_emission
            score_u = up + lp_stay + lp_emission
            score_l = left + lp_skip

            max_score = score_d
            from_where = FROM_D

            if score_u > max_score:
                max_score = score_u

            if max_score == score_u:
                from_where = FROM_U

            if score_l > max_score:
                max_score = score_l

            if max_score == score_l:
                from_where = FROM_L

            # //fprintf(stderr, "[orgfill] up: %.2lf diag: %.2lf left: %.2lf\n", up, diag, left)
            # fprintf(stderr, "[orgfill] e: %d k: %d s: %.2lf f: %d emit: %.2lf\n", event_idx, kmer_idx, max_score, from, lp_emission)
            viterbi_matrix[row][col] = max_score
            backtrack_matrix[row][col] = from_where
            fills += 1


    # // Initialize by finding best alignment between an event and the last kmer
    curr_k_idx = n_kmers - 1
    curr_event_idx = 0
    max_score = -np.infty
    for row in range(n_rows):
        col = curr_k_idx + 1
        ei = row + min_event_idx_by_kmer[curr_k_idx]
        s = viterbi_matrix[row][col] + (n_events - ei - 1) * lp_trim
        if s > max_score and ei < n_events:
            max_score = s
            curr_event_idx = ei

    sum_emission = 0
    n_aligned_events = 0
    moves = 0
    out = []
    while curr_k_idx >= 0:
        # // emit alignment
        out.append((curr_k_idx, curr_event_idx))
        event_mean = events['mean'][curr_event_idx]
        # print(event_idx, event_mean)
        kmer = nucleotide_seq[curr_k_idx:curr_k_idx+k]

        kmer_emission = model.log_event_mean_gaussian_probability_match(event_mean, kmer)
        sum_emission += kmer_emission
        if moves == 0:
            events['model_state'][curr_event_idx] = kmer
            events['p_model_state'][curr_event_idx] = np.exp(kmer_emission)

        # kmer_rank = alphabet->kmer_rank(sequence.substr(curr_k_idx, k).c_str(), k);
        # sum_emission += log_probability_match_r9(read, pore_model, kmer_rank, curr_event_idx, strand_idx);
        n_aligned_events += 1
        # // update indices using backtrack pointers
        row = curr_event_idx - min_event_idx_by_kmer[curr_k_idx]
        col = curr_k_idx + 1

        from_where = backtrack_matrix[row][col]
        if from_where == FROM_D:
            moves += 1
            events['move'][curr_event_idx] = moves

            curr_k_idx -= 1
            curr_event_idx -= 1
            moves = 0
        elif from_where == FROM_U:
            events['move'][curr_event_idx] = 0
            print("Moves when skip {}".format(moves))
            curr_event_idx -= 1
            moves = 0
        else:
            moves += 1
            curr_k_idx -= 1

    events['move'][0] = 0
    out = out[::-1]
    avg_log_emission = sum_emission / n_aligned_events
    spanned = out[0][0] == 0 and out[-1][0] == n_kmers - 1
    if avg_log_emission < min_average_log_emission or not spanned:
        print(spanned, avg_log_emission)
        print("We failed... fuck")
    return events[out[0][1]:out[-1][1]+1], sum_emission


def adaptive_banded_simple_event_align(events, model, nucleotide_seq, debug=False):
    """Generate a banded alignment between events and a nucleotide sequence using the adapted banded approach
    source: https://www.biorxiv.org/content/biorxiv/early/2017/04/25/130633.full.pdf /  nanopolish

    :param events: event table with required fields: ('start', 'length', 'mean', 'stdv', 'model_state', 'move', 'p_model_state')
    :param model: SignalHmm model
    :param nucleotide_seq: nucleotide sequence to match up
    :param debug: boolean debug option
    """
    # initialize helper functions
    def move_down(curr_band):
        """Move location in dp matrix down one aka add one to event_index"""
        return EventKmerPair(event_idx=curr_band.event_idx + 1, kmer_idx=curr_band.kmer_idx)

    def move_right(curr_band):
        """Move location in dp matrix right one aka add one to kmer_index"""
        return EventKmerPair(event_idx=curr_band.event_idx, kmer_idx=curr_band.kmer_idx + 1)

    def event_kmer_to_band(ei, ki):
        return (ei + 1) + (ki + 1)

    def band_event_to_offset(bi, ei):
        """Get event index of a specific event pair in the 'band_lower_left' pairs and subtract the event offset"""
        return band_lower_left[bi].event_idx - ei

    def band_kmer_to_offset(bi, ki):
        """Subtract kmer index from the specific event pair in the 'band_lower_left' pairs """
        return ki - band_lower_left[bi].kmer_idx

    def is_offset_valid(offset1):
        """Check if the offset is greater than zero and smaller than the bandwidth"""
        return 0 <= offset1 < bandwidth

    def event_at_offset(bi, offset1):
        """Get kmer index minus offset for a band within the 'band_lower_left' array of event/kmer pairs"""
        return band_lower_left[bi].event_idx - offset1

    def kmer_at_offset(bi, offset1):
        """Get kmer index plus offset for a band within the 'band_lower_left' array of event/kmer pairs"""
        return band_lower_left[bi].kmer_idx + offset1

    check_numpy_table(events, req_fields=('start', 'length', 'mean', 'stdv', 'model_state', 'move', 'p_model_state'))
    assert isinstance(model, SignalHmm), "Input model needs to be SignalHmm"
    k = model.kmer_length
    # strand_idx = 0
    # how to deal with 2d?
    n_events = len(events)
    n_kmers = len(nucleotide_seq) - k + 1

    # backtrack markers
    FROM_D = 0
    FROM_U = 1
    FROM_L = 2

    # # // qc
    min_average_log_emission = -5.0
    max_gap_threshold = 50

    # // banding
    bandwidth = 100
    half_bandwidth = int(bandwidth / 2)
    # setting a tiny skip penalty helps keep the true alignment within the adaptive band
    # this was empirically determined (From Nanopolish)

    # transition penalties
    events_per_kmer = float(n_kmers) / n_events
    p_stay = 1 - (1 / (events_per_kmer + 1))

    # transitions
    epsilon = 1e-10
    lp_skip = np.log(epsilon)
    lp_stay = np.log(p_stay)
    lp_step = np.log(1.0 - np.exp(lp_skip) - np.exp(lp_stay))
    lp_trim = np.log(0.01)

    n_events = len(events)
    n_kmers = len(nucleotide_seq) - k + 1

    n_rows = n_events + 1
    n_cols = n_kmers + 1
    n_bands = n_rows + n_cols

    bands = pd.DataFrame(np.zeros([bandwidth, n_bands])) * -np.infty
    trace = pd.DataFrame(np.zeros([bandwidth, n_bands]))

    # Keep track of the event/kmer index for the lower left corner of the band
    # these indices are updated at every iteration to perform the adaptive banding
    # Only the first two bands have their coordinates initialized, the rest are computed adaptively
    EventKmerPair = namedtuple('EventKmerPair', ['event_idx', 'kmer_idx'])
    band_lower_left = [namedtuple('EventKmerPair', ['event_idx', 'kmer_idx']) for _ in range(n_bands)]

    # initialize range of first two bands
    band_lower_left[0].event_idx = half_bandwidth - 1
    band_lower_left[0].kmer_idx = -1 - half_bandwidth
    band_lower_left[1] = move_down(band_lower_left[0])

    # band 0: score zero in the central cell
    start_cell_offset = band_kmer_to_offset(0, -1)
    assert(is_offset_valid(start_cell_offset)), "Offset is outside the bounds [0, {}]: {}".format(bandwidth, start_cell_offset)
    assert(band_event_to_offset(0, -1) == start_cell_offset), "Event offset is not correct:" \
                                                              " {} != {}".format(band_event_to_offset(0, -1),
                                                                                 start_cell_offset)
    bands[0][start_cell_offset] = 0.0

    # band 1: first event is trimmed
    first_trim_offset = band_event_to_offset(1, 0)
    negative_one = kmer_at_offset(1, first_trim_offset)
    assert(negative_one == -1), "Kmer offset is not correct: {} != {}".format(negative_one, -1)
    assert(is_offset_valid(start_cell_offset)), "Offset is outside the bounds [0, {}]: {}".format(bandwidth, offset)
    bands[1][first_trim_offset] = lp_trim
    trace[1][first_trim_offset] = FROM_U

    fills = 0
    # fill in remaining bands
    for band_idx in range(2, n_bands):
        print(band_idx)
        # Determine placement of this band according to Suzuki's adaptive algorithm
        # When both ll and ur are out-of-band (ob) we alternate movements
        # otherwise we decide based on scores

        ll = bands[band_idx - 1][0]
        ur = bands[band_idx - 1][bandwidth - 1]
        ll_ob = ll == -np.infty
        ur_ob = ur == -np.infty

        if ll_ob and ur_ob:
            right = band_idx % 2 == 1
        else:
            right = ll < ur  # Suzuki's rule

        if right:
            band_lower_left[band_idx] = move_right(band_lower_left[band_idx - 1])
        else:
            band_lower_left[band_idx] = move_down(band_lower_left[band_idx - 1])

        # If the trim state is within the band, fill it in here
        trim_offset = band_kmer_to_offset(band_idx, -1)
        if is_offset_valid(trim_offset):
            event_idx = event_at_offset(band_idx, trim_offset)
            if 0 <= event_idx < n_events:
                bands[band_idx][trim_offset] = lp_trim * (event_idx + 1)
                trace[band_idx][trim_offset] = FROM_U
            else:
                bands[band_idx][trim_offset] = -np.infty
        else:
            "This happened!"

        # Get the offsets for the first and last event and kmer
        # We restrict the inner loop to only these values
        kmer_min_offset = band_kmer_to_offset(band_idx, 0)
        kmer_max_offset = band_kmer_to_offset(band_idx, n_kmers)
        event_min_offset = band_event_to_offset(band_idx, n_events - 1)
        event_max_offset = band_event_to_offset(band_idx, -1)

        min_offset = max(kmer_min_offset, event_min_offset)
        min_offset = max(min_offset, 0)

        max_offset = min(kmer_max_offset, event_max_offset)
        max_offset = min(max_offset, bandwidth)

        for offset in range(min_offset, max_offset):
            event_idx = event_at_offset(band_idx, offset)
            kmer_idx = kmer_at_offset(band_idx, offset)

            # kmer_rank = kmer_ranks[kmer_idx]

            offset_up   = band_event_to_offset(band_idx - 1, event_idx - 1)
            offset_left = band_kmer_to_offset(band_idx - 1, kmer_idx - 1)
            offset_diag = band_kmer_to_offset(band_idx - 2, kmer_idx - 1)

            if debug:
                # verify loop conditions
                assert(0 <= kmer_idx < n_kmers)
                assert(0 <= event_idx < n_events)
                assert(offset_diag == band_event_to_offset(band_idx - 2, event_idx - 1))
                assert(offset_up - offset_left == 1)
                assert(0 <= offset < bandwidth)

            if is_offset_valid(offset_up):
                up = bands[band_idx - 1][offset_up]
            else:
                up = -np.infty
            if is_offset_valid(offset_left):
                left = bands[band_idx - 1][offset_left]
            else:
                left = -np.infty
            if is_offset_valid(offset_diag):
                diag = bands[band_idx - 2][offset_diag]
            else:
                diag = -np.infty

            event_mean = events['mean'][event_idx]
            kmer = nucleotide_seq[kmer_idx:kmer_idx+k]
            lp_emission = model.log_event_mean_gaussian_probability_match(event_mean, kmer)

            score_d = diag + lp_step + lp_emission
            score_u = up + lp_stay + lp_emission
            score_l = left + lp_skip

            max_score = score_d
            from_where = FROM_D

            if score_u > max_score:
                max_score = score_u

            if max_score == score_u:
                from_where = FROM_U

            if score_l > max_score:
                max_score = score_l

            if max_score == score_l:
                from_where = FROM_L

            if debug:
                print("[adafill] offset-up: %d offset-diag: %d offset-left: %d\n", offset_up, offset_diag, offset_left)
                print("[adafill] up: %.2lf diag: %.2lf left: %.2lf\n", up, diag, left)
                print("[adafill] bi: %d o: %d e: %d k: %d s: %.2lf f: %d emit: %.2lf\n", band_idx, offset, event_idx, kmer_idx, max_score, from_where, lp_emission)

            bands[band_idx][offset] = max_score
            trace[band_idx][offset] = from_where
            fills += 1


    # Backtrack to compute alignment
    sum_emission = 0
    n_aligned_events = 0

    max_score = -np.infty
    curr_event_idx = 0
    curr_kmer_idx = n_kmers - 1
    # Find best score between an event and the last k-mer. after trimming the remaining events
    for event_idx in range(n_events):
        band_idx1 = event_kmer_to_band(event_idx, curr_kmer_idx)
        # assert(band_idx < bands.size())
        offset = band_event_to_offset(band_idx1, event_idx)
        if is_offset_valid(offset):
            s = bands[band_idx1][offset] + (n_events - event_idx) * lp_trim
            if s > max_score:
                max_score = s
                curr_event_idx = event_idx

    if debug:
        print("[adaback] ei: %d ki: %d s: %.2f\n", curr_event_idx, curr_kmer_idx, max_score)

    out = []
    curr_gap = 0
    max_gap = 0
    moves = 0
    while curr_kmer_idx >= 0 and curr_event_idx >= 0:
        # emit alignment
        out.append((curr_kmer_idx, curr_event_idx))
        if debug:
            print("[adaback] ei: %d ki: %d\n", curr_event_idx, curr_kmer_idx)
        # qc stats
        event_mean = events['mean'][curr_event_idx]
        # print(event_idx, event_mean)
        kmer = nucleotide_seq[curr_kmer_idx:curr_kmer_idx+k]

        kmer_emission = model.log_event_mean_gaussian_probability_match(event_mean, kmer)
        sum_emission += kmer_emission
        if moves == 0:
            events['model_state'][curr_event_idx] = kmer
            events['p_model_state'][curr_event_idx] = np.exp(kmer_emission)

        n_aligned_events += 1

        band_idx = event_kmer_to_band(curr_event_idx, curr_kmer_idx)
        offset = band_event_to_offset(band_idx, curr_event_idx)
        assert(band_kmer_to_offset(band_idx, curr_kmer_idx) == offset)

        from_where = trace[band_idx][offset]

        if from_where == FROM_D:
            moves += 1
            events['move'][curr_event_idx] = moves
            moves = 0

            curr_kmer_idx -= 1
            curr_event_idx -= 1
            curr_gap = 0

        elif from_where == FROM_U:
            events['move'][curr_event_idx] = 0
            print("Moves when skip {}".format(moves))
            moves = 0
            curr_event_idx -= 1
            curr_gap = 0

        else:
            moves += 1
            curr_kmer_idx -= 1
            curr_gap += 1
            max_gap = max(curr_gap, max_gap)
    events['move'][0] = 0
    # QC results
    out = out[::-1]
    avg_log_emission = sum_emission / n_aligned_events
    spanned = out[0][0] == 0 and out[-1][0] == n_kmers - 1
    if avg_log_emission < min_average_log_emission or not spanned or max_gap > max_gap_threshold:
        print(spanned, avg_log_emission)
        print("We failed... fuck")
    # //fprintf(stderr, "ada\t%s\t%s\t%.2lf\t%zu\t%.2lf\t%d\t%d\t%d\n", read.read_name.substr(0, 6).c_str(), failed ? "FAILED" : "OK", events_per_kmer, sequence.size(), avg_log_emission, curr_event_idx, max_gap, fills);

    return events[out[0][1]:out[-1][1]+1], sum_emission


def main():
    # fast5_path = "/Users/andrewbailey/CLionProjects/nanopolish/ucsc_run5_20170922_directRNA/ucsc_run5_20170922_directRNA_fast5/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_744_ch_111_strand.fast5"
    # nucleotide_seq = "CCUGAAAGCAAUACCUGAUGGAGGCAGCAACAAAGUGUUCCUGGCCAAGUAACCUCGAGAGCUACUUUGACCGUCUGUCUAUCAGGAUGAGAUCGCUGGUGCAUUGAAGGCCUACGAGAAAAUUUUACUGAGGCCACCCAGAACUUCAACACCAAAAGAUGACAGACUACGCCAAGAGGUGAGUGUCCUGGGCCCAACAACUACGGAUAGUUUUUGCCAGCCAGCAGAAGCCGGACACCAUUCCCACAGAACUGGCCAAACGGGUUCGAGUUAUGCCGGCAGCUGGAGAUGAAACCGAUCGUCUGAGCCCCGGGCACUGGUGGGCGGGCAGGGUCUACAAACAGUUCCGCAAGGUCCAAAGGUGGACGUCCAUCCUAAAGCCAAGC"
    # TODO RNA model files are 3' to 5' but in nanopolish they are 5' to 3'

    # fast5_path = "/Users/andrewbailey/CLionProjects/nanopore-RNN/test_files/minion-reads/rna_reads/one_rna/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_36_ch_218_strand.fast5"
    # nucleotide_seq = "UGCACAUUUAGCAUUGGCUGCUGUUGGUUUAGGCCGUCCGAACCUGGGCAGGAAGAUGUGGCCGCAAAGAAGACGAAAAGUCUGGAGUCGAUCAACUCUAGGCCUUAUUGCUAUGAAAGUGGGAAGUACGUCCUGGGUACAAGCAGACUCUGAAGAUGAUGAUCAGACACCAAGGGCAAGCGAAAUUGGUCAUUCGCUAACAACUGCCCAGCUUUGAGGAAAAUCUGAAAUAGAGUACUAUGUUAUGUUGGCUAAAACUGGUGCACAUCACUACGGUGGCAAUAAUAUUGAACUGCUGGGCACAGCAUCGGAAAGCUACUACAGAGUGCGCCAUUGGCUAUCAUUGAUCCAGGGGUGACUUGACCAUUAGAAGCUGCCAGAAAGACUGGUGAAAGUAAACCACACAAAAUUUUCAGCAAACUUCUAAACCUGCAUAAAAAUUCUUUAAUAAAUUCUGCUUGUUAAAAUUCCUCCAUCCUCCAUUCAUCCAUAUUAUCAUAUCAUAUCCCUUACCUAUCCUACAAAAUCCAA"
    #
    # nucleotide_seq_3_to_5 = nucleotide_seq[::-1].replace("U", "T")
    # nucleotide_seq_3_to_5 = nucleotide_seq.replace("U", "T")

    # create Fast5 object
    # RNA_MINKNOW = dict(window_lengths=(7, 14), thresholds=(2.5, 9.0), peak_height=1.0)
    # event_table, f5fh = create_minknow_events_from_fast5(fast5_path, **RNA_MINKNOW)
    # fast5_path = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/tests/minion_test_reads/canonical_ecoli_R9/miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_strand.fast5"
    #
    # model_file = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/models/testModelR9p4_5mer_acgt_RNA.model"

    # print(event_table)
    # fastq = f5fh.get_fastq(analysis="Basecall_1D", section="template")
    # sequence = fastq.split()[1]
    #
    # model_types = ["threeState", "threeStateHdp"]
    # model = get_model(model_types[0], model_file)
    # create events
    # events, sum_emission = simple_banded_event_align(event_table, model, nucleotide_seq_3_to_5)
    # events, sum_emission = simple_banded_event_align(event_table, model, nucleotide_seq_3_to_5)

    # embed
    # name = "SimpleBandedAlignment_00{}"
    # keys = ["log(total_probability)", "time_stamp"]
    # values = [sum_emission, TimeStamp().posix_date()]
    # attributes = dict(zip(keys, values))
    # # f5fh = Fast5(fast5_path, read='r+')
    # f5fh.set_new_event_table(name, events, attributes, overwrite=False)

    #DNA
    model_file = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/models/testModelR9_5mer_acgt_template.model"
    model_types = ["threeState", "threeStateHdp"]
    model = get_model(model_types[0], model_file)

    fast5_path = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/tests/minion_test_reads/1D/LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch92_read1108_strand.fast5"
    event_table, f5fh = create_minknow_events_from_fast5(fast5_path)

    fastq = f5fh.get_fastq(analysis="Basecall_1D", section="template")
    nucleotide_seq = fastq.split('\n')[1]

    name = "AdaptiveBandedAlignment_00{}"
    events, sum_emission = adaptive_banded_simple_event_align(event_table[:30], model, nucleotide_seq[:25], debug=False)
    # embed
    keys = ["log(total_probability)", "time_stamp"]
    values = [sum_emission, TimeStamp().posix_date()]
    attributes = dict(zip(keys, values))
    f5fh.set_new_event_table(name, events, attributes, overwrite=True)


if __name__ == '__main__':
    main()
