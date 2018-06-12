#!/usr/bin/env python
"""Re-segment raw read and re-label the data"""
########################################################################
# File: event_detection.py
#  executable: event_detection.py
#
# Author: Andrew Bailey
# History: 10/6/17 Created
########################################################################


from __future__ import print_function
import sys
import os
import collections
import re
import numpy as np
import traceback
from collections import defaultdict
from timeit import default_timer as timer
from signalalign.utils.pyporeParsers import SpeedyStatSplit
from signalalign.fast5 import Fast5
from signalalign.utils.filters import minknow_event_detect, scrappie_event_detect
from py3helpers.utils import check_numpy_table, list_dir, TimeStamp, change_np_field_type, merge_dicts
from py3helpers.seq_tools import create_fastq_line, check_fastq_line, ReverseComplement, pairwise_alignment_accuracy

EVENT_DETECT_SPEEDY = "speedy"
EVENT_DETECT_MINKNOW = "minknow"
EVENT_DETECT_SCRAPPIE = "scrappie"

def create_speedy_event_table(signal, sampling_freq, start_time, min_width=5, max_width=80, min_gain_per_sample=0.008,
                              window_width=800):
    """Create new event table using SpeedyStatSplit Event detection

    :param signal: list or array of signal in pA for finding events
    :param sampling_freq: sampling frequency of ADC in Hz
    :param start_time: start time from fast5 file (time in seconds * sampling frequency)
    :param min_width: param for SpeedyStatSplit
    :param max_width: param for SpeedyStatSplit
    :param min_gain_per_sample: param for SpeedyStatSplit
    :param window_width: param for SpeedyStatSplit
    :return: Table of events without model state or move information
    """
    assert np.sign(start_time) == 1, "Start time has to be positive: {}".format(start_time)
    assert type(signal[0]) is np.float64, "Signal needs to be in pA. Not ADC counts"

    # define speedy stat split
    parser = SpeedyStatSplit(min_width=min_width, max_width=max_width,
                             min_gain_per_sample=min_gain_per_sample,
                             window_width=window_width, sampling_freq=sampling_freq)
    # parse events
    events = parser.parse(np.asarray(signal, dtype=np.float64))
    num_events = len(events)
    # create empty event table
    event_table = get_empty_event_table(num_events)
    # set events into event table
    for i, event in enumerate(events):
        event_table['start'][i] = event.start / sampling_freq + (start_time / sampling_freq)
        event_table['raw_start'][i] = event.start
        event_table['length'][i] = event.duration / sampling_freq
        event_table['raw_length'][i] = event.duration
        event_table['mean'][i] = event.mean
        event_table['stdv'][i] = event.std

    return event_table


def create_minknow_event_table(signal, sampling_freq, start_time,
                               window_lengths=(3, 6), thresholds=(1.4, 9.0), peak_height=0.2):
    """Create new event table using minknow_event_detect event detection

    :param signal: list or array of signal in pA for finding events
    :param sampling_freq: sampling frequency of ADC in Hz
    :param start_time: start time from fast5 file (time in seconds * sampling frequency)
    :param window_lengths: t-test windows for minknow_event_detect
    :param thresholds: t-test thresholds for minknow_event_detect
    :param peak_height: peak height param for minknow_event_detect
    :return: Table of events without model state or move information
    """
    assert np.sign(start_time) == 1, "Start time has to be positive: {}".format(start_time)
    assert type(signal[0]) is np.float64, "Signal needs to be in pA. Not ADC counts"
    events = minknow_event_detect(np.asarray(signal, dtype=float), sample_rate=sampling_freq,
                                  get_peaks=False, window_lengths=window_lengths,
                                  thresholds=thresholds, peak_height=peak_height)
    num_events = len(events)
    event_table = get_empty_event_table(num_events)
    for i, event in enumerate(events):
        event_table['start'][i] = event["start"] + (start_time / sampling_freq)
        event_table['length'][i] = event["length"]
        event_table['mean'][i] = event["mean"]
        event_table['stdv'][i] = event["stdv"]
        event_table['raw_start'][i] = np.round(event["start"] * sampling_freq)
        event_table['raw_length'][i] = np.round(event["length"] * sampling_freq)

    return event_table


def create_scrappie_event_table(fast5_location, sampling_freq):
    """Create new event table using minknow_event_detect event detection

    :param signal: list or array of signal in pA for finding events
    :param sampling_freq: sampling frequency of ADC in Hz
    :return: Table of events without model state or move information
    """
    events = scrappie_event_detect(fast5_location)
    assert events is not None and len(events)>0, "Could not use scrappie to detect events in {}".format(fast5_location)
    num_events = len(events)
    event_table = get_empty_event_table(num_events - 1)

    for i, event in enumerate(events):
        # drop last event (cannot compute length without "next" event)
        if i == num_events - 1: break

        event_table['start'][i] = event["start"] / sampling_freq
        event_table['length'][i] = 0.0
        event_table['mean'][i] = event["mean"]
        event_table['stdv'][i] = event["stdv"]
        event_table['raw_start'][i] = np.round(event["start"])
        event_table['raw_length'][i] = 0

        if i != 0:
            event_table['length'][i-1] = event_table['start'][i] - event_table['start'][i - 1]
            event_table['raw_length'][i-1] = event_table['raw_start'][i] - event_table['raw_start'][i - 1]

    return event_table


def get_empty_event_table(num_events):
    """
    Gets default event table returned by create_XXX_event_table
    :param num_events: number of events to initialize
    :return: empty numpy array with appropriate columns
    """
    return np.empty(num_events, dtype=[('start', float), ('length', float),
                                              ('mean', float), ('stdv', float),
                                              ('model_state', 'S5'), ('move', '<i4'),
                                              ('raw_start', int), ('raw_length', int),
                                              ('p_model_state', float)])


def create_anchor_kmers(new_events, old_events):
    """
    Create anchor kmers for new event table.

    Basically, grab kmer and move information from previous event table and
    pull events covering the same time span into new event table.
    :param new_events: new event table
    :param old_events: event table from Fast5 file
    :return New event table
    """
    num_old_events = len(old_events)
    check_numpy_table(new_events, req_fields=('start', 'length', 'mean', 'stdv', 'model_state', 'move', 'p_model_state'))
    check_numpy_table(old_events, req_fields=('start', 'length', 'mean', 'stdv', 'model_state', 'move', 'p_model_state'))
    # index of old events
    old_indx = 0
    # start index to trim new_events for those with data from old_events
    start_index = 0
    end_index = len(new_events)
    # personal tracker for dealing with how the segmentation algorithm is working
    most_moves = 0
    # tracking overlaped events
    selected_overlap = False
    check_overlap = False
    homopolymer = False
    # keep track of events passed
    last_left_over = 0
    for i, event in enumerate(new_events):
        # skip events that occur before labels from old events
        if old_events[0]["start"] <= event["start"]:
            # time of old event in new event for a given kmer
            time = []
            probs = []
            moves = []
            kmers = []
            # new event's start and end
            current_event_start = round(event["start"], 7)
            current_event_end = round(current_event_start + event["length"], 7)
            # if first event or event start is after current old_event start.
            if old_indx != num_old_events:
                prev_kmer = str()
                num_loops = 0
                # print(round(old_events[old_indx]["start"], 7), old_events[old_indx]["length"], current_event_end, round(old_events[old_indx]["start"], 7) < current_event_end)
                while round(old_events[old_indx]["start"], 7) < current_event_end and old_indx != num_old_events:
                    # print("INSIDE LOOP", round(old_events[old_indx]["start"], 7), old_events[old_indx]["length"], current_event_end, round(old_events[old_indx]["start"], 7) < current_event_end)
                    # deal with bad event files and final event
                    if old_indx == num_old_events-1:
                        old_event_end = round(old_events[old_indx]["start"] + old_events[old_indx]["length"], 7)
                    else:
                        old_event_end = round(old_events[old_indx+1]["start"], 7)
                    old_event_start = round(old_events[old_indx]["start"], 7)
                    old_kmer = bytes.decode(old_events[old_indx]["model_state"])
                    # homopolymers or stays should be tracked together
                    if old_kmer == prev_kmer:
                        if len(set(old_kmer)) == 1:
                            if not homopolymer and selected_overlap and num_loops <= 1:
                                moves[index] = 0
                            homopolymer = True
                        else:
                            homopolymer = False
                        index = kmers.index(old_kmer)
                        probs[index] = max(probs[index], old_events[old_indx]["p_model_state"])
                        moves[index] += old_events[old_indx]["move"]
                    else:
                        # add new kmer
                        index = len(time)
                        kmers.append(old_kmer)
                        probs.append(old_events[old_indx]["p_model_state"])
                        moves.append(old_events[old_indx]["move"])
                        time.append(0)
                        homopolymer = False
                    prev_kmer = old_kmer
                    # if old event passes through current event calculate correct time in current event
                    # deal with old events ending after the new event end
                    if old_event_end > current_event_end:
                        time[index] += current_event_end - old_event_start
                        new_check_overlap = True
                        break
                    # check if entire old event is within the new event or not
                    else:
                        if old_event_start < current_event_start:
                            time[index] += old_event_end - current_event_start
                        else:
                            time[index] += old_event_end - old_event_start
                        # if old_event_end != current_event_end:
                        old_indx += 1
                        new_check_overlap = False
                    num_loops += 1
                    # break loop at end of old events
                    if old_indx == num_old_events:
                        break
            else:
                end_index = i
            num_kmers = len(kmers)
            # select index of best kmer to assign
            if num_kmers == 1:
                best_index = 0
                left_over = 0
            elif num_kmers > 1:
                # select on time in new event only
                best_index = time.index(max(time))
                # if there are several old events in a new event, track how many
                if new_check_overlap:
                    left_over = sum(moves[best_index+1:-1])
                else:
                    left_over = sum(moves[best_index+1:])
            else:
                # end of possible alignments
                end_index = i
                break
            # if previous old event overlapped into current new event
            # check if old event is going to be assigned twice
            if selected_overlap and best_index == 0 and check_overlap:
                if homopolymer:
                    move = moves[best_index]
                else:
                    move = 0
            elif selected_overlap and best_index != 0 and check_overlap:
                move = min(5, moves[best_index] + last_left_over)
            else:
                move = min(5, moves[best_index]+sum(moves[:best_index])+last_left_over)
                if most_moves < moves[best_index]+sum(moves[:best_index])+last_left_over:
                    most_moves = moves[best_index]+sum(moves[:best_index])+last_left_over
            # print(kmers, moves, left_over, moves[best_index], sum(moves[:best_index]), last_left_over, move)
            # if new overlap
            if new_check_overlap:
                # new overlapped event will be tracked on next new_event so we drop a left_over count
                left_over = max(0, left_over-1)
                if most_moves < left_over-1:
                    most_moves = left_over-1

                # check if we currently selected an overlapping old event
                if best_index == num_kmers-1:
                    selected_overlap = True
                else:
                    selected_overlap = False
            else:
                selected_overlap = False

            kmer = kmers[best_index]
            prob = probs[best_index]
            # assign event probs, move and model state
            event["p_model_state"] = prob
            event["move"] = move
            event["model_state"] = kmer
            check_overlap = new_check_overlap
            last_left_over = left_over
            new_check_overlap = False
            homopolymer = False
        else:
            # skip event since the
            start_index = i + 1
    # print(most_moves)
    return new_events[start_index:end_index]


def check_event_table_time(event_table):
    """Check if event table has correct math for start and length timing for each event

    :param event_table: event table with "start" and "length" columns
    """
    check_numpy_table(event_table, req_fields=('start', 'length'))

    prev_end = event_table[0]["start"] + event_table[0]["length"]
    for event in event_table[1:]:
        if prev_end != event["start"]:
            return False
        prev_end = event["start"]+event["length"]

    return True


def resegment_reads(fast5_path, params=None, speedy=False, overwrite=True, analysis_path="ReSegmentBasecall_000"):
    """Re-segment and create anchor alignment from previously base-called fast5 file
    :param fast5_path: path to fast5 file
    :param params: event detection parameters
    :param speedy: boolean option for speedyStatSplit or minknow
    :param overwrite: overwrite a previous event re-segmented event table
    :param analysis_path: name of key where events table will be placed (Analyses/'name'/Events)
    :return True when completed
    """
    assert os.path.isfile(fast5_path), "File does not exist: {}".format(fast5_path)
    # create Fast5 object and sanity check
    f5fh = Fast5(fast5_path, read='r+')
    if not f5fh.has_basecall_data():
        f5fh.close()
        return None

    # gather previous event detection
    old_event_table = f5fh.get_basecall_data()

    read_id = bytes.decode(f5fh.raw_attributes['read_id'])
    sampling_freq = f5fh.sample_rate
    start_time = f5fh.raw_attributes['start_time']

    # get params
    if params is None: params = get_default_event_detection_params(
        EVENT_DETECT_SPEEDY if speedy else EVENT_DETECT_MINKNOW)

    # pick event detection algorithm
    signal = f5fh.get_read(raw=True, scale=True)
    if speedy:
        event_table = create_speedy_event_table(signal, sampling_freq, start_time, **params)
        params = merge_dicts([params, {"event_detection": "speedy_stat_split"}])
    else:
        event_table = create_minknow_event_table(signal, sampling_freq, start_time, **params)
        params = merge_dicts([params, {"event_detection": "minknow_event_detect"}])

    # metadata
    keys = ["nanotensor version", "time_stamp"]
    values = ["0.2.0", TimeStamp().posix_date()]
    attributes = merge_dicts([params, dict(zip(keys, values)), f5fh.raw_attributes])

    # do resegmentation
    if f5fh.is_read_rna():
        old_event_table = index_to_time(old_event_table, sampling_freq=sampling_freq, start_time=start_time)
    new_event_table = create_anchor_kmers(new_events=event_table, old_events=old_event_table)

    # get destination in fast5
    #todo find latest location? ie: save_event_table_and_fastq(..)
    destination = analysis_path

    f5fh.set_event_table(destination, new_event_table, attributes, overwrite=overwrite)

    # gather new sequence
    sequence = sequence_from_events(new_event_table)
    if f5fh.is_read_rna():
        sequence = ReverseComplement().reverse(sequence)
        sequence = sequence.replace("T", "U")
    quality_scores = '!'*len(sequence)
    fastq = create_fastq_line(read_id+" :", sequence, quality_scores)

    # set fastq
    f5fh.set_fastq(destination, fastq)
    return f5fh


def get_default_event_detection_params(event_detection_strategy):
    if event_detection_strategy == EVENT_DETECT_SPEEDY:
        return dict(min_width=5, max_width=80, min_gain_per_sample=0.008, window_width=800)
    elif event_detection_strategy == EVENT_DETECT_MINKNOW:
        return dict(window_lengths=(5, 10), thresholds=(2.0, 1.1), peak_height=1.2)
    elif event_detection_strategy == EVENT_DETECT_SCRAPPIE:
        return {}
    else:
        return None


def generate_events_and_alignment(fast5_path, nucleotide_sequence, nucleotide_qualities=None,
                                  event_detection_params=None, event_detection_strategy=None,
                                  save_to_fast5=True, overwrite=False,
                                  analysis_identifier=Fast5.__default_basecall_1d_analysis__, ):

    assert os.path.isfile(fast5_path), "File does not exist: {}".format(fast5_path)

    # create Fast5 object
    f5fh = Fast5(fast5_path, read='r+')
    read_id = bytes.decode(f5fh.raw_attributes['read_id'])
    sampling_freq = f5fh.sample_rate
    start_time = f5fh.raw_attributes['start_time']
    success = False

    # event detection prep
    if event_detection_strategy is None:
        event_detection_strategy = EVENT_DETECT_MINKNOW
    if event_detection_params is None:
        event_detection_params = get_default_event_detection_params(event_detection_strategy)

    # detect events
    if event_detection_strategy == EVENT_DETECT_SPEEDY:
        signal = f5fh.get_read(raw=True, scale=True)
        event_table = create_speedy_event_table(signal, sampling_freq, start_time, **event_detection_params)
        event_detection_params = merge_dicts([event_detection_params, {"event_detection": "speedy_stat_split"}])
    elif event_detection_strategy == EVENT_DETECT_MINKNOW:
        signal = f5fh.get_read(raw=True, scale=True)
        event_table = create_minknow_event_table(signal, sampling_freq, start_time, **event_detection_params)
        event_detection_params = merge_dicts([event_detection_params, {"event_detection": "minknow_event_detect"}])
    elif event_detection_strategy == EVENT_DETECT_SCRAPPIE:
        event_table = create_scrappie_event_table(fast5_path, sampling_freq)
        event_detection_params = merge_dicts([event_detection_params, {"event_detection": "scrappie_event_detect"}])
    else:
        raise Exception("PROGRAMMER ERROR: unknown resegment strat {}: expected {}"
                        .format(event_detection_strategy, [EVENT_DETECT_SPEEDY, EVENT_DETECT_MINKNOW, EVENT_DETECT_SCRAPPIE]))

    # gather attributes
    keys = ["nanotensor version", "time_stamp"]
    values = ["0.2.0", TimeStamp().posix_date()]
    attributes = merge_dicts([event_detection_params, dict(zip(keys, values)), f5fh.raw_attributes])

    # do the alignment
    # todo do_alignment(events, nucleotide_sequence)
    # success = evaluate_success()

    # save to fast5 (if appropriate)
    saved_location = None
    if save_to_fast5:
        fastq = create_fastq_line(read_id, nucleotide_sequence,
                                  "*" if nucleotide_qualities is None else nucleotide_qualities)
        saved_location = save_event_table_and_fastq(f5fh, event_table, fastq, attributes=attributes,
                                                    overwrite=overwrite, analysis_identifier=analysis_identifier)

    # close
    f5fh.close()

    return success, event_table, saved_location


def save_event_table_and_fastq(f5fh, event_table, fastq_line, attributes=None, overwrite=False,
                                    analysis_identifier=Fast5.__default_basecall_1d_analysis__):
    # get destination root
    if overwrite:
        try:
            # get current identifier
            destination = f5fh.get_analysis_latest(analysis_identifier)
            f5fh.delete(destination, ignore=True)
        except IndexError:
            # doesn't exist, get initial location (ie Basecall_000)
            destination = f5fh.get_analysis_new(analysis_identifier)
    else:
        # get incremented directory for identifier
        destination = f5fh.get_analysis_new(analysis_identifier)

    # set event table (if present)
    if event_table is not None:
        f5fh.set_event_table(destination, event_table, attributes)

    # set nucleotide sequence (if present)
    if fastq_line is not None:
        f5fh.set_fastq(destination, fastq_line)

    return destination


def index_to_time(basecall_events, sampling_freq=0, start_time=0):
    """Convert RNA basecall read start and length from indexes to time stamps

    :param basecall_events: basecall events from albacore/metricore basecalled event table
    :param sampling_freq: sampling frequency of experiment
    :param start_time: start time of experiment via fasta5 file
    """
    check_numpy_table(basecall_events, req_fields=('start', 'length'))
    assert basecall_events["start"].dtype is np.dtype('uint64'), "Event start should be np.int32 type: {}"\
        .format(basecall_events["start"].dtype)
    assert sampling_freq != 0, "Must set sampling frequency"
    assert start_time != 0, "Must set start time"

    event_table = change_np_field_type(basecall_events, 'start', float)
    event_table = change_np_field_type(event_table, 'length', float)
    event_table["start"] = (event_table["start"] / sampling_freq) + (start_time / sampling_freq)
    event_table["length"] = event_table["length"] / float(sampling_freq)
    return event_table


def time_to_index(event_table, sampling_freq=0, start_time=0):
    """Convert start and lengths from time to raw signal indexes

    :param event_table: basecall events from albacore/metricore basecalled event table
    :param sampling_freq: sampling frequency of experiment
    :param start_time: start time of experiment via fasta5 file
    """
    check_numpy_table(event_table, req_fields=('start', 'length'))
    assert event_table["start"].dtype is not np.dtype('uint64'), "Event start should not be np.int32 type: {}" \
        .format(event_table["start"].dtype)
    assert sampling_freq != 0, "Must set sampling frequency"
    assert start_time != 0, "Must set start time"

    event_table["start"] = np.round((event_table["start"] - (start_time / float(sampling_freq))) * sampling_freq)
    event_table["length"] = np.round(event_table["length"] * sampling_freq)
    event_table = change_np_field_type(event_table, 'start', int)
    event_table = change_np_field_type(event_table, 'length', int)

    return event_table


def sequence_from_events(events):
    """Get new read from event table with 'model_state' and 'move' fields

    :param events: event table with 'model_state' and 'move' fields

    """
    check_numpy_table(events, req_fields=("model_state", "move"))
    bases = []
    for i, event in enumerate(events):
        if i == 0:
            bases.extend([chr(x) for x in event['model_state']])

        else:
            if event['move'] > 0:
                bases.append(bytes.decode
                             (event['model_state'][-event['move']:]))
    sequence = ''.join(bases)
    return sequence


def get_resegment_accuracy(fast5handle, section="template"):
    """Get accuracy comparison between original sequence and resegmented generated sequence

    :param fast5handle: Fast5 object with re-segemented read
    """
    assert isinstance(fast5handle, Fast5), "fast5handle needs to be a Fast5 instance"
    # get fastqs
    resegment_fastq = fast5handle.get_fastq(analysis="ReSegmentBasecall", section=section)
    original_fastq = fast5handle.get_fastq(analysis="Basecall_1D", section=section)[:-1]
    # make sure the underlying assumption that we can split on newline is ok
    check_fastq_line(resegment_fastq)
    check_fastq_line(original_fastq)
    # get sequence
    resegment_seq = resegment_fastq.split('\n')[1]
    original_seq = original_fastq.split('\n')[1]
    return pairwise_alignment_accuracy(original_seq, resegment_seq, soft_clip=True)


def create_minknow_events_from_fast5(fast5_path, window_lengths=(3, 6), thresholds=(1.4, 9.0), peak_height=0.2):
    """Create events with ('start', 'length', 'mean', 'stdv', 'model_state', 'move', 'p_model_state') fields from
        fast5 file. The 'model_state', 'move' and 'p_model_state' are all empty

    :param fast5_path: path to fast5 file
    :param window_lengths: Length 2 list of window lengths across
        raw data from which `t_stats` are derived
    :param thresholds: Length 2 list of thresholds on t-statistics
    :param peak_height: Absolute height a peak in signal must rise below
        previous and following minima to be considered relevant
    """
    assert os.path.isfile(fast5_path), "File does not exist: {}".format(fast5_path)
    f5fh = Fast5(fast5_path, read='r+')
    signal = f5fh.get_read(raw=True, scale=True)
    # read_id = bytes.decode(f5fh.raw_attributes['read_id'])
    sampling_freq = f5fh.sample_rate
    start_time = f5fh.raw_attributes['start_time']
    #
    event_table = create_minknow_event_table(signal, sampling_freq, start_time, window_lengths=window_lengths,
                                             thresholds=thresholds, peak_height=peak_height)

    return event_table, f5fh


def main():
    """Main docstring"""
    start = timer()

    dna_reads = "/Users/andrewbailey/CLionProjects/nanopore-RNN/test_files/minion-reads/canonical/"
    rna_reads = "/Users/andrewbailey/CLionProjects/nanopore-RNN/test_files/minion-reads/rna_reads"

    dna_minknow_params = dict(window_lengths=(5, 10), thresholds=(2.0, 1.1), peak_height=1.2)
    dna_speedy_params = dict(min_width=5, max_width=80, min_gain_per_sample=0.008, window_width=800)
    rna_minknow_params = dict(window_lengths=(5, 10), thresholds=(2.0, 1.1), peak_height=1.2)
    rna_speedy_params = dict(min_width=5, max_width=40, min_gain_per_sample=0.008, window_width=800)


    rna_minknow_params = dict(window_lengths=(5, 10), thresholds=(1.9, 1.0), peak_height=1.2)
    rna_speedy_params = dict(min_width=5, max_width=40, min_gain_per_sample=0.008, window_width=800)
    dna_minknow_params = dict(window_lengths=(5, 10), thresholds=(2.0, 1.1), peak_height=1.2)
    dna_speedy_params = dict(min_width=5, max_width=80, min_gain_per_sample=0.008, window_width=800)

    rna_files = list_dir(rna_reads, ext='fast5')
    dna_files = list_dir(dna_reads, ext='fast5')
    print("MAX RNA SKIPS: Speedy")
    for fast5_path in rna_files:
        print(fast5_path)
        f5fh = resegment_reads(fast5_path, rna_speedy_params, speedy=True, overwrite=True)
        print(get_resegment_accuracy(f5fh))
        # f5fh = resegment_reads(fast5_path, rna_minknow_params, speedy=False, overwrite=True)
        # print(test_resegment_accuracy(f5fh))

    print("MAX RNA SKIPS: Minknow")
    for fast5_path in rna_files:
        f5fh = resegment_reads(fast5_path, rna_minknow_params, speedy=False, overwrite=True)
        print(get_resegment_accuracy(f5fh))

    print("MAX DNA SKIPS: speedy")
    for fast5_path in dna_files:
            print(fast5_path)
            f5fh = resegment_reads(fast5_path, dna_speedy_params, speedy=True, overwrite=True)
            print(get_resegment_accuracy(f5fh))
    print("MAX DNA SKIPS:Minknow")
    for fast5_path in dna_files:
        f5fh = resegment_reads(fast5_path, dna_minknow_params, speedy=False, overwrite=True)
        print(get_resegment_accuracy(f5fh))

        # print(fast5_path)
    stop = timer()
    print("Running Time = {} seconds".format(stop - start), file=sys.stderr)


if __name__ == "__main__":
    main()

    raise SystemExit
