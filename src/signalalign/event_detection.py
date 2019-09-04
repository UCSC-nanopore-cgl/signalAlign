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
import subprocess
import numpy as np
import h5py
import traceback
import tempfile
from numpy.lib.recfunctions import append_fields
from shutil import which
from contextlib import closing
from collections import defaultdict
from timeit import default_timer as timer
# from signalalign.utils.pyporeParsers import SpeedyStatSplit
from signalalign.fast5 import Fast5
from signalalign.utils.sequenceTools import get_full_nucleotide_read_from_alignment, reverse_complement
from signalalign.utils.filters import minknow_event_detect
from py3helpers.utils import check_numpy_table, list_dir, TimeStamp, change_np_field_type, merge_dicts
from py3helpers.seq_tools import create_fastq_line, check_fastq_line, ReverseComplement, pairwise_alignment_accuracy


EVENT_KMERALIGN_TMP = "KmerEventAlign_tmp"


def check_event_table_time(event_table, min_difference = 0.0000001):
    """Check if event table has correct math for start and length timing for each event

    :param event_table: event table with "start" and "length" columns
    """
    check_numpy_table(event_table, req_fields=('start', 'length'))

    prev_end = event_table[0]["start"] + event_table[0]["length"]
    for event in event_table[1:]:
        if prev_end - event["start"] > min_difference:
            return False
        prev_end = event["start"] + event["length"]

    return True



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


def add_raw_start_and_raw_length_to_events(basecall_events, sampling_freq, start_time):
    check_numpy_table(basecall_events, req_fields=('start', 'length'))
    assert basecall_events["start"].dtype == np.dtype('float64'), "Event start should be float64 type: {}" \
        .format(basecall_events["start"].dtype)
    assert basecall_events["length"].dtype == np.dtype('float64'), "Event length should be float64 type: {}" \
        .format(basecall_events["length"].dtype)
    assert sampling_freq >= 0, "Invalid sampling frequency: {}".format(sampling_freq)
    assert start_time != 0, "Invalid start time: {}".format(start_time)

    calc_raw_start = lambda x: np.uint64(np.round((basecall_events[x]["start"] -
                                                   (start_time / float(sampling_freq))) * sampling_freq))
    calc_raw_length = lambda x: np.uint64(np.round(basecall_events[x]["length"] * sampling_freq))
    raw_starts = list(map(calc_raw_start, range(len(basecall_events))))
    raw_lengths = list(map(calc_raw_length, range(len(basecall_events))))
    basecall_events = append_fields(basecall_events, "raw_start", raw_starts, usemask=False)
    basecall_events = append_fields(basecall_events, "raw_length", raw_lengths, usemask=False)
    return basecall_events


def add_start_and_length_to_events(basecall_events, sampling_freq, start_time):
    check_numpy_table(basecall_events, req_fields=('raw_start', 'raw_length'))
    assert basecall_events["raw_start"].dtype == np.dtype('uint64'), "Event raw_start should be uint64 type: {}" \
        .format(basecall_events["raw_start"].dtype)
    assert basecall_events["raw_length"].dtype == np.dtype('uint64'), "Event raw_length should be uint64 type: {}" \
        .format(basecall_events["raw_length"].dtype)
    assert sampling_freq >= 0, "Invalid sampling frequency: {}".format(sampling_freq)
    assert start_time != 0, "Invalid start time: {}".format(start_time)

    calc_start = lambda x: np.float64((basecall_events[x]["raw_start"] / float(sampling_freq)) +
                                      (start_time / float(sampling_freq)))
    calc_length = lambda x: np.float64(basecall_events[x]["raw_length"] / float(sampling_freq))
    starts = list(map(calc_start, range(len(basecall_events))))
    lengths = list(map(calc_length, range(len(basecall_events))))
    basecall_events = append_fields(basecall_events, "start", starts, usemask=False)
    basecall_events = append_fields(basecall_events, "length", lengths, usemask=False)
    return basecall_events


def load_from_raw(np_handle, alignment_file, model_file_location, path_to_bin="./",
                  nucleotide_sequence=None, analysis_identifier=None, write_failed_alignments=False, rna=False):
    """Load a nanopore read from raw signal and an alignment file. Need a model to create banded alignment.
    :param np_handle: NanoporeRead class object
    :param alignment_file: sam/bam file
    :param model_file_location: path to model file
    :param path_to_bin: bath to signalAlign bin where executables are stored
    :param nucleotide_sequence: nucleotide sequence (needed if no alignment file is available)
    :param analysis_identifier: identifier for storage of event table and fastq
    :param write_failed_alignments: still write alignments that failed quality checks
    :param rna: boolean option to assume rna read
    :return: path to events in fast5 file or -1 if the task fails
    """
    assert os.path.isfile(model_file_location), \
        "Model_file_location must be a real path to a SignalAlign HMM model file"
    assert os.path.exists(path_to_bin), \
        "path_to_bin must exist"
    if not os.path.isfile(str(alignment_file)) and nucleotide_sequence is None:
        nucleotide_sequence = np_handle.get_template_read(initalize_bypass=True)
        assert nucleotide_sequence, "alignment_file must be a real path a SAM/BAM alignment file, or " \
                                    "nucleotide_sequence must be specified (retrieval attempted from fast5). " \
                                    "alignment_file: {}, nucleotide_sequence:{}".format(alignment_file, nucleotide_sequence)

    # check if file is open
    if not np_handle.open():
        return False
    # grab read id
    read_id = np_handle.read_label

    # get nucleotides and qualities
    if nucleotide_sequence is None:
        # get/build nucleotide sequence from alignment file (accounting for hardclipping)
        nucleotide_sequence, nucleotide_qualities, _, _, _ = \
            get_full_nucleotide_read_from_alignment(alignment_file, read_id)
        if nucleotide_sequence is None:
            print("[load_from_raw] nucleotides for {} not found in {}".format(read_id, alignment_file), file=sys.stderr)
            return False
    else:
        nucleotide_qualities = None
    if nucleotide_qualities is None:
        nucleotide_qualities = "!" * len(nucleotide_sequence)

    # get fastq (this is saved with the event table)
    fastq = create_fastq_line(read_id, nucleotide_sequence, nucleotide_qualities)

    # get temp location
    tmp_root = np_handle.fastFive.get_analysis_new(EVENT_KMERALIGN_TMP)
    tmp_dest = np_handle.fastFive.get_analysis_events_path_new(EVENT_KMERALIGN_TMP)
    assert tmp_dest.startswith(tmp_root), "Invalid analysis path management"
    file_name = np_handle.filename
    np_handle.close()
    tmp_directory = tempfile.mkdtemp()
    # run the c code which does the required stuff
    status = run_kmeralign_exe(file_name, nucleotide_sequence, model_file_location, tmp_dest, path_to_bin,
                               write_failed_alignments=write_failed_alignments, tmp_directory=tmp_directory,
                               rna=rna)
    os.removedirs(tmp_directory)
    # alignment succeeded, save it to the appropriate location
    if status:
        np_handle.open()
        if analysis_identifier is None: analysis_identifier = Fast5.__default_basecall_1d_analysis__
        # get attrs
        keys = ["signalAlign version", "time_stamp"]
        values = ["0.2.0", TimeStamp().posix_date()]
        attributes = merge_dicts([dict(zip(keys, values)), np_handle.fastFive.raw_attributes])
        # get events (and delete tmp location)
        events = np_handle.fastFive.get_custom_analysis_events(EVENT_KMERALIGN_TMP)
        np_handle.fastFive.delete(tmp_root, ignore=False)
        # save events and fastq
        saved_loc = save_event_table_and_fastq(np_handle.fastFive, events, fastq, attributes,
                                               analysis_identifier=analysis_identifier)
        return saved_loc

    # alignment failed, remove offending location (if it exists) and report
    else:
        print("[load_from_raw] error performing kmeralign", file=sys.stderr)
        np_handle.open()
        np_handle.fastFive.delete(tmp_root, ignore=True)
        return False


def load_from_raw2(np_handle, aligned_segment, model_file_location, path_to_bin="./",
                   analysis_identifier=None, write_failed_alignments=False, rna=False):
    """Load a nanopore read from raw signal and an alignment file. Need a model to create banded alignment.
    :param np_handle: NanoporeRead class object
    :param aligned_segment: pysam aligned_segment object
    :param model_file_location: path to model file
    :param path_to_bin: bath to signalAlign bin where executables are stored
    :param analysis_identifier: identifier for storage of event table and fastq
    :param write_failed_alignments: still write alignments that failed quality checks
    :param rna: boolean option for rna reads
    :return: path to events in fast5 file or -1 if the task fails
    """
    assert os.path.isfile(model_file_location), \
        "Model_file_location must be a real path to a SignalAlign HMM model file"
    assert os.path.exists(path_to_bin), \
        "path_to_bin must exist"
    # check if file is open
    if not np_handle.open():
        return False
    # grab read id
    read_id = np_handle.read_label

    # get nucleotides and qualities
    nucleotide_sequence = aligned_segment.query_sequence.upper()
    nucleotide_qualities = aligned_segment.qual

    # check for reverse mapping
    if aligned_segment.is_reverse:
        nucleotide_sequence = reverse_complement(nucleotide_sequence, reverse=True, complement=True)
        if nucleotide_qualities is not None and len(nucleotide_qualities) != 0:
            nucleotide_qualities = ''.join(reversed(list(nucleotide_qualities)))

    if nucleotide_qualities is None:
        nucleotide_qualities = "!" * len(nucleotide_sequence)

    # get fastq (this is saved with the event table)
    fastq = create_fastq_line(read_id, nucleotide_sequence, nucleotide_qualities)

    # get temp location
    tmp_root = np_handle.fastFive.get_analysis_new(EVENT_KMERALIGN_TMP)
    tmp_dest = np_handle.fastFive.get_analysis_events_path_new(EVENT_KMERALIGN_TMP)
    assert tmp_dest.startswith(tmp_root), "Invalid analysis path management"
    file_name = np_handle.filename
    np_handle.close()
    tmp_directory = tempfile.mkdtemp()
    # run the c code which does the required stuff
    status = run_kmeralign_exe(file_name, nucleotide_sequence, model_file_location, tmp_dest, path_to_bin,
                               write_failed_alignments=write_failed_alignments, tmp_directory=tmp_directory,
                               rna=rna)
    os.removedirs(tmp_directory)
    # alignment succeeded, save it to the appropriate location
    if status:
        np_handle.open()
        if analysis_identifier is None: analysis_identifier = Fast5.__default_basecall_1d_analysis__
        # get attrs
        keys = ["signalAlign version", "time_stamp"]
        values = ["0.2.0", TimeStamp().posix_date()]
        attributes = merge_dicts([dict(zip(keys, values)), np_handle.fastFive.raw_attributes])
        # get events (and delete tmp location)
        events = np_handle.fastFive.get_custom_analysis_events(EVENT_KMERALIGN_TMP)
        np_handle.fastFive.delete(tmp_root, ignore=False)
        # save events and fastq
        saved_loc = save_event_table_and_fastq(np_handle.fastFive, events, fastq, attributes,
                                               analysis_identifier=analysis_identifier)
        return saved_loc

    # alignment failed, remove offending location (if it exists) and report
    else:
        print("[load_from_raw] error performing kmeralign", file=sys.stderr)
        np_handle.open()
        np_handle.fastFive.delete(tmp_root, ignore=True)
        return False


def run_kmeralign_exe(fast5_path, nuc_sequence, model_file, dest, path_to_bin="./", tmp_directory=None,
                      delete_tmp_fasta=True, write_failed_alignments=False, rna=False):
    """Run kmerEventAlign. Generates alignment file
    :param fast5_path: path to fast5
    :param nuc_sequence: nucleotide sequence to align
    :param model_file: signal align model file
    :param dest: location to place events.
                    ex.  passing "Analyses/Basecalled_1D_template" will create "Analyses/Basecalled_1D_template/Events"
    :param path_to_bin: path to SignalAlign bin with executables
    :param write_failed_alignments: write failed alignments
    :param tmp_directory: if set, a fasta will be written here and used in kmerEventAlign invocation
    :param delete_tmp_fasta: if set and fasta is written, fasta will be deleted afterwards
    :param rna: if set assume is rna read
    :return:
    """
    executable = os.path.join(path_to_bin, "kmerEventAlign")
    fasta_location = None
    try:
        cmd = [executable, '-f', fast5_path, '-m', model_file, '-p', dest]
        if write_failed_alignments:
            cmd.append('-w')
        if tmp_directory is None:
            cmd.extend(['-N', nuc_sequence])
        else:
            fasta_location = os.path.join(tmp_directory, "{}.fa".format(
                ''.join(os.path.basename(fast5_path).split('.')[:-1])))
            with open(fasta_location, 'w') as fa_out:
                fa_out.write(">{}\n".format(os.path.basename(fast5_path)))
                fa_out.write("{}\n".format(nuc_sequence))
            cmd.extend(['-n', fasta_location])
        if rna:
            cmd.append('--rna')
        subprocess.check_call(cmd)
        status = True
    except Exception as e:
        print("Exception in run_kmeralign: {}".format(e))
        status = False
    finally:
        if fasta_location is not None and os.path.isfile(fasta_location) and delete_tmp_fasta:
            os.remove(fasta_location)

    return status
