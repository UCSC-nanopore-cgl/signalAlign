#!/usr/bin/env python
"""Class and methods to deal with aligned signal to reference"""
########################################################################
# File: alignedsignal.py
#  executable: alignedsignal.py
#
# Author: Andrew Bailey
# History: Created 03/09/18
########################################################################

import sys
import os
import subprocess
import numpy as np
from timeit import default_timer as timer
from collections import defaultdict, namedtuple
from py3helpers.utils import check_numpy_table
from py3helpers.seq_tools import ReferenceHandler, ReverseComplement, get_minimap_alignment, \
    initialize_aligned_segment_wrapper
from signalalign.fast5 import Fast5
from signalalign.mea_algorithm import maximum_expected_accuracy_alignment, mea_slow, \
    mea_slower, create_random_prob_matrix, get_mea_params_from_events, match_events_with_signalalign
from signalalign.event_detection import time_to_index
from itertools import islice


class AlignedSignal(object):
    """Labeled nanopore signal data"""

    def __init__(self, scaled_signal):
        """Initialize the scaled signal and label

        :param scaled_signal: scaled signal to pA
        """
        self.scaled_signal = None
        self.raw_signal = None
        self._add_scaled_signal(scaled_signal)
        self.signal_length = len(self.scaled_signal)
        self.minus_strand = None
        # label can be used for neural network training with all signal continuously labelled
        self.label = defaultdict()
        # predictions can have multiple labels for different sections of current
        self.prediction = defaultdict()
        # guides are sections that we are confident in (guide alignments)
        self.guide = defaultdict(defaultdict)

    def add_raw_signal(self, signal):
        """Add raw signal to class

        :param signal: raw current signal in ADC counts
        """
        assert int(signal[0]) == signal[0], "Raw signal are always integers"
        assert len(signal) == len(self.scaled_signal) and len(signal) == self.signal_length, \
            "Raw signal must be same size as scaled signal input:{} != scale:{}".format(signal, self.scaled_signal)
        self.raw_signal = signal

    def _add_scaled_signal(self, signal):
        """Add scaled signal to class

        :param signal: normalized current signal to pA
        """
        if type(signal) is np.ndarray:
            signal = signal.tolist()
        assert type(signal[0]) == float, "scaled signal must be a float"
        self.scaled_signal = signal

    def add_label(self, label, name, label_type, guide_name=None):
        """Add labels to class.

        :param label: label numpy array with required fields ['raw_start', 'raw_length', 'reference_index',
                                                              'kmer', 'posterior_probability']
        :param name: name of the label for signal
        :param label_type: type of label  :['label', 'prediction', 'guide']
        :param guide_name: must pass your own label via guide_name if label_type is guide
        """
        assert label_type in ['label', 'prediction', 'guide'], \
            "{} not in ['label', 'prediction', 'guide']: Must select an acceptable type".format(label_type)
        check_numpy_table(label, req_fields=('raw_start', 'raw_length', 'reference_index',
                                            'kmer', 'posterior_probability'))

        # label.sort(order=['raw_start'], kind='mergesort')
        # check the labels are in the correct format
        assert min(label["raw_start"]) >= 0, "Raw start cannot be less than 0"
        assert 0 <= max(label["posterior_probability"]) <= 1, \
            "posterior_probability must be between zero and one {}".format(row["posterior_probability"])
        if label_type == 'guide':
            assert guide_name is not None, "If label_type is 'guide', you must pass in a guide_name"
        # make sure last label can actually index the signal correctly
        try:
            self.scaled_signal[label[-1]["raw_start"]:label[-1]["raw_start"] + label[-1]["raw_length"]]
        except IndexError:
            raise IndexError("labels are longer than signal")

        label1 = np.sort(label, order=['raw_start'], kind='mergesort')

        # infer strand alignment of read
        if label1[0]["reference_index"] >= label1[-1]["reference_index"]:
            minus_strand = True
        else:
            minus_strand = False
        if self.minus_strand is not None:
            if label[0]["raw_start"] != label[-1]["raw_start"]:
                assert self.minus_strand == minus_strand, "New label has different strand direction, check label"
        else:
            self.minus_strand = minus_strand

        # set label with the specified name
        if label_type == 'label':
            self.label[name] = label
        elif label_type == 'prediction':
            self.prediction[name] = label
        elif label_type == 'guide':
            self.guide[guide_name][name] = label

    def generate_label_mapping(self, name, scaled=True):
        """Create a generator of the mapping between the signal and the label

        :param name: name of mapping to create label mapping
        :param scaled: boolean option for returning scaled or unscaled signal
        """
        assert name in self.label.keys(), "{} is not in labels dataset: {}".format(name, self.label.keys())
        label = self.label[name]
        len_label = len(label)
        if scaled:
            signal = self.scaled_signal
        else:
            assert self.raw_signal is not None, "Must set raw signal in order to generate raw signal alignments"
            signal = self.raw_signal
        for i, segment in enumerate(label):
            start = segment["raw_start"]
            if i < len_label - 1:
                end = label[i + 1]["raw_start"]
            else:
                end = segment["raw_start"] + segment["raw_length"]
            yield signal[start:end], segment['kmer'], segment['posterior_probability'], segment['reference_index']


class CreateLabels(Fast5):
    """Create an Aligned Signal object from a fast5 file with """

    def __init__(self, fast5_path):
        """Initialize fast5 object and keep track of AlignedSignal object
        :param fast5_path: path to fast5 file
        :param reference: path to reference so we can run signalAlign
        :param alignment_file: alignment file to run signalAlign or for alignment info
        """
        self.fast5_path = fast5_path
        super(CreateLabels, self).__init__(fast5_path)
        self.aligned_signal = self._initialize()
        self.kmer_index = 2
        self.rna = self.is_read_rna()

    def _initialize(self):
        """Initialize AlignedSignal class by adding the raw and scaled signal"""
        scaled_signal = self.get_read(raw=True, scale=True)
        raw_signal = self.get_read(raw=True, scale=False)
        # add raw signal information to AlignedSignal
        aligned_signal = AlignedSignal(scaled_signal)
        aligned_signal.add_raw_signal(raw_signal)
        return aligned_signal

    def add_mea_labels(self):
        """Gather mea_alignment labels information from fast5 file."""
        mea_alignment = self.get_signalalign_events(mea=True)
        # rna reference positions are on 5' edge aka right side of kmer
        if self.rna:
            mea_alignment["reference_index"] -= self.kmer_index
        else:
            mea_alignment["reference_index"] += self.kmer_index

        self.aligned_signal.add_label(mea_alignment, name="mea_signalalign", label_type='label')
        return True

    def add_signal_align_predictions(self):
        """Create prediction using probabilities from full output format from signalAlign"""
        sa_events = self.get_signalalign_events()
        # cut out duplicates
        sa_events = np.unique(sa_events)
        events = self.get_basecall_data()
        predictions = match_events_with_signalalign(sa_events=sa_events, event_detections=events)
        # rna reference positions are on 5' edge aka right side of kmer
        if self.rna:
            predictions["reference_index"] -= self.kmer_index
        else:
            predictions["reference_index"] += self.kmer_index
        self.aligned_signal.add_label(predictions, name="full_signalalign", label_type='prediction')
        return True

    def add_basecall_alignment(self, sam=None):
        """Add the original basecalled event table and add the 'guide' alignment labels to signal_label handle
        :param sam: correctly formatted SAM string
        :return: True if the correct labels are added to the AlignedSignal internal class object
        """
        if not sam:
            sam = self.get_signalalign_events(sam=True)

        events = self.get_basecall_data()
        try:
            check_numpy_table(events, req_fields=('raw_start', 'raw_length'))
        # if events do not have raw_start or raw_lengths
        except KeyError:
            events = time_to_index(events,
                                   sampling_freq=self.sample_rate,
                                   start_time=self.raw_attributes["start_time"])

        cigar_labels = create_labels_from_guide_alignment(events=events, sam_string=sam,
                                                          kmer_index=self.kmer_index)
        for i, block in enumerate(cigar_labels):
            self.aligned_signal.add_label(block, name="basecalled_alignment{}".format(i), label_type='guide',
                                          guide_name="basecall")
        return True

    def add_basecall_event_table_alignment(self, event_table_path=None):
        """Add guide alignment labels to signal_label handle"""
        assert event_table_path in self, "Event table path {}: is not in fast5.".format(event_table_path)
        test_sam = self.get_signalalign_events(sam=True)
        events = np.array(self[event_table_path])
        cigar_labels = create_labels_from_guide_alignment(events=events, sam_string=test_sam,
                                                          kmer_index=self.kmer_index)
        for i, block in enumerate(cigar_labels):
            # print(block)
            self.aligned_signal.add_label(block, name="{}{}".format(event_table_path, i), label_type='guide',
                                          guide_name=event_table_path)
        return True

    def add_nanoraw_labels(self, reference):
        """Add nanoraw labels to signal_label handle"""
        events, corr_start_rel_to_raw = self.get_corrected_events()
        events["start"] += corr_start_rel_to_raw
        sequence = ''.join([bytes.decode(x) for x in events['base']])
        hit = get_minimap_alignment(reference, sequence, preset='map-ont')
        cigar_label = np.zeros(len(sequence), dtype=[('raw_start', int), ('raw_length', int), ('reference_index', int),
                                                     ('posterior_probability', float), ('kmer', 'S5')])
        # assign labels
        cigar_label['raw_start'] = events["start"]
        cigar_label['raw_length'] = events["length"]
        if hit.strand:
            reference_map = list(range(hit.r_st, hit.r_en))[::-1]
        else:
            reference_map = list(range(hit.r_st, hit.r_en))
        cigar_label['reference_index'] = reference_map
        cigar_label['kmer'] = events['base']
        cigar_label['posterior_probability'] = [1 for _ in range(hit.r_st, hit.r_en)]
        self.aligned_signal.add_label(cigar_label, name="nanoraw", label_type='label')
        return True

    def add_eventalign_labels(self):
        """Add eventalign labels"""
        section = "template"
        ea_events = self.get_eventalign_events(section=section)
        events = self.get_basecall_data(section=section)
        sampling_freq = self.sample_rate
        start_time = self.raw_attributes['start_time']
        events = time_to_index(events, sampling_freq=sampling_freq, start_time=start_time)
        lables = match_events_with_eventalign(events=ea_events, event_detections=events)
        if self.rna:
            lables["reference_index"] -= self.kmer_index
        else:
            lables["reference_index"] += self.kmer_index

        self.aligned_signal.add_label(lables, name='eventAlign', label_type='label')


def create_labels_from_guide_alignment(events, sam_string, rna=False, reference_path=None, kmer_index=2,
                                       one_ref_indexing=False):
    """Create labeled signal from a guide alignment with only matches being reported

    :param events: path to fast5 file
    :param sam_string: sam alignment string
    :param rna: if read is rna, reverse again
    :param reference_path: if sam_string has MDZ field the reference sequence can be inferred, otherwise, it is needed
    :param kmer_index: index of the kmer to select for reference to event mapping
    :param one_ref_indexing: boolean zero or 1 based indexing for reference
    """
    # test if the required fields are in structured numpy array
    try:
        check_numpy_table(events, req_fields=('raw_start', 'model_state', 'p_model_state', 'raw_length', 'move'))
    except IndexError:
        assert basecall_events["start"].dtype is np.dtype('uint64'), "Event 'start' should be np.int32 type: {}" \
            .format(basecall_events["start"].dtype)
        events = time_to_index(events, sampling_freq=self.sample_rate, start_time=self.raw_attributes['start_time'])

    assert type(one_ref_indexing) is bool, "one_ref_indexing must be a boolean"

    psam_h = initialize_aligned_segment_wrapper(sam_string, reference_path=reference_path)

    # create an indexed map of the events and their corresponding bases
    bases, base_raw_starts, base_raw_lengths, probs = index_bases_from_events(events, kmer_index=kmer_index)

    # check if string mapped to reverse strand
    if psam_h.alignment_segment.is_reverse:
        probs = probs[::-1]
        base_raw_starts = base_raw_starts[::-1]
        # rna reads go 3' to 5' so we dont need to reverse if it mapped to reverse strand
        if not rna:
            bases = ReverseComplement().reverse(''.join(bases))
    # reverse if it mapped to forward strand and RNA
    elif rna:
        bases = ReverseComplement().reverse(''.join(bases))

    # all 'matches' and 'mismatches'
    matches_map = psam_h.seq_alignment.matches_map
    # zero indexed reference start
    ref_start = psam_h.alignment_segment.reference_start + one_ref_indexing
    # set labels
    raw_start = []
    raw_length = []
    reference_index = []
    kmer = []
    posterior_probability = []
    cigar_labels = []
    prev = matches_map[0].reference_index
    for i, alignment in enumerate(matches_map):
        if i == 0 or alignment.reference_index == prev + 1:
            raw_start.append(base_raw_starts[alignment.query_index])
            raw_length.append(base_raw_lengths[alignment.query_index])
            reference_index.append(alignment.reference_index + ref_start)
            kmer.append(alignment.reference_base)
            posterior_probability.append(probs[alignment.query_index])
        else:
            # initialize labels
            cigar_label = np.zeros(len(raw_start),
                                   dtype=[('raw_start', int), ('raw_length', int), ('reference_index', int),
                                          ('posterior_probability', float), ('kmer', 'S5')])
            # assign labels
            cigar_label['raw_start'] = raw_start
            cigar_label['raw_length'] = raw_length
            cigar_label['reference_index'] = reference_index
            cigar_label['kmer'] = kmer
            cigar_label['posterior_probability'] = posterior_probability
            # add to other blocks
            cigar_labels.append(cigar_label)
            # reset trackers
            raw_start = [base_raw_starts[alignment.query_index]]
            raw_length = [base_raw_lengths[alignment.query_index]]
            reference_index = [alignment.reference_index + ref_start]
            kmer = [alignment.reference_base]
            posterior_probability = [probs[alignment.query_index]]
        # keep track of reference positions
        prev = alignment.reference_index

    # catch the last label
    cigar_label = np.zeros(len(raw_start), dtype=[('raw_start', int), ('raw_length', int), ('reference_index', int),
                                                  ('posterior_probability', float), ('kmer', 'S5')])
    # assign labels
    cigar_label['raw_start'] = raw_start
    cigar_label['raw_length'] = raw_length
    cigar_label['reference_index'] = reference_index
    cigar_label['kmer'] = kmer
    cigar_label['posterior_probability'] = posterior_probability
    # add to other blocks
    cigar_labels.append(cigar_label)

    return cigar_labels


def index_bases_from_events(events, kmer_index=2):
    """Map basecalled sequence to events from a table with required fields

    :param kmer_index: index of kmer to create map
    :param events: original base-called events with required fields
    """

    check_numpy_table(events, req_fields=('raw_start', 'model_state', 'p_model_state', 'raw_length', 'move'))
    assert len(events[0]['model_state']) > kmer_index, \
        "Selected too big of a kmer_index len(kmer) !> kmer_index, {} !> {} ".format(len(events[0]['model_state']),
                                                                                     kmer_index)
    probs = []
    base_raw_starts = []
    bases = []
    base_raw_lengths = []
    for i, event in enumerate(events):
        if i == 0:
            # initialize with first kmer
            base_raw_starts.extend([event['raw_start'] for _ in event['model_state'][:kmer_index + 1]])
            probs.extend([event['p_model_state'] for _ in event['model_state'][:kmer_index + 1]])
            bases.extend([chr(x) for x in event['model_state'][:kmer_index + 1]])
            base_raw_lengths.extend([event['raw_length'] for _ in event['model_state'][:kmer_index + 1]])
        else:
            # if there was a move, gather the information for each base by index
            if event['move'] > 0:
                char_moves = bytes.decode(event['model_state'][kmer_index:kmer_index + event['move']])
                for x in range(event['move']):
                    try:
                        base_raw_starts.append(event['raw_start'])
                        probs.append(event['p_model_state'])
                        bases.append(char_moves[x])
                        base_raw_lengths.append(event['raw_length'])
                    except IndexError:
                        print(event["model_state"])
                        print(event['move'])
                        print(char_moves)
                        print(kmer_index)
    # gather last bases for the last event
    base_raw_starts.extend([event['raw_start'] for _ in event['model_state'][kmer_index + 1:]])
    probs.extend([event['p_model_state'] for _ in event['model_state'][kmer_index + 1:]])
    bases.extend([chr(x) for x in event['model_state'][kmer_index + 1:]])
    base_raw_lengths.extend([event['raw_length'] for _ in event['model_state'][kmer_index + 1:]])

    # the index of each corresponds to the index of the final sequence
    return bases, base_raw_starts, base_raw_lengths, probs


def get_eventalign_events(fast5_dir, reference, output_dir, threads=1, overwrite=False):
    """Get nanopolish eventalign events"""
    eventalign_output_path, eventalign_fofn_path = call_eventalign_script(fast5_dir,
                                                                          reference,
                                                                          output_dir,
                                                                          threads=threads,
                                                                          overwrite=overwrite)
    fast5_files = []
    with open(eventalign_fofn_path, 'r') as fofn:
        for line in fofn:
            fast5_files.append(line.split('\t')[1][:-1])
    # gather event data
    dtype = [('contig', 'S10'), ('position', int),
             ('reference_kmer', 'S6'), ('read_index', int),
             ('strand', 'S1'), ('event_index', int),
             ('event_level_mean', float), ('event_stdv', float),
             ('event_length', float), ('model_kmer', 'S6'),
             ('model_mean', float), ('model_stdv', float),
             ('standardized_level', float)]
    with open(eventalign_output_path, 'r') as event_align:
        read_number = 0
        eventalign_data_template = []
        eventalign_data_complement = []
        event_align.readline()
        for line in event_align:
            data = line.split('\t')

            if int(data[3]) != read_number:
                t = np.array(eventalign_data_template, dtype=dtype)
                c = np.array(eventalign_data_complement, dtype=dtype)
                yield t, c, fast5_files[read_number]
                read_number = int(data[3])
                eventalign_data_template = []
                eventalign_data_complement = []

            data[1] = int(data[1])
            data[3] = int(data[3])
            data[5] = int(data[5])
            data[6] = float(data[6])
            data[7] = float(data[7])
            data[8] = float(data[8])
            data[10] = float(data[10])
            data[11] = float(data[11])
            data[12] = float(data[12])
            if str(data[4]) == 't':
                eventalign_data_template.append(tuple(data))
            else:
                eventalign_data_complement.append(tuple(data))

            # print(int(data[3]), read_number)
        t = np.array(eventalign_data_template, dtype=dtype)
        c = np.array(eventalign_data_complement, dtype=dtype)
        assert t or c, "Check reference genome, no alignment generated for any read: {}".format(reference)
        yield t, c, fast5_files[read_number]


def call_eventalign_script(fast5_dir, reference, output_dir, threads=1, overwrite=False):
    """Call eventalign script from scripts folder"""
    # call_eventalign.sh -f ../test_files/minion-reads/canonical/ -t 1 -r ~/data/example_references/ecoli_k12_mg1655.fa -o ~/CLionProjects/nanopolish/dnacanonical/
    call_eventalign_exe = "call_eventalign.sh"
    eventalign_output_path = os.path.join(output_dir, "eventalign.txt")
    eventalign_fofn_path = os.path.join(output_dir, "all_files.fastq.index.readdb")
    if not os.path.exists(eventalign_output_path) or overwrite:
        subprocess.call([call_eventalign_exe, '-f', fast5_dir, '-t', str(threads), '-r', reference, '-o', output_dir])
    return eventalign_output_path, eventalign_fofn_path


def embed_eventalign_events(fast5_dir, reference, output_dir, threads=1, overwrite=False):
    """Call eventalign and embed events"""
    event_generator = get_eventalign_events(fast5_dir,
                                            reference,
                                            output_dir,
                                            threads=threads,
                                            overwrite=overwrite)
    attributes = None
    for template, complement, fast5path in event_generator:
        print(fast5path)
        print("template", template)
        if template or complement:
            handle = Fast5(fast5path, read='r+')
            handle.set_eventalign_table(template=template, complement=complement, meta=attributes, overwrite=True)
        else:
            print("{} did not align".format(fast5path))
    return True


def match_events_with_eventalign(events=None, event_detections=None, minus=False, rna=False):
    """Match event index with event detection data to label segments of signal for each kmer

    # RNA is sequenced 3'-5'
    # reversed for fasta/q sequence
    # if mapped to reverse strand
    # reverse reverse complement = complement

    # DNA is sequenced 5'-3'
    # if mapped to reverse strand
    # reverse complement

    :param events: events table reference_index', 'event_index', 'aligned_kmer', 'posterior_probability
    :param event_detections: event detection event table
    :param minus: boolean option to for minus strand mapping
    :param rna: boolean for RNA read
    """
    assert events is not None, "Must pass signal alignment events"
    assert event_detections is not None, "Must pass event_detections events"

    check_numpy_table(events, req_fields=('position', 'event_index',
                                            'reference_kmer'))

    check_numpy_table(event_detections, req_fields=('start', 'length'))

    label = np.zeros(len(events), dtype=[('raw_start', int), ('raw_length', int), ('reference_index', int),
                                            ('posterior_probability', float), ('kmer', 'S6')])

    label['raw_start'] = [event_detections[x]["start"] for x in events["event_index"]]
    label['raw_length'] = [event_detections[x]["length"] for x in events["event_index"]]
    label['reference_index'] = events["position"]

    def convert_to_str(string):
        """Helper function to catch bytes as strings"""
        if type(string) is str:
            return string
        else:
            return bytes.decode(string)

    flip = ReverseComplement()
    if minus:
        if rna:
            kmers = [flip.complement(convert_to_str(x)) for x in events["reference_kmer"]]
        else:
            kmers = [flip.reverse_complement(convert_to_str(x)) for x in events["reference_kmer"]]
    else:
        if rna:
            kmers = [flip.reverse(convert_to_str(x)) for x in events["reference_kmer"]]
        else:
            kmers = events["reference_kmer"]
    label['kmer'] = kmers
    label['posterior_probability'] = np.ones(len(events))
    # np.sort(label, order='raw_start', kind='mergesort')

    return label


def main():
    """Main docstring"""
    start = timer()



    # sam = "/Users/andrewbailey/CLionProjects/nanopore-RNN/signalAlign/bin/test_output/tempFiles_alignment/tempFiles_miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_strand/temp_sam_file_5048dffc-a463-4d84-bd3b-90ca183f488a.sam"\
    # dna_read = "/Users/andrewbailey/CLionProjects/nanopore-RNN/test_files/minion-reads/canonical/miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch100_read280_strand.fast5"
    # dna_read = "/Users/andrewbailey/CLionProjects/nanopore-RNN/nanotensor/tests/test_files/minion-reads/canonical/miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch100_read280_strand.fast5"
    # dna_read2 = "/Users/andrewbailey/CLionProjects/nanopore-RNN/test_files/minion-reads/canonical/miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_strand.fast5"
    # dna_read3 = "/Users/andrewbailey/CLionProjects/nanopore-RNN/test_files/minion-reads/canonical/over_run/miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_strand.fast5"
    # dna_read4 = "/Users/andrewbailey/CLionProjects/nanopore-RNN/test_files/minion-reads/canonical/consortium_r94_human_dna/rel3-fast5-chr1.part04/DEAMERNANOPORE_20161206_FNFAB49164_MN16450_sequencing_run_MA_821_R9_4_NA12878_12_06_16_71094_ch190_read404_strand.fast5"
    # reference = "/Users/andrewbailey/CLionProjects/nanopore-RNN/test_files/reference-sequences/ecoli_k12_mg1655.fa"
    # reference2 = "/Users/andrewbailey/CLionProjects/nanopore-RNN/test_files/reference-sequences/fake_rna.fa"
    # out_ref = "/Users/andrewbailey/CLionProjects/nanopore-RNN/test_files/reference-sequences/fake_rna_reversed.fa"

    # ReverseComplement().convert_write_fasta(reference2, out_ref, complement=False, reverse=True)
    # rh = ReferenceHandler(reference)
    # seq = rh.get_sequence(chromosome_name="Chromosome", start=623200, stop=623216)
    # print(seq)
    # print("CCACGGGTCCGTCTGG")
    # print("Reference")
    # print(seq)
    # print(ReverseComplement().complement(seq))
    # print("Query")
    # print("CCACGGGTCCGTCTGG")
    # print(ReverseComplement().complement("CCACGGGTCCGTCTGG"))
    # seq = rh.get_sequence(chromosome_name="Chromosome", start=623200-5, stop=623200)
    # print(seq)

    # fast5_dir = "/Users/andrewbailey/CLionProjects/nanopore-RNN/nanotensor/tests/test_files/minion-reads/canonical/"
    # output_dir = "/Users/andrewbailey/data/test_event_align_output"
    # fast5_dir2 = "/Users/andrewbailey/CLionProjects/nanopore-RNN/test_files/minion-reads/canonical/consortium_r94_human_dna/rel3-fast5-chr1.part04/"
    #
    # embed_eventalign_events(fast5_dir2, reference, output_dir, threads=1, overwrite=True)
    #
    # f5handle = Fast5(dna_read4)
    # section = "template"
    # ea_events = f5handle.get_eventalign_events(section=section)
    # print(ea_events.dtype)
    # print(ea_events)
    #
    # events = f5handle.get_basecall_data(section=section)
    # sampling_freq = f5handle.sample_rate
    # start_time = f5handle.raw_attributes['start_time']
    #
    # events = time_to_index(events, sampling_freq=sampling_freq, start_time=start_time)
    # lables = match_events_with_eventalign(events=ea_events, event_detections=events)
    # handle = Fast5(dna_read2)
    # events, corr_start_rel_to_raw = handle.get_corrected_events()
    # events["start"] += corr_start_rel_to_raw
    # sequence = ''.join([bytes.decode(x) for x in events['base']])
    # hit = get_minimap_alignment(reference, sequence, preset='map-ont')
    # print(hit.r_st)
    # print(hit.r_en)
    # print(hit.strand)
    # if str(hit.cigar_str) == str(len(sequence))+'M':
    #     print("Yes")
    # print(events.dtype)
    # cigar_label = np.zeros(len(sequence), dtype=[('raw_start', int), ('raw_length', int), ('reference_index', int),
    #                                               ('posterior_probability', float), ('kmer', 'S5')])
    # # assign labels
    # cigar_label['raw_start'] = events["start"]
    # cigar_label['raw_length'] = events["length"]
    # cigar_label['reference_index'] = list(range(hit.r_st, hit.r_en))
    # cigar_label['kmer'] = events['base']
    # cigar_label['posterior_probability'] = [1 for _ in range(hit.r_st, hit.r_en)]
    #
    # print(cigar_label)
    # events = handle.get_resegment_basecall()

    # kmer_length = 5
    # def make_map(events):
    #     event_map = [0]
    #     previous_prob = 0
    #     for i, line in islice(enumerate(events), 1, None):
    #         print(i, line)
    #         move = line['move']
    #         this_prob = line['p_model_state']
    #         if move == 1:
    #             event_map.append(i)
    #         if move > 1:
    #             for skip in range(move - 1):
    #                 event_map.append(i - 1)
    #             event_map.append(i)
    #         if move == 0:
    #             if this_prob > previous_prob:
    #                 event_map[-1] = i
    #         previous_prob = this_prob
    #     final_event_index = [event_map[-1]]
    #     padding = final_event_index * (kmer_length - 1)
    #     event_map = event_map + padding
    #     return event_map
    # # print(events)
    # print(len(seq))
    # print(len(make_map(events)))

    stop = timer()
    print("Running Time = {} seconds".format(stop - start), file=sys.stderr)


if __name__ == "__main__":
    main()
    raise SystemExit
