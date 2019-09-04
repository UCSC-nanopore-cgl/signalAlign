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
import pandas as pd
from collections import defaultdict
from py3helpers.utils import check_numpy_table, merge_lists
from py3helpers.seq_tools import ReverseComplement, initialize_aligned_segment_wrapper
from signalalign.fast5 import Fast5
from signalalign.mea_algorithm import create_label_from_events
from signalalign.event_detection import add_raw_start_and_raw_length_to_events
from signalalign.variantCaller import MarginalizeVariants


class AlignedSignal(object):
    """Labeled nanopore signal data"""

    def __init__(self, scaled_signal, rna=False):
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
        # place to store event starts
        self.raw_starts = None
        self.variant_calls = defaultdict()
        self.rna = rna

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

    def add_raw_starts(self, raw_starts):
        """Add event start indices to the class"""
        self.raw_starts = np.asarray(raw_starts)

    def add_label(self, label, name, label_type, guide_name=None, check_strand=True):
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
            "posterior_probability must be between zero and one {}".format(max(label["posterior_probability"]))
        if label_type == 'guide':
            assert guide_name is not None, "If label_type is 'guide', you must pass in a guide_name"
        # make sure last label can actually index the signal correctly
        try:
            self.scaled_signal[label[-1]["raw_start"]:label[-1]["raw_start"] + label[-1]["raw_length"]]
        except IndexError:
            raise IndexError("labels are longer than signal")

        label1 = np.sort(label, order=['raw_start'], kind='mergesort')
        if check_strand:
            self.check_strand_mapping(label1)
        # set label with the specified name
        if label_type == 'label':
            self.label[name] = label1
        elif label_type == 'prediction':
            self.prediction[name] = label1
        elif label_type == 'guide':
            self.guide[guide_name][name] = label1

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

    def check_strand_mapping(self, data):
        """Check to see if alignment to reverse strand

        :param data: numpy table with 'reference_index' field
        """
        check_numpy_table(data, req_fields=('reference_index', 'raw_start'))
        # infer strand alignment of read
        if data[0]["reference_index"] >= data[-1]["reference_index"]:
            minus_strand = True
        else:
            minus_strand = False

        # rna is read 3' - 5'
        if self.rna:
            minus_strand = not minus_strand

        if self.minus_strand is not None:
            if data[0]["raw_start"] != data[-1]["raw_start"]:
                assert self.minus_strand == minus_strand, "New label has different strand direction, check label"
        else:
            self.minus_strand = minus_strand

    def add_variant_call(self, variants, name):
        """Add variant call data to the class"""
        self.variant_calls[name] = variants
        return 0


class CreateLabels(Fast5):
    """Create an Aligned Signal object from a fast5 file with """

    def __init__(self, fast5_path, kmer_index=2, rna=False):
        """Initialize fast5 object and keep track of AlignedSignal object
        :param fast5_path: path to fast5 file
        :param kmer_index: nuc index of kmer for mapping
        :param rna: boolean option
        """
        self.fast5_path = fast5_path
        super(CreateLabels, self).__init__(fast5_path)
        self.kmer_index = kmer_index
        self.rna = self.is_read_rna() or rna
        self.aligned_signal = self._initialize()
        self.has_guide_alignment = False

    def _initialize(self):
        """Initialize AlignedSignal class by adding the raw and scaled signal"""
        scaled_signal = self.get_read(raw=True, scale=True)
        raw_signal = self.get_read(raw=True, scale=False)
        # add raw signal information to AlignedSignal
        aligned_signal = AlignedSignal(scaled_signal, rna=self.rna)
        aligned_signal.add_raw_signal(raw_signal)
        return aligned_signal

    def add_variant_data(self, number=None, complement=False):
        """Add variant data to aligned signal"""
        if number is not None:
            assert type(number) is int, "Number must be an integer"
            path = self.__default_signalalign_events__.format(number)
            variant_data = self.get_signalalign_events(variant=True, override_path=path, complement=complement)
            # print("mea_alignment path: {}".format(path))
            name = "variantCaller_{}_{}".format(number, "complement" if complement else "template")
        else:
            variant_data = self.get_signalalign_events(variant=True)
            name = "variantCaller"

        mv_h = MarginalizeVariants(variant_data, variants=[x.decode() for x in set(variant_data["base"])],
                                   read_name=self.get_read_id())
        data = mv_h.get_data()
        if complement:
            data = data[data["strand"] == "c"]
        else:
            data = data[data["strand"] == "t"]

        final_data = match_ref_position_with_raw_start_band(data, variant_data)

        self.aligned_signal.add_variant_call(final_data, name=name)
        return final_data

    def fix_sa_reference_indexes(self, data):
        """Fix reference indexes based on kmer length and kmer index"""
        check_numpy_table(data, req_fields=('reference_index', 'raw_start', 'kmer'))
        kmer_len = len(data[0]["kmer"])

        if self.aligned_signal.minus_strand:
            data["reference_index"] += ((kmer_len - 1) - self.kmer_index)
        else:
            data["reference_index"] += self.kmer_index
        return data

    def add_mea_labels(self, number=None, complement=False):
        """Gather mea_alignment labels information from fast5 file.
        :param number: integer representing which signal align predictions to plot
        :param complement: option to look for mea_complement data
        """
        if number is not None:
            assert type(number) is int, "Number must be an integer"
            path = self.__default_signalalign_events__.format(number)
            mea_alignment = self.get_signalalign_events(mea=True, override_path=path, complement=complement)
            # print("mea_alignment path: {}".format(path))

            name = "mea_signalalign_{}_{}".format(number, "complement" if complement else "template")
        else:
            mea_alignment = self.get_signalalign_events(mea=True)
            name = "mea_signalalign"
        # rna reference positions are on 5' edge aka right side of kmer
        if self.aligned_signal.minus_strand is None:
            self.aligned_signal.check_strand_mapping(mea_alignment)

        mea_alignment = self.fix_sa_reference_indexes(mea_alignment)

        self.aligned_signal.add_label(mea_alignment, name=name, label_type='label', check_strand=not complement)
        return mea_alignment

    def add_signal_align_predictions(self, number=None, add_basecall=False, two_d=False):
        """Create prediction using probabilities from full output format from signalAlign
        :param number: integer representing which signal align predictions to plot
        :param add_basecall: if set to true, will add basecalled event table
        :param two_d: boolean option to visualize 2d reads
        """
        if number is not None:
            assert type(number) is int, "Number must be an integer"
            path = self.__default_signalalign_events__.format(number)
            sa_events = self.get_signalalign_events(override_path=path)
            basecall_events_path = self.get_signalalign_basecall_path(override_path=path)
            # print("sa_events path: {}".format(path))
            # print("basecall_events_path: {}".format(basecall_events_path))
            sam = self.get_signalalign_events(sam=True, override_path=path)
            name = "full_signalalign_{}".format(number)
        else:
            sa_events = self.get_signalalign_events()
            basecall_events_path = self.get_signalalign_basecall_path()
            name = "full_signalalign"
            sam = None
        # cut out duplicates
        if add_basecall:
            try:
                events = self[basecall_events_path][()]
            except:
                raise ValueError('Could not retrieve basecall_1D data from events')
            self.add_basecall_alignment_prediction(my_events=events,
                                                   number=basecall_events_path[24],
                                                   sam=sam)
        if not two_d:
            sa_events = sa_events[sa_events["strand"] == b"t"]
        predictions = create_label_from_events(sa_events)

        predictions = self.fix_sa_reference_indexes(predictions)

        self.aligned_signal.add_label(predictions, name=name, label_type='prediction', check_strand=not two_d)
        return predictions

    def add_basecall_alignment_prediction(self, sam=None, number=None, add_mismatches=False, my_events=None, trim=None):
        """Add the original basecalled event table and add matches and missmatches to 'prediction'
            alignment labels to signal_label handle
        :param sam: correctly formatted SAM string
        :param number: integer representing which signal align predictions to plot
        :param add_mismatches: boolean option to add mismatches to labels
        :param my_events: if you want to pass the event table in directly
        :param trim: trim both sides of continuous matches to show anchor pairs
        :return: True if the correct labels are added to the AlignedSignal internal class object
        """
        if not sam:
            sam = self.get_signalalign_events(sam=True)
        matches_name = "matches_guide_alignment"
        mismatches_name = "mismatches_guide_alignment"

        if my_events is not None:
            events = my_events
            if number is not None:
                matches_name = "matches_guide_alignment_{}".format(number)
                mismatches_name = "mismatches_guide_alignment_{}".format(number)

        else:
            if number is not None:
                events = self.get_basecalled_data_by_number(number)
                matches_name = "matches_guide_alignment_{}".format(number)
                mismatches_name = "mismatches_guide_alignment_{}".format(number)
            else:
                events = self.get_basecall_data()
        try:
            check_numpy_table(events, req_fields=('raw_start', 'raw_length'))
        # if events do not have raw_start or raw_lengths
        except KeyError:
            events = add_raw_start_and_raw_length_to_events(events, self.sample_rate, self.raw_attributes["start_time"])

        matches, mismatches, raw_starts = match_cigar_with_basecall_guide(events=events, sam_string=sam, rna=self.rna,
                                                                          kmer_index=self.kmer_index)

        if trim is not None:
            matches = trim_matches(matches, trim=trim)

        self.aligned_signal.add_label(matches, name=matches_name, label_type='prediction')
        if add_mismatches:
            self.aligned_signal.add_label(mismatches, name=mismatches_name, label_type='prediction')
        self.aligned_signal.add_raw_starts(raw_starts)
        self.has_guide_alignment = True
        return matches, mismatches

    def get_basecalled_data_by_number(self, number):
        """Get basecalled event data by the number"""
        assert type(number) is int, "Number must be an integer"
        events_path = self.__default_template_1d_basecall_events__.format(number)
        try:
            events = self[events_path][()]
        except:
            raise ValueError('Could not retrieve basecall_1D data from {}'.format(events_path))
        try:
            check_numpy_table(events, req_fields=('raw_start', 'raw_length'))
        # if events do not have raw_start or raw_lengths
        except KeyError:
            events = add_raw_start_and_raw_length_to_events(events, self.sample_rate, self.raw_attributes["start_time"])
        return events

    def add_tombo_labels(self, number=0):
        """Add nanoraw labels to signal_label handle"""
        events, corr_start_rel_to_raw = self.get_corrected_events(number)
        attributes = self.get_corrected_events_attr()
        events["start"] += corr_start_rel_to_raw
        sequence = ''.join([bytes.decode(x) for x in events['base']])
        cigar_label = np.zeros(len(sequence), dtype=[('raw_start', int), ('raw_length', int), ('reference_index', int),
                                                     ('posterior_probability', float), ('kmer', 'S5')])
        # assign labels
        cigar_label['raw_start'] = events["start"]
        cigar_label['raw_length'] = events["length"]
        reference_map = list(range(attributes["mapped_start"], attributes["mapped_end"]))
        if attributes["mapped_strand"] != "+":
            reference_map = reference_map[::-1]
        cigar_label['reference_index'] = reference_map
        cigar_label['kmer'] = events['base']
        cigar_label['posterior_probability'] = [1 for _ in range(len(reference_map))]
        self.aligned_signal.add_label(cigar_label, name="tombo_{}".format(number), label_type='label')
        return cigar_label

    # def add_eventalign_labels(self):
    #     """Add eventalign labels"""
    #     section = "template"
    #     ea_events = self.get_eventalign_events(section=section)
    #     events = self.get_basecall_data(section=section)
    #     sampling_freq = self.sample_rate
    #     start_time = self.raw_attributes['start_time']
    # TODO fix time to index
    #     events = time_to_index(events, sampling_freq=sampling_freq, start_time=start_time)

    #     lables = match_events_with_eventalign(events=ea_events, event_detections=events)
    #     if self.rna:
    #         lables["reference_index"] -= self.kmer_index
    #     else:
    #         lables["reference_index"] += self.kmer_index
    #
    #     self.aligned_signal.add_label(lables, name='eventAlign', label_type='label')


def get_distance_from_guide_alignment(data, guide_data, reference_index_key="position", minus_strand=False):
    """Calculate the distance of input data alignment to the guide alignment.
    :param data: input data with at least "raw_start", "raw_length", and reference_index_key fields
    :param guide_data: guide alignmnet data
    :param reference_index_key: key to grab reference index from data
    :param minus_strand: boolean option if data is aligned to minus strand
    :return: modified data with "guide_delta" field
    """
    variant_data = data.sort_values(by=reference_index_key)
    if minus_strand:
        guide_data = guide_data[::-1]

    distance_to_guide = []
    variant_index = 0
    len_variant_data = len(variant_data)
    v_position = variant_data.iloc[variant_index][reference_index_key]
    for i, guide in enumerate(guide_data.itertuples()):
        if getattr(guide, "reference_index") >= v_position:
            if getattr(guide, "reference_index") == v_position:
                guide_index = i
            else:
                guide_index = i-1

            v_position_middle = (variant_data.iloc[variant_index]["raw_start"] +
                                 (variant_data.iloc[variant_index]["raw_length"] / 2))
            guide_middle_position = np.round((guide_data.iloc[guide_index]["raw_start"] + (guide_data.iloc[guide_index]["raw_length"] / 2)))

            distance_to_guide.append(v_position_middle - guide_middle_position)
            variant_index += 1
            if variant_index < len_variant_data:
                v_position = variant_data.iloc[variant_index][reference_index_key]
            else:
                break

    distance = pd.DataFrame(distance_to_guide, columns=['guide_delta'])
    final_data = pd.concat([variant_data, distance], axis=1)
    return final_data


def match_ref_position_with_raw_start_band(aggregate_reference_position, per_event_data):
    """Match up the reference position from aggregated probability table and the per event data"""
    final_data = []
    for position in aggregate_reference_position["position"]:

        # get the start and length of total number of events which had bases aligned to this position
        pos_data = per_event_data[per_event_data["reference_position"] == position]
        min_raw_start = min(pos_data["raw_start"])
        last_event_length = pos_data[pos_data["raw_start"] == max(pos_data["raw_start"])]["raw_length"][0]
        total_length = max(pos_data["raw_start"]) - min(pos_data["raw_start"]) + last_event_length
        final_data.append(
            merge_lists([aggregate_reference_position[aggregate_reference_position["position"] == position].values.tolist()[0],
                         [min_raw_start, total_length]]))
    final_data = pd.DataFrame(final_data, columns=merge_lists([aggregate_reference_position.columns,
                                                               ["raw_start", "raw_length"]]))
    return final_data


def trim_matches(matches, trim):
    """Trim continuous matches by trim on both sides of the block of matches.

    if trim = 2
        10M -> 6M

    :param matches:

    """
    # mask = np.ones(len(arr), dtype=bool)
    # mask[[0,2,4]] = False
    # result = arr[mask,...]

    # deal with simple case
    if trim == 0:
        return matches
    else:
        mask = np.ones(len(matches), dtype=bool)
        remove_indices = []
        curr_ref_index = matches[0]["reference_index"]-1
        # track where we are in match block
        length = 0
        for i, match in enumerate(matches):
            if match["reference_index"] == curr_ref_index + 1:
                length += 1
                curr_ref_index = match["reference_index"]
                continue

            if length < 2*trim:
                for x in range(i-length, i):
                    remove_indices.append(x)
            else:
                for x in range(i-length, i-length+trim):
                    remove_indices.append(x)
                for x in range(i-trim, i):
                    remove_indices.append(x)
            length = 1
            curr_ref_index = match["reference_index"]

        if length < 2*trim:
            for x in range(i-length+1, i+1):
                remove_indices.append(x)
        else:
            for x in range(i-length+1, i-length+trim+1):
                remove_indices.append(x)
            for x in range(i-trim+1, i+1):
                remove_indices.append(x)

        mask[remove_indices] = False
        matches = matches[mask]
        return matches


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
                        pass
                        # print(event["model_state"])
                        # print(event['move'])
                        # print(char_moves)
                        # print(kmer_index)
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


def match_cigar_with_basecall_guide(events, sam_string, kmer_index, rna=False, reference_path=None,
                                    one_ref_indexing=False):
    """Create labeled signal from a guide alignment with only matches being reported

    :param events: path to fast5 file
    :param sam_string: sam alignment string
    :param rna: if read is rna, reverse again
    :param reference_path: if sam_string has MDZ field the reference sequence can be inferred, otherwise, it is needed
    :param kmer_index: index of the kmer to select for reference to event mapping
    :param one_ref_indexing: boolean zero or 1 based indexing for reference
    """
    check_numpy_table(events, req_fields=('raw_start', 'model_state', 'p_model_state', 'raw_length', 'move'))
    assert type(one_ref_indexing) is bool, "one_ref_indexing must be a boolean"

    psam_h = initialize_aligned_segment_wrapper(sam_string, reference_path=reference_path)

    # create an indexed map of the events and their corresponding bases
    _, base_raw_starts, base_raw_lengths, probs = index_bases_from_events(events, kmer_index=kmer_index)
    if rna:
        # events are 3'-5', swap to correct for alignment file
        base_raw_starts = base_raw_starts[::-1]
        base_raw_lengths = base_raw_lengths[::-1]
        probs = probs[::-1]
    # all 'matches' and 'mismatches'
    matches_map = psam_h.seq_alignment.matches_map
    ref_len = len(psam_h.get_reference_sequence())
    # zero indexed reference start
    ref_start = psam_h.alignment_segment.reference_start + one_ref_indexing
    # set labels
    matches_raw_start = []
    matches_raw_length = []
    matches_reference_index = []
    matches_kmer = []
    matches_posterior_probability = []

    mismatches_raw_start = []
    mismatches_raw_length = []
    mismatches_reference_index = []
    mismatches_kmer = []
    mismatches_posterior_probability = []

    for i, alignment in enumerate(matches_map):
        if alignment.query_base == alignment.reference_base:
            matches_raw_start.append(base_raw_starts[alignment.query_index])
            matches_raw_length.append(base_raw_lengths[alignment.query_index])
            matches_kmer.append(alignment.reference_base)
            matches_posterior_probability.append(probs[alignment.query_index])
            if psam_h.alignment_segment.is_reverse:
                matches_reference_index.append((ref_start + ref_len - 1) - alignment.reference_index)
            else:
                matches_reference_index.append(ref_start + alignment.reference_index)
        else:
            mismatches_raw_start.append(base_raw_starts[alignment.query_index])
            mismatches_raw_length.append(base_raw_lengths[alignment.query_index])
            mismatches_kmer.append(alignment.reference_base)
            mismatches_posterior_probability.append(probs[alignment.query_index])
            if psam_h.alignment_segment.is_reverse:
                mismatches_reference_index.append((ref_start + ref_len - 1) - alignment.reference_index)
            else:
                mismatches_reference_index.append(ref_start + alignment.reference_index)

    matches = np.zeros(len(matches_raw_start), dtype=[('raw_start', int), ('raw_length', int), ('reference_index', int),
                                                      ('posterior_probability', float), ('kmer', 'S5')])
    mismatches = np.zeros(len(mismatches_raw_start),
                          dtype=[('raw_start', int), ('raw_length', int), ('reference_index', int),
                                 ('posterior_probability', float), ('kmer', 'S5')])
    # assign labels
    matches['raw_start'] = matches_raw_start
    matches['raw_length'] = matches_raw_length
    matches['reference_index'] = matches_reference_index
    matches['kmer'] = matches_kmer
    matches['posterior_probability'] = matches_posterior_probability

    mismatches['raw_start'] = mismatches_raw_start
    mismatches['raw_length'] = mismatches_raw_length
    mismatches['reference_index'] = mismatches_reference_index
    mismatches['kmer'] = mismatches_kmer
    mismatches['posterior_probability'] = mismatches_posterior_probability

    # trim extra event alignments
    return matches[kmer_index:len(matches)-kmer_index], mismatches[kmer_index:len(matches)-kmer_index], events["raw_start"]


def get_highest_single_digit_integer_from_string(in_string):
    """Get the highest single digit integer from a string """
    items = set(in_string)
    ints = set("1234567890")
    items &= ints
    max_i = -1
    if len(items) == 0:
        return None
    else:
        for i in items:
            if int(i) > max_i:
                max_i = int(i)
        return max_i


if __name__ == "__main__":
    print("This is a library of functions")
    raise SystemExit
