#!/usr/bin/env python
"""Plot labelled read"""
########################################################################
# File: plot_labelled_read.py
#  executable: plot_labelled_read.py
#
# Author: Andrew Bailey
# History: Created 02/24/18
########################################################################

from __future__ import print_function
import sys
import os
import platform
import colorsys
import matplotlib as mpl

if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
if platform.system() == "Darwin":
    mpl.use("macosx")
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from matplotlib.collections import LineCollection

import numpy as np
import pandas as pd
from collections import defaultdict
from timeit import default_timer as timer
from argparse import ArgumentParser
from signalalign.alignedsignal import AlignedSignal, CreateLabels
from signalalign.nanoporeRead import NanoporeRead
from py3helpers.utils import list_dir, find_substring_indices
from py3helpers.seq_tools import initialize_aligned_segment_wrapper, ReverseComplement

MEA = 'm'
SA_FULL = 's'
BASECALL = 'b'
EVENT_INDEX = 'e'
SA_ALIGNMENT_DIFF = 'd'
ABS_SA_ALIGNMENT_DIFF = 'D'
RAW_START = 'r'


def parse_args():
    parser = ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument('--dir', '-d', action='store',
                       dest='dir', default=None,
                       help="Path to directory of fast5 files to plot")

    group.add_argument('--f5_path', action='store', default=None,
                       dest='f5_path',
                       help="Path to a specific fast5 file to plot")

    parser.add_argument('--basecall', nargs='+',
                        dest='basecall', required=False,
                        help="Option to plot the most recent basecalled data aligned to reference")

    parser.add_argument('--rna', action='store_true', default=False,
                        dest='rna', required=False,
                        help="Option store true if RNA reads")

    parser.add_argument('--basecall_scatter', action='store_true', default=False,
                        dest='basecall_scatter', required=False,
                        help="Plot basecalled segments as scatter plot")

    parser.add_argument('--mea', nargs='+',
                        dest='mea', required=False,
                        help="Option to plot the maximum expected accuracy alignment from signalalign")

    parser.add_argument('--tombo', nargs='+',
                        dest='tombo', required=False,
                        help="Option to plot the alignment from Tombo (RawGenomeCorrected)")

    parser.add_argument('--sa_full', nargs='+',
                        dest='sa_full', required=False,
                        help="Option to plot all of the posterior probabilities from the signalalign output")

    parser.add_argument('--variant', nargs='+',
                        dest='variant', required=False,
                        help="Path to variant data tsv file")

    parser.add_argument('--output_dir', action='store',
                        dest='output_dir', required=False, type=str,
                        help="If set, will write out to output directory otherwise the graph is just shown")

    parser.add_argument('--plot_alpha', action='store_true', default=False, dest="plot_alpha", required=False,
                        help="If set, will plot probability information as the alpha associated with each data point")

    parser.add_argument('--two_d', action='store_true', default=False, dest="two_d", required=False,
                        help="Plot's complement mea")

    parser.add_argument('--trim', action='store', default=0, dest="trim", required=False, type=int,
                        help="If set, will trim runs of matches")

    parser.add_argument('--kmers', action='store', nargs='+', default=None, dest="kmers", required=False, type=str,
                        help="Highlight kmers in reference sequence")

    args = parser.parse_args()
    return args


class PlotSignal(object):
    """Handles alignments between nanopore signals, events and reference sequences"""

    def __init__(self, aligned_signal):
        """Plot alignment

        :param fname: path to fast5 file
        """
        assert isinstance(aligned_signal, AlignedSignal), "aligned_signal must be an instance of AlignedSignal"
        self.signal_h = aligned_signal
        self.names = []
        self.alignments = []
        self.predictions = []
        self.variant_calls = []
        self.event_starts = None
        self.guide_alignments = defaultdict(list)

    def get_alignments(self):
        """Format alignments from AlignedSignal to make it possible to plot easily"""
        for name, label in self.signal_h.label.items():
            self.names.append(name)
            self.alignments.append([label['raw_start'], label['reference_index']])

        # predictions can have multiple alignments for each event and are not a path
        for name, prediction in self.signal_h.prediction.items():
            self.names.append(name)
            self.predictions.append([prediction['raw_start'], prediction['raw_length'], prediction['reference_index'],
                                     prediction['posterior_probability']])
        # guides are continuous paths but are an accumulation of multiple match groups
        for whole_guide_name, guide in self.signal_h.guide.items():
            for _, sub_guide in guide.items():
                self.guide_alignments[whole_guide_name].append([sub_guide['raw_start'], sub_guide['reference_index']])
            # gather tail ends of alignments
            self.names.append(whole_guide_name)

        # variant calling data
        for name, v_c in self.signal_h.variant_calls.items():
            self.names.append(name)
            self.variant_calls.append(v_c)


        # add raw starts
        if self.signal_h.raw_starts is not None:
            self.event_starts = self.signal_h.raw_starts

    def plot_alignment(self, save_fig_path=None, plot_alpha=False, plot_lines=False, kmer_info=None):
        """Plot the alignment between events and reference with the guide alignment and mea alignment
        :param kmer_info: plot kmer data if provided.

                    data structure : [[[ref_start, length], [ref_start2, length2] ], {pos: base, pos2: base2}]]

        :param plot_lines: boolean option, if set to true, lines will be plotted instead of points
        :param save_fig_path: if set, will write image to path
        :param plot_alpha: boolean:
                        if set to false will ignore probability information associated with data points
        """
        if save_fig_path:
            assert os.path.exists(os.path.dirname(save_fig_path)), \
                "Output directory does not exist: {}".format(save_fig_path)

        self.get_alignments()

        plt.figure(figsize=(6, 8))
        panel1 = plt.axes([0.1, 0.22, .8, .7])
        panel2 = plt.axes([0.1, 0.09, .8, .1], sharex=panel1)
        panel1.set_xlabel('Events')
        panel1.set_ylabel('Reference')
        panel1.xaxis.set_label_position('top')
        panel1.invert_yaxis()
        panel1.xaxis.tick_top()
        panel1.grid(color='black', linestyle='-', linewidth=1)

        handles = list()
        colors2 = ['blue', 'green', 'red', 'yellow', 'magenta', 'deepskyblue', 'purple', 'lime']
        colors = [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                  [1.0, 0.0, 1.0], [1.0, 0.5, 0.0], [0.0, 1.0, 1.0], [0.0, 0.0, 0.0]]

        # colors = get_spaced_colors(len(self.names) + 1)
        color_selection = 0
        # plot signal alignments
        for i, alignment in enumerate(self.alignments):
            handle, = panel1.plot(alignment[0], alignment[1] - (i * 0.5), color=colors[color_selection], alpha=1)
            color_selection += 1
            handles.append(handle)

        # plot predictions (signal align outputs)
        for i, prediction in enumerate(self.predictions):
            rgba_colors = np.tile(np.array(colors[color_selection]), (len(prediction[0]), 1))
            if plot_alpha and "full_signalalign" in self.names[color_selection]:
                rgba_colors = np.insert(rgba_colors, 3, prediction[3].tolist(), axis=1)
            else:
                rgba_colors = np.insert(rgba_colors, 3, [1 for _ in range(len(prediction[3].tolist()))], axis=1)

            # rgba_colors[:, 3] = prediction[2].tolist()
            if plot_lines:
                lines = list(zip(list(zip(prediction[0].tolist(), [x + (i * 0.05) for x in prediction[2].tolist()])),
                                 list(zip([x + y for x, y in zip(prediction[0].tolist(), prediction[1].tolist())],
                                          [x + (i * .05) for x in prediction[2].tolist()]))))
                line_segments = LineCollection(lines,
                                               linestyles='solid',
                                               colors=rgba_colors)
                handle = panel1.add_collection(line_segments)
            else:
                handle = panel1.scatter(prediction[0].tolist(), [x + (i * .1) for x in prediction[2].tolist()],
                                        marker='.',
                                        c=rgba_colors)
            handles.append(handle)
            color_selection += 1

        # plot original basecalled events
        if self.guide_alignments:
            for i, (name, guide) in enumerate(self.guide_alignments.items()):
                for j, sub_guide in enumerate(guide):
                    handle, = panel1.plot(sub_guide[0], sub_guide[1], c=colors[color_selection])
                handles.append(handle)
                color_selection += 1

        colors2_index = -1
        if self.variant_calls is not []:
            for i, data in enumerate(self.variant_calls):
                variants = set(data.columns) - {'forward_mapped', 'raw_start', 'position', 'raw_length', 'read_name', 'contig', 'strand'}
                colors2_index += 1
                color = colors2[colors2_index]
                for j, position in data.iterrows():
                    rectangle = plt.Rectangle((position["raw_start"], position["position"]-5), position["raw_length"],
                                              11, alpha=0.4, color=color)

                    panel1.add_patch(rectangle)
                    variant_string = ""
                    for variant_base in variants:
                        variant_string += "p({} = {}) = {}\n".format(position["position"], variant_base,
                                                                     np.round(position[variant_base], 3))

                    panel1.text(position["raw_start"] + 10, position["position"]-6, "{}".format(variant_string))

        # if kmers provided plot them
        if kmer_info:
            kmer_blocks = kmer_info[0]
            positions = kmer_info[1]
            for i, kmer_block in enumerate(kmer_blocks):
                color = colors2[colors2_index]
                colors2_index += 1
                for ref_start, length in kmer_block:
                    rectangle = plt.Rectangle((0, ref_start), len(self.signal_h.scaled_signal), length,
                                              alpha=0.4, color=color)
                    panel1.add_patch(rectangle)

            def format_func(value, tick_number):
                # find number of multiples of pi/2
                if positions[value] != 0:
                    return positions[value]
                else:
                    return value

            panel1.yaxis.set_major_formatter(plt.FuncFormatter(format_func))



        panel2.set_xlabel('Time')
        panel2.set_ylabel('Current (pA)')
        # if self.event_starts is not None:
        #     panel2.set_xticks(self.event_starts, minor=True)

        panel2.axvline(linewidth=0.1, color="k")

        handle, = panel2.plot(self.signal_h.scaled_signal, color="black", lw=0.4)
        handles.append(handle)
        self.names.append("scaled_signal")

        box = panel1.get_position()
        panel1.set_position([box.x0, box.y0 + box.height * 0.1,
                             box.width, box.height * 0.9])

        # Put a legend below current axis
        panel1.legend(handles, self.names, loc='upper center', bbox_to_anchor=(0.5, -0.05),
                      fancybox=True, shadow=True, ncol=5)

        # panel1.legend(handles, self.names, loc='upper right')
        if save_fig_path:
            plt.savefig(save_fig_path)
        else:
            plt.show()
            pass


def plot_basecalled_sequences(events1, events2=None, signal=None, save_fig_path=False):
    """Plot the basecalled event table with the kmers on the y axis and signal on x axis"""
    #
    handles = list()
    names = []
    plt.figure(figsize=(6, 8))
    panel1 = plt.axes([0.1, 0.22, .8, .7])
    panel2 = plt.axes([0.1, 0.09, .8, .1], sharex=panel1)
    panel1.set_xlabel('Events')
    panel1.set_ylabel('Basecalled Sequence')
    panel1.xaxis.set_label_position('top')
    panel1.invert_yaxis()
    panel1.xaxis.tick_top()
    panel1.grid(color='black', linestyle='-', linewidth=1)

    sequence = NanoporeRead.sequence_from_events(events1)
    sequence2 = NanoporeRead.sequence_from_events(events2)
    assert sequence2 == sequence, "Basecalled sequences are different"
    event_map = NanoporeRead.make_event_map(events1, kmer_length=len(events1[0]["model_state"]))
    # for y, x in enumerate([events1[x]["raw_start"] for x in event_map]):
    #     panel1.text(x, y, sequence[y])

    handle = panel1.scatter([events1[x]["raw_start"] for x in event_map], range(len(event_map)),
                            marker='.',
                            c='blue')
    handles.append(handle)
    names.append("Original Event Table")
    if events2 is not None:
        event_map = NanoporeRead.make_event_map(events2, kmer_length=len(events2[0]["model_state"]))

        handle2 = panel1.scatter([events2[x]["raw_start"] for x in event_map], range(len(event_map)),
                                 marker='.',
                                 c='red')
        handles.append(handle2)
        names.append("Load_from_raw Event Table")

    if signal is not None:
        handle, = panel2.plot(signal, color="black", lw=0.4)

        handles.append(handle)
        names.append("signal")

    panel2.set_xlabel('Time')
    panel2.set_ylabel('Current (pA)')

    panel1.legend(handles, names, loc='upper center', bbox_to_anchor=(0.5, -0.05),
                  fancybox=True, shadow=True, ncol=5)

    if save_fig_path:
        plt.savefig(save_fig_path)
    else:
        plt.show()
        pass


def get_spaced_colors(n):
    """Create n evenly spaced RGB colors
    :param n: number of colors needed

    source: https://stackoverflow.com/questions/876853/generating-color-ranges-in-python
    """
    hsv_tuples = [(x * 1.0 / n, 0.5, 0.5) for x in range(n)]
    rgb_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), hsv_tuples)
    return list(rgb_tuples)


def analyze_event_skips(mea_events, sa_full_events, basecall_events,
                        generate_plot=True):
    # prep
    raw_start = lambda x: x['raw_start'] if 'raw_start' in x else x[0]
    event_rank = lambda x: x['posterior_probability'] if 'posterior_probability' in x else x[3]
    aligned_position = lambda x: x['reference_index'] if 'reference_index' in x else x[2]

    # build datastructure holding all events by raw start
    raw_start_to_event = dict()

    def add_event(event, type):
        key = raw_start(event)
        if key not in raw_start_to_event:
            raw_start_to_event[key] = dict()
        if type in raw_start_to_event[key]:
            if event_rank(raw_start_to_event[key][type]) > event_rank(event):
                return
        raw_start_to_event[key][type] = event

    list(map(lambda x: add_event(x, MEA), mea_events))
    list(map(lambda x: add_event(x, SA_FULL), sa_full_events))
    list(map(lambda x: add_event(x, BASECALL), basecall_events))

    # enumerating events
    all_event_keys = list(raw_start_to_event.keys())
    all_event_keys.sort()
    raw_start_to_event_idx = {s: i for i, s in enumerate(all_event_keys)}
    event_idx_to_raw_start = {i: s for i, s in enumerate(all_event_keys)}
    max_event_idx = max(event_idx_to_raw_start.keys())

    # finding distance to basecall
    sa_event_summary = dict()
    for i, raw_start in enumerate(all_event_keys):
        # skip events only aligned in orig basecall
        if SA_FULL not in raw_start_to_event[raw_start]: continue

        # get sa_full event
        sa_full = raw_start_to_event[raw_start][SA_FULL]

        # aligned distance is easy
        if BASECALL in raw_start_to_event[raw_start]:
            difference = aligned_position(sa_full) - aligned_position(raw_start_to_event[raw_start][BASECALL])

        # if they're not aligned, take the mean reference pos between the prev and next alignment
        else:
            last_aligned_pos = None
            j = i
            while j >= 0:
                j_rs = event_idx_to_raw_start[j]
                if BASECALL in raw_start_to_event[j_rs]:
                    last_aligned_pos = aligned_position(raw_start_to_event[j_rs][BASECALL])
                    break
                j -= 1

            next_aligned_pos = None
            j = i
            while j <= max_event_idx:
                j_rs = event_idx_to_raw_start[j]
                if BASECALL in raw_start_to_event[j_rs]:
                    next_aligned_pos = aligned_position(raw_start_to_event[j_rs][BASECALL])
                    break
                j += 1

            all_reference_pos = list(filter(lambda x: x is not None, [last_aligned_pos, next_aligned_pos]))
            avg_reference_pos = int(np.mean(all_reference_pos))
            difference = aligned_position(sa_full) - avg_reference_pos

        sa_event_summary[raw_start] = {
            SA_FULL: sa_full,
            MEA: MEA in raw_start_to_event[raw_start],
            EVENT_INDEX: raw_start_to_event_idx[raw_start],
            SA_ALIGNMENT_DIFF: difference,
            ABS_SA_ALIGNMENT_DIFF: abs(difference),
            RAW_START: raw_start
        }

    # get all summary info
    summaries = list(sa_event_summary.values())
    summaries.sort(key=lambda x: x[RAW_START])

    # separate based on mea and not mea
    summaries__mea_event_idx = list(map(lambda x: x[EVENT_INDEX], list(filter(lambda x: x[MEA], summaries))))
    summaries__mea_difference = list(map(lambda x: x[SA_ALIGNMENT_DIFF], list(filter(lambda x: x[MEA], summaries))))
    summaries__no_mea_event_idx = list(map(lambda x: x[EVENT_INDEX], list(filter(lambda x: not x[MEA], summaries))))
    summaries__no_mea_difference = list(
        map(lambda x: x[SA_ALIGNMENT_DIFF], list(filter(lambda x: not x[MEA], summaries))))

    # plot it if appropriate
    if generate_plot is not None and generate_plot:
        plt.plot(summaries__mea_event_idx, summaries__mea_difference, color='blue', label='MEA')
        plt.scatter(summaries__no_mea_event_idx, summaries__no_mea_difference, color='red', label='Not MEA', alpha=.5,
                    s=10)
        plt.xlabel("Event Index")
        plt.ylabel("Reference Alignment Difference: (SignalAlignRefPos - CigarStringRefPos)")
        plt.legend()

        if type(generate_plot) == str:
            plt.savefig(generate_plot)
            plt.close()
        else:
            plt.show()

    return summaries


def generate_kmer_info_for_plotting(kmers, aligned_segment):
    """Create the correct data structure to bed fed into the plot_alignment function"""
    ref_seq = aligned_segment.get_reference_sequence()
    ref_start = aligned_segment.alignment_segment.reference_start
    almost_final_dict = defaultdict(int)
    all_blocks = []
    for kmer in kmers:
        blocks = [[x + ref_start, len(kmer)] for x in
                  find_substring_indices(ref_seq,
                                         kmer,
                                         offset=0)]
        all_blocks.append(blocks)
        for block in blocks:
            for i, base in enumerate(kmer):
                almost_final_dict[block[0] + i] = base

    kmer_info = [all_blocks, almost_final_dict]
    return kmer_info


def main(args=None):
    """Plot event to reference labelled ONT nanopore reads """
    start = timer()

    args = args if args is not None else parse_args()

    if args.dir:
        f5_locations = list_dir(args.dir, ext="fast5")
    else:
        f5_locations = [args.f5_path]

    for f5_path in f5_locations:
        try:
            main_plot = True
            save_fig_path = None
            cl_handle = CreateLabels(f5_path, kmer_index=2, rna=args.rna)
            if args.output_dir:
                save_fig_path = "{}.png".format(os.path.join(args.output_dir,
                                                             os.path.splitext(os.path.basename(f5_path))[0]))

            mea_list = list()
            if args.mea:
                for number in args.mea:
                    mea = cl_handle.add_mea_labels(number=int(number))
                    mea_list.append(mea)
                    if args.two_d:
                        mea = cl_handle.add_mea_labels(number=int(number), complement=True)
                        mea_list.append(mea)

            sa_full_list = list()
            if args.sa_full:
                for number in args.sa_full:
                    sa_full = cl_handle.add_signal_align_predictions(number=int(number),
                                                                     add_basecall=True,
                                                                     two_d=args.two_d)
                    sa_full_list.append(sa_full)

            if args.tombo:
                for number in args.tombo:
                    tombo_full = cl_handle.add_tombo_labels(number=int(number))

            basecall_list = list()
            if args.basecall:
                events_list = []
                for number in args.basecall:
                    matches, mismatches = cl_handle.add_basecall_alignment_prediction(number=int(number),
                                                                                      trim=args.trim)
                    basecall = list()
                    basecall.extend(matches)
                    basecall.extend(mismatches)
                    basecall.sort(key=lambda x: x['raw_start'])
                    basecall_list.append(basecall)
                    events = cl_handle.get_basecalled_data_by_number(int(number))
                    events_list.append(events)
                if len(events_list) > 1:
                    main_plot = False
                    print("Plotting {}".format(args.f5_path))
                    plot_basecalled_sequences(events_list[0], events2=events_list[1],
                                              signal=cl_handle.aligned_signal.scaled_signal)

            if args.variant is not None:
                for number in args.variant:
                    cl_handle.add_variant_data(number=int(number), complement=False)
                    if args.two_d:
                        cl_handle.add_variant_data(number=int(number), complement=True)

            max_analysis_index = min(map(lambda x: len(x), [mea_list, sa_full_list])) if args.basecall else 0
            for i in range(max_analysis_index):
                print("Analyzing events for index {}".format(i))
                analyze_event_skips(mea_list[i], sa_full_list[i], basecall_list[i])

            kmer_info = None
            if args.kmers:
                aligned_segment = initialize_aligned_segment_wrapper(cl_handle.get_signalalign_events(sam=True))
                kmer_info = generate_kmer_info_for_plotting(args.kmers, aligned_segment)

            if main_plot:
                print(args.kmers)
                print("Plotting {}".format(f5_path))
                ps = PlotSignal(cl_handle.aligned_signal)
                ps.plot_alignment(save_fig_path=save_fig_path, plot_alpha=args.plot_alpha, kmer_info=kmer_info)
        except Exception as e:
            print("FAILED: {}: {}".format(e, f5_path))

    stop = timer()
    print("Running Time = {} seconds".format(stop - start), file=sys.stderr)


if __name__ == "__main__":
    main()
    raise SystemExit
