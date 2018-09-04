#!/usr/bin/env python
"""Plot speeds of maximum expected accuracy methods"""
########################################################################
# File: plot_mea_speeds.py
#  executable: plot_mea_speeds.py
#
# Author: Andrew Bailey
# History: Created 02/24/18
########################################################################

from __future__ import print_function
import sys
import os
import matplotlib
import colorsys

# matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from matplotlib.collections import LineCollection

import numpy as np

from collections import defaultdict
from timeit import default_timer as timer
from argparse import ArgumentParser
from signalalign.alignedsignal import AlignedSignal, CreateLabels
from py3helpers.utils import list_dir


def parse_args():
    parser = ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument('--dir', '-d', action='store',
                       dest='dir', default=None,
                       help="Path to directory of fast5 files to plot")

    group.add_argument('--f5_path', action='store', default=None,
                       dest='f5_path',
                       help="Path to a specific fast5 file to plot")

    parser.add_argument('--basecall', action='store_true', default=False,
                        dest='basecall', required=False,
                        help="Option to plot the most recent basecalled data aligned to reference")

    parser.add_argument('--mea', action='store_true', default=False,
                        dest='mea', required=False,
                        help="Option to plot the maximum expected accuracy alignment from signalalign")

    parser.add_argument('--sa_full', action='store_true', default=False,
                        dest='sa_full', required=False,
                        help="Option to plot all of the posterior probabilities from the signalalign output")

    parser.add_argument('--output_dir', action='store',
                        dest='output_dir', required=False, type=str,
                        help="If set, will write out to output directory otherwise the graph is just shown")

    parser.add_argument('--plot_alpha', action='store_true', default=False, dest="plot_alpha", required=False,
                        help="If set, will plot probability information as the alpha associated with each data point")

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
        # add raw starts
        if self.signal_h.raw_starts is not None:
            self.event_starts = self.signal_h.raw_starts

    def plot_alignment(self, save_fig_path=None, plot_alpha=False, plot_lines=False):
        """Plot the alignment between events and reference with the guide alignment and mea alignment
        :param plot_lines: boolean option, if set to true, lines will be plotted instead of points
        :param save_fig_path: if set, will write image to path
        :param plot_alpha: boolean option, if set to false will ignore probability information associated with data points
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
        colors = get_spaced_colors(len(self.names) + 1)
        color_selection = 1
        # plot signal alignments
        for i, alignment in enumerate(self.alignments):
            handle, = panel1.plot(alignment[0], alignment[1], color=colors[color_selection], alpha=0.8)
            color_selection += 1
            handles.append(handle)

        # plot predictions (signal align outputs)
        for i, prediction in enumerate(self.predictions):
            rgba_colors = np.tile(np.array(colors[color_selection]), (len(prediction[0]), 1))
            if plot_alpha or self.names[color_selection-1] == "full_signalalign":
                rgba_colors = np.insert(rgba_colors, 3, prediction[3].tolist(), axis=1)
            # rgba_colors[:, 3] = prediction[2].tolist()
            if plot_lines:
                lines = list(zip(list(zip(prediction[0].tolist(), [x+(i*0.01) for x in prediction[2].tolist()])),
                                 list(zip([x + y for x, y in zip(prediction[0].tolist(), prediction[1].tolist())],
                                          [x+(i*.01) for x in prediction[2].tolist()]))))
                line_segments = LineCollection(lines,
                                               linestyles='solid',
                                               colors=rgba_colors)
                handle = panel1.add_collection(line_segments)
            else:
                handle = panel1.scatter(prediction[0].tolist(), [x+(i*.01) for x in prediction[2].tolist()], marker='.',
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

        panel2.set_xlabel('Time')
        panel2.set_ylabel('Current (pA)')
        if self.event_starts is not None:
            panel2.set_xticks(self.event_starts, minor=True)

        panel2.axvline(linewidth=0.1, color="k")

        handle, = panel2.plot(self.signal_h.scaled_signal, color="black", lw=0.4)
        handles.append(handle)
        self.names.append("scaled_signal")

        box = panel1.get_position()
        panel1.set_position([box.x0, box.y0 + box.height * 0.1,
                             box.width, box.height * 0.9])

        # Put a legend below current axis
        panel1.legend(handles, self.names, loc='upper center', bbox_to_anchor=(0.5, -0.05),
                      fancybox=True, shadow=True, ncol=len(colors))

        # panel1.legend(handles, self.names, loc='upper right')
        if save_fig_path:
            plt.savefig(save_fig_path)
        else:
            plt.show()


def get_spaced_colors(n):
    """Create n evenly spaced RGB colors
    :param n: number of colors needed

    source: https://stackoverflow.com/questions/876853/generating-color-ranges-in-python
    """
    hsv_tuples = [(x*1.0/n, 0.5, 0.5) for x in range(n)]
    rgb_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), hsv_tuples)
    return list(rgb_tuples)


def main():
    """Plot event to reference labelled ONT nanopore reads """
    start = timer()

    args = parse_args()
    assert args.mea or args.sa_full or args.basecall, "--mea, --sa_full or --basecall must be set."

    if args.dir:
        for f5_path in list_dir(args.dir, ext="fast5"):
            save_fig_path = None
            cl_handle = CreateLabels(f5_path)
            if args.output_dir:
                save_fig_path = "{}.png".format(os.path.join(args.output_dir,
                                                             os.path.splitext(os.path.basename(f5_path))[0]))
            if args.mea:
                cl_handle.add_mea_labels()
            if args.sa_full:
                cl_handle.add_signal_align_predictions()
            if args.basecall:
                cl_handle.add_basecall_alignment()

            print("Plotting {}".format(args.f5_path))
            ps = PlotSignal(cl_handle.aligned_signal)
            ps.plot_alignment(save_fig_path=save_fig_path, plot_alpha=args.plot_alpha)

    else:
        save_fig_path = None
        cl_handle = CreateLabels(args.f5_path)
        if args.output_dir:
            save_fig_path = "{}.png".format(os.path.join(args.output_dir,
                                                         os.path.splitext(os.path.basename(args.f5_path))[0]))
        if args.mea:
            cl_handle.add_mea_labels()
        if args.sa_full:
            cl_handle.add_signal_align_predictions()
        if args.basecall:
            cl_handle.add_basecall_alignment_prediction()

        print("Plotting {}".format(args.f5_path))
        ps = PlotSignal(cl_handle.aligned_signal)
        ps.plot_alignment(save_fig_path=save_fig_path, plot_alpha=False, plot_lines=True)

    stop = timer()
    print("Running Time = {} seconds".format(stop - start), file=sys.stderr)


if __name__ == "__main__":
    main()
    raise SystemExit
