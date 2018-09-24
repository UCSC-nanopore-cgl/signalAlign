#!/usr/bin/env python
"""Test if there is a major discrepancy between load from raw and original basecalled event table"""
########################################################################
# File: verify_load_from_raw.py
#  executable: verify_load_from_raw.py
#
# Author: Andrew Bailey
# History: 9/14/17 Created
########################################################################

from __future__ import print_function
import sys
import os
import colorsys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from matplotlib.collections import LineCollection
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

    parser.add_argument('--basecall', nargs='+',
                        dest='basecall', required=True,
                        help="Must select the two basecalled sections to compare")

    parser.add_argument('--output_dir', action='store',
                        dest='output_dir', required=False, type=str,
                        help="If set, will write out to output directory otherwise the graph is just shown")

    args = parser.parse_args()
    return args


def get_raw_start_delta(original_matches, new_matches):
    """Compare the distance between two sets of matched alignment information
    :param original_matches: matches from original basecalled table
    :param new_matches: matches from load from raw or new basecalled table
    :return: np array of raw_start deltas
    """
    assert (original_matches["reference_index"] == new_matches["reference_index"]).all, \
        "Reference index is different between original " \
        "matches and new matches. {} != {}".format(original_matches["reference_index"], new_matches["reference_index"])
    deltas = original_matches["raw_start"] - new_matches["raw_start"]

    return deltas


def plot_deltas_vs_raw_start(deltas, raw_starts, save_fig_path=None):
    """Simple scatter plot of deltas compared to the raw start
    :param deltas: array to plot on y axis
    :param raw_starts: array to match with deltas for x axis placement
    :param save_fig_path: if set, will save figure to location
    """
    plt.figure(figsize=(6, 8))
    panel1 = plt.axes([0.1, 0.1, 0.8, 0.8])
    panel1.set_title("Delta vs Raw starts")
    panel1.set_xlabel('Raw Start')
    panel1.set_ylabel('Delta (Original - load_from_raw)')
    panel1.grid(color='black', linestyle='-', linewidth=1)
    panel1.scatter(raw_starts, deltas, color='red', label='Delta (og-new)', alpha=.5, s=10)
    panel1.legend()
    if save_fig_path:
        plt.savefig(save_fig_path)
    else:
        plt.show()
        pass


def plot_deltas_as_distribution(deltas, save_fig_path=None, hist=False):
    """Simple density plot of deltas

    :param deltas: input array to plot
    :param hist: boolean option to make histogram instead of density
    :param save_fig_path: if set, will save figure to path instead of showing
    """
    plt.figure(figsize=(6, 8))
    panel1 = plt.axes([0.1, 0.1, 0.8, 0.8])
    panel1.set_xlabel('Delta (Original - load_from_raw)')
    panel1.set_ylabel('Density')
    panel1.grid(color='black', linestyle='-', linewidth=1)
    if hist:
        panel1.set_title("Delta Histogram")
        panel1.hist(deltas, bins=len(set(deltas)), color='blue', edgecolor='black')
    else:
        panel1.set_title("Delta Density")
        sns.distplot(deltas, hist=False, kde=True,
                     color='darkblue',
                     hist_kws={'edgecolor': 'black'},
                     kde_kws={'linewidth': 3})

    panel1.legend()
    if save_fig_path:
        plt.savefig(save_fig_path)
    else:
        plt.show()
        pass


def main(args=None):
    """Go through fast5 files and if there are two basecalled tables, compare and flag if different"""
    start = timer()
    # get args
    args = args if args is not None else parse_args()
    # get fast5s
    if args.dir:
        f5_locations = list_dir(args.dir, ext="fast5")
    else:
        f5_locations = [args.f5_path]

    assert len(args.basecall) == 2, "Must select two basecalled sections for comparison. " \
                                    "You selected {} basecalled tables".format(len(args.basecall))
    # loop through fast5s
    for f5_path in f5_locations:
        try:
            # grab basecalled data
            cl_handle = CreateLabels(f5_path)
            og_matches, og_mismatches = cl_handle.add_basecall_alignment_prediction(number=int(args.basecall[0]))
            new_matches, new_mismatches = cl_handle.add_basecall_alignment_prediction(number=int(args.basecall[1]))
            # get and plot deltas
            deltas = get_raw_start_delta(og_matches, new_matches)
            mm_deltas = get_raw_start_delta(og_mismatches, new_mismatches)
            all_deltas = np.append(deltas, mm_deltas)
            # if max(all_deltas) != 0:
                # plot_deltas_vs_raw_start(all_deltas, np.append(og_matches["raw_start"], og_mismatches['raw_start']))
                # plot_deltas_as_distribution(all_deltas, hist=True)
            percent_nonzero = np.count_nonzero(all_deltas)/len(all_deltas)
            if percent_nonzero > 0.3:
                print("Check read {}".format(f5_path))
        except KeyError as err:
            print("{}: {}".format(err, f5_path))
            continue

    stop = timer()
    print("Running Time = {} seconds".format(stop - start), file=sys.stderr)


if __name__ == "__main__":
    main()
    raise SystemExit
