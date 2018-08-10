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
from collections import defaultdict
from timeit import default_timer as timer
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
from signalalign.alignedsignal import AlignedSignal, CreateLabels


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
        self.guide_alignments = defaultdict(list)

    def get_alignments(self):
        """Format alignments from AlignedSignal to make it possible to plot easily"""
        for name, label in self.signal_h.label.items():
            self.names.append(name)
            self.alignments.append([label['raw_start'], label['reference_index']])

        # predictions can have multiple alignments for each event and are not a path
        for name, prediction in self.signal_h.prediction.items():
            self.names.append(name)
            self.predictions.append([prediction['raw_start'], prediction['reference_index'],
                                     prediction['posterior_probability']])
        # TODO make it possible to add multiple guide alignments
        for whole_guide_name, guide in self.signal_h.guide.items():
            for _, sub_guide in guide.items():
                self.guide_alignments[whole_guide_name].append([sub_guide['raw_start'], sub_guide['reference_index']])
            # gather tail ends of alignments
            self.names.append(whole_guide_name)

    def plot_alignment(self):
        """Plot the alignment between events and reference with the guide alignment and mea alignment

        signal: normalized signal
        posterior_matrix: matrix with col = reference , rows=events
        guide_cigar: guide alignment before signal align
        events: doing stuff
        mea_alignment: final alignment between events and reference
        """

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
        colors = get_spaced_colors(len(self.names)+1)
        color_selection = 1
        # plot signal alignments
        for i, alignment in enumerate(self.alignments):
            handle, = panel1.plot(alignment[0], alignment[1], color=colors[color_selection][:3], alpha=0.8)
            color_selection += 1
            handles.append(handle)

        # plot predictions (signal align outputs)
        for i, prediction in enumerate(self.predictions):
            rgba_colors = np.zeros((len(prediction[0]), 4))
            rgba_colors[:, 3] = prediction[2].tolist()
            handle = panel1.scatter(prediction[0].tolist(), prediction[1].tolist(), marker='.', c=rgba_colors)
            handles.append(handle)

        # plot original basecalled events
        if self.guide_alignments:
            for i, (name, guide) in enumerate(self.guide_alignments.items()):
                for j, sub_guide in enumerate(guide):
                    handle, = panel1.plot(sub_guide[0], sub_guide[1], c=colors[color_selection])
                handles.append(handle)
                color_selection += 1

        panel2.set_xlabel('Time')
        panel2.set_ylabel('Current (pA)')
        if len(self.predictions) > 0:
            panel2.set_xticks(self.predictions[0][0], minor=True)
        elif len(self.alignments) > 0:
            panel2.set_xticks(self.alignments[0][0], minor=True)
        else:
            panel2.set_xticks(self.signal_h.scaled_signal, minor=True)

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

        plt.show()


def get_spaced_colors(n):
    max_value = 16581375 #255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]
    return [(int(i[:2], 16)/255, int(i[2:4], 16)/255, int(i[4:], 16)/255, 1) for i in colors]


def main():
    """Main docstring"""
    start = timer()
    # sam = "/Users/andrewbailey/CLionProjects/nanopore-RNN/signalAlign/bin/`output/tempFiles_alignment/tempFiles_miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_strand/temp_sam_file_5048dffc-a463-4d84-bd3b-90ca183f488a.sam"\

    rna_read = "/Users/andrewbailey/CLionProjects/nanopore-RNN/test_files/minion-reads/rna_reads/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_36_ch_218_strand.fast5"
    # dna_read = "/Users/andrewbailey/CLionProjects/nanopore-RNN/test_files/minion-reads/canonical/miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch100_read280_strand.fast5"
    dna_read = "/Users/andrewbailey/CLionProjects/nanopore-RNN/nanotensor/tests/test_files/minion-reads/canonical/miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch100_read280_strand.fast5"
    dna_read2 = "/Users/andrewbailey/CLionProjects/nanopore-RNN/test_files/minion-reads/canonical/miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_strand.fast5"
    # dna_read3 = "/Users/andrewbailey/CLionProjects/nanopore-RNN/test_files/minion-reads/canonical/over_run/miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_strand.fast5"

    reference = "/Users/andrewbailey/CLionProjects/nanopore-RNN/test_files/reference-sequences/ecoli_k12_mg1655.fa"

    read_dir = "/Users/andrewbailey/data/directRNA_withEvents/0"
    for read in os.listdir(read_dir):
        path = os.path.join(read_dir, read)
        print(path)
        cl_handle = CreateLabels(path)
        cl_handle.add_guide_alignment()
        cl_handle.add_mea_labels()
        cl_handle.add_basecall_event_table_alignment("Analyses/Basecall_1D_002/BaseCalled_template/Events")

        ps = PlotSignal(cl_handle.aligned_signal)
        print("Plotting {}".format(path))
        ps.plot_alignment()


    # rna_read = "/Users/andrewbailey/data/directRNA_withEvents/0/DEAMERNANOPORE_20170922_FAH26525_MN16450_mux_scan_MA_821_R94_NA12878_mRNA_09_22_17_34495_read_114_ch_43_strand.fast5"
    # cl_handle = CreateLabels(rna_read)
    # cl_handle.add_guide_alignment()
    # cl_handle.add_mea_labels()
    # cl_handle.add_basecall_event_table_alignment("Analyses/Basecall_1D_002/BaseCalled_template/Events")
    #
    # ps = PlotSignal(cl_handle.aligned_signal)
    # print("Plotting {}".format(rna_read))
    # ps.plot_alignment()

    # test = CreateLabels(rna_read)
    # test.add_guide_alignment()
    # test.add_mea_labels()
    # test.add_signal_align_predictions()
    # test.add_nanoraw_labels(reference)
    # test.add_eventalign_labels()
    # ps = PlotSignal(test.aligned_signal)
    # print("Plotting")
    # ps.plot_alignment()


    stop = timer()
    print("Running Time = {} seconds".format(stop - start), file=sys.stderr)


if __name__ == "__main__":
    main()
    raise SystemExit
