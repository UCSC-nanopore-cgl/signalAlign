#!/usr/bin/env python
"""Filter reads based on minimum average base quality and alignment"""
########################################################################
# File: filter_reads.py
#  executable: filter_reads.py
#
# Author: Andrew Bailey
# History: 10/15/18 Created
########################################################################

import os
import pysam
import shutil
import numpy as np
import pandas as pd
from contextlib import closing
from collections import defaultdict, namedtuple
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from signalalign.fast5 import Fast5
from signalalign.nanoporeRead import NanoporeRead
from signalalign.validateSignalAlignment import flag_large_gaps
from signalalign.visualization.plot_labelled_read import analyze_event_skips, CreateLabels
from signalalign.utils.sequenceTools import get_full_nucleotide_read_from_alignment
from py3helpers.utils import list_dir
from py3helpers.seq_tools import AlignmentSegmentWrapper


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--alignment_file', '-a', action='store',
                        dest='alignment_file', required=True, type=str, default=None,
                        help="Sam/Bam file with all alignment data")

    parser.add_argument('--fast5_dir', '-d', action='store', default=None,
                        dest='fast5_dir', required=True, type=str,
                        help="Directory of all fast5 files")

    parser.add_argument('--embedded', '-e', action='store_true', default=False,
                        dest='check_embedded', required=False,
                        help="Directory of all fast5 files")

    parser.add_argument('--quality_threshold', '-q', action='store', default=7,
                        dest='quality_threshold', required=False, type=float,
                        help="Minimum average base quality threshold. Default = 7")

    parser.add_argument('--output_dir', '-o', action='store',
                        dest='output_dir', required=True,
                        help="Directory to place plots and summary data")

    parser.add_argument('--max_reads', '-m', action='store', default=100,
                        dest='max_reads', required=False,
                        help="Maximum number of reads to sample")

    parser.add_argument('--from_pickle', '-p', action='store',
                        dest='from_pickle', required=False,
                        help="If we want to generate plots from a summary file")

    args = parser.parse_args()
    return args


def get_summary_info_table(read_ids):
    row_names = read_ids
    column_names = ["q_score_average", "seq_start_time",
                    "basecalled_accuracy", "no_mapping",
                    "chimera_mapping", "soft_clipped_percentage",
                    "num_secondary_mappings",
                    "num_flagged_gaps", "avg_flagged_gap_size", 'pass',
                    'other_errors']

    d = pd.DataFrame(0.0, columns=column_names, index=row_names)
    return d


def get_alignment_summary_info(fast5s, alignment_file, pass_threshold=7, gap_size=10, verbose=False, max_reads=100):
    """Filter fast5 files based on a quality threhsold and if there is an alignment"""
    # collect for every read
    fast5_dict = defaultdict()
    # loop through fast5s
    for fast5_path in fast5s:
        assert os.path.exists(fast5_path), "fast5 path does not exist: {}".format(fast5_path)
        f5h = NanoporeRead(fast5_path)
        f5h._initialize_metadata()
        read_name = f5h.read_label
        fast5_dict[read_name] = fast5_path
    print("Created read_id to fast5_path mapping")
    # summary data stored here
    mapped_reads = get_summary_info_table(list(fast5_dict.keys()))
    # grab aligned segment
    number = 0
    counter = 0
    reads_seen = []

    with closing(pysam.AlignmentFile(alignment_file, 'rb' if alignment_file.endswith("bam") else 'r')) as aln:
        for aligned_segment in aln.fetch():
            if counter > max_reads:
                break
            try:
                read_name = aligned_segment.qname.split("_")[0]
                fast5_path = fast5_dict[read_name]
                if read_name not in reads_seen:
                    reads_seen.append(read_name)
                    counter += 1
                print(fast5_path)
                cl_handle = CreateLabels(fast5_path, kmer_index=2)
                seq_start_time = cl_handle.raw_attributes['start_time']
                q_score_average = 0
                if aligned_segment.query_qualities is None:
                    print("Alignment done with fasta instead of fastq so read qualities will not be reported")
                else:
                    q_score_average = np.mean(aligned_segment.query_qualities)

                mapped_reads["q_score_average"][read_name] = q_score_average
                mapped_reads["seq_start_time"][read_name] = seq_start_time

                if aligned_segment.is_secondary or aligned_segment.is_unmapped \
                        or aligned_segment.is_supplementary or aligned_segment.has_tag("SA"):

                    if aligned_segment.is_secondary:
                        mapped_reads["num_secondary_mappings"][read_name] += 1
                    if aligned_segment.is_unmapped:
                        mapped_reads["no_mapping"][read_name] = 1
                    if aligned_segment.is_supplementary or aligned_segment.has_tag("SA"):
                        mapped_reads["chimera_mapping"][read_name] += 1
                else:
                    soft_clipped_percentage = \
                        1 - float(len(aligned_segment.query_alignment_sequence)) / len(aligned_segment.query_sequence)
                    mapped_reads["soft_clipped_percentage"][read_name] = soft_clipped_percentage

                    handle = AlignmentSegmentWrapper(aligned_segment)
                    handle.initialize()

                    accuracy = handle.alignment_accuracy()
                    mapped_reads["basecalled_accuracy"][read_name] = accuracy
                    try:
                        mea = cl_handle.add_mea_labels(number=int(number))
                        sa_full = cl_handle.add_signal_align_predictions(number=int(number), add_basecall=True)
                        all_basecall_data = []
                        for name, basecall_data in cl_handle.aligned_signal.prediction.items():
                            if "guide" in name:
                                all_basecall_data.extend(basecall_data)

                        alignment_summary = analyze_event_skips(mea, sa_full, all_basecall_data, generate_plot=False)
                        flagged_gaps_summary = flag_large_gaps(alignment_summary, gap_size, verbose=verbose)
                        counter = 0
                        total_distance = 0
                        for gap in flagged_gaps_summary:
                            if gap["mea_peak_distance"] > 10:
                                counter += 1
                                total_distance += gap["mea_peak_distance"]
                        if counter > 0:
                            mapped_reads["num_flagged_gaps"][read_name] = counter
                            mapped_reads["avg_flagged_gap_size"][read_name] = float(total_distance) / counter

                        if mapped_reads["q_score_average"][read_name] > pass_threshold:
                            mapped_reads["pass"][read_name] = 1

                    except KeyError:
                        mapped_reads["other_errors"][read_name] = 1
            except Exception as e:
                print(e)

        return mapped_reads


def print_summary_information(summary_pd):
    """Print out important summary stats"""
    passing_reads = summary_pd[summary_pd["pass"] == 1]
    if len(passing_reads) > 0:

        print("Fraction of passing reads: {}".format(len(passing_reads) / len(summary_pd)))

        pass_reads_w_gaps = passing_reads[passing_reads["num_flagged_gaps"] > 0]
        fraction_pass_flagged = len(pass_reads_w_gaps) / len(passing_reads)
        print("Fraction of passing reads which had at least one flagged gap: {}".format(fraction_pass_flagged))

        secondary_mappings_passed = passing_reads[passing_reads["num_secondary_mappings"] > 0]
        print("Fraction of passing reads which had secondary mappings: {}".format(len(secondary_mappings_passed)/len(passing_reads)))
    else:
        print("None of the reads passed")

    failed_reads = summary_pd[summary_pd["pass"] == 0]
    if len(failed_reads) > 0:

        non_mapped_reads = failed_reads[failed_reads["no_mapping"] > 0]
        print("Fraction of failed reads which did not map to reference: {}".format(len(non_mapped_reads)/len(failed_reads)))

        chimera_mapped_reads = failed_reads[failed_reads["chimera_mapping"] > 0]
        print("Fraction of failed reads which had supplemental or chimeric mappings: {}".format(len(chimera_mapped_reads)/len(failed_reads)))

        secondary_mappings = failed_reads[failed_reads["num_secondary_mappings"] > 0]
        print("Fraction of failed reads which had secondary mappings: {}".format(len(secondary_mappings)/len(failed_reads)))

        other_errors = failed_reads[failed_reads["other_errors"] > 0]
        print("Fraction of failed reads which had other errors: {}".format(len(other_errors)/len(failed_reads)))

    else:
        print("None of the reads failed!")


def plot_summary_information(summary_pd, output_dir=None):
    """Plot summary information regarding a sequencing run

     :param summary_pd: the pandas dataframe created by get_alignment_summary_info
     """
    if output_dir:
        output_file = os.path.join(output_dir, "accuracy_vs_qscore.jpg")
    else:
        output_file = None
    plot_accuracy_vs_qscore(summary_pd, save_fig_path=output_file)
    if output_dir:
        output_file = os.path.join(output_dir, "accuracy_vs_qscore.jpg")
    else:
        output_file = None
    plot_mapping_errors_vs_qscore(summary_pd, save_fig_path=output_file)

    return True


def plot_mapping_errors_vs_qscore(summary_pd, save_fig_path=None):
    """Plot how q-score correlates with mapping errors"""
    failed_reads = summary_pd[summary_pd["pass"] == 0]
    if len(failed_reads) > 0:
        q_scores = failed_reads["q_score_average"]
        unmapped_reads = q_scores[failed_reads["no_mapping"] != 0]
        chimera_reads = q_scores[failed_reads["chimera_mapping"] != 0]
        multiple_mapping = q_scores[failed_reads["num_secondary_mappings"] != 0]
        bins = np.linspace(0, max(q_scores), num=30)

        if len(unmapped_reads) > 0:
            plt.hist(unmapped_reads, bins=bins, color='blue', label="unmapped_reads")
        if len(chimera_reads) > 0:
            plt.hist(q_scores, bins=bins, color='red', label='chimeric_reads')
        if len(multiple_mapping) > 0:
            plt.scatter(q_scores, bins=bins, color='green', label='multiple_mapping_locations')

        plt.xlabel("Average Quality Score (phred)")
        plt.ylabel("Basecalling Accuracy: Matches/ Alignment Length")
        plt.legend()

        if save_fig_path:
            plt.savefig(save_fig_path)
        else:
            plt.show()

    return True


def plot_accuracy_vs_qscore(summary_pd, save_fig_path=None):
    """Plot how basecalling accuracy maps with quality scores"""
    accuracy = summary_pd['basecalled_accuracy']
    q_scores = summary_pd['q_score_average']

    plt.scatter(q_scores, accuracy, color='blue', label='accuracy vs avg. quality score')
    plt.xlabel("Average Quality Score (phred)")
    plt.ylabel("Basecalling Accuracy: Matches/ Alignment Length")
    plt.legend()

    if save_fig_path:
        plt.savefig(save_fig_path)
    else:
        plt.show()

    return True


def main():
    args = parse_args()

    fast5s = list_dir(args.fast5_dir, ext='fast5')
    assert len(fast5s) > 0, "Check fast5_dir. No files with fast5 extension found: {}".format(args.fast5_dir)
    assert os.path.isdir(args.output_dir), "{} is not a directory".format(args.output_dir)
    # Make output dir if it doesn't exist
    if args.from_pickle is not None:
        print("Loading sequencing summary info from pickle")
        summary_pd = pd.read_pickle(os.path.join(args.output_dir, "summary_info.pkl"))
    else:
        print("Creating alignment summary info")
        summary_pd = get_alignment_summary_info(fast5s, args.alignment_file,
                                                pass_threshold=7, max_reads=args.max_reads)
        summary_pd.to_pickle(os.path.join(args.output_dir, "summary_info.pkl"))

    print("Printing summary stats and creating plots")
    print_summary_information(summary_pd)
    plot_summary_information(summary_pd, output_dir=args.output_dir)


if __name__ == "__main__":
    main()
    raise SystemExit
