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
import sys
from contextlib import closing
from collections import defaultdict, namedtuple
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from timeit import default_timer as timer
from signalalign.fast5 import Fast5
from signalalign.filter_reads import parse_readdb
from signalalign.nanoporeRead import NanoporeRead
from signalalign.validateSignalAlignment import flag_large_gaps
from signalalign.visualization.plot_labelled_read import analyze_event_skips, CreateLabels
from signalalign.utils.sequenceTools import get_full_nucleotide_read_from_alignment
from signalalign.utils import multithread
from py3helpers.utils import list_dir
from py3helpers.seq_tools import AlignmentSegmentWrapper, sam_string_to_aligned_segment

COLUMN_NAMES = ["q_score_average", "seq_start_time",
                "basecalled_accuracy", "no_mapping",
                "chimera_mapping", "soft_clipped_percentage",
                "num_secondary_mappings",
                "num_flagged_gaps", "avg_flagged_gap_size", 'pass',
                'other_errors', 'seen', "map_q"]


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--alignment_file', '-a', action='store',
                        dest='alignment_file', required=False, type=str, default=None,
                        help="Sam/Bam file with all alignment data")

    parser.add_argument('--fast5_dir', '-d', action='store', default=None,
                        dest='fast5_dir', required=False, type=str,
                        help="Directory of all fast5 files")

    parser.add_argument('--quality_threshold', '-q', action='store', default=7,
                        dest='quality_threshold', required=False, type=float,
                        help="Minimum average base quality threshold. Default = 7")

    parser.add_argument('--readdb', '-b', action='store',
                        dest='readdb', required=False, type=str,
                        help="Path to readdb file for easy mapping")

    parser.add_argument('--output_dir', '-o', action='store',
                        dest='output_dir', required=True,
                        help="Directory to place plots and summary data")

    parser.add_argument('--max_reads', '-m', action='store', default=100,
                        dest='max_reads', required=False, type=int,
                        help="Maximum number of reads to sample")

    parser.add_argument('--number', '-n', action='store', default=0,
                        dest='number', required=False, type=int,
                        help="Which signalalignnment to grab from the fast5 file")

    parser.add_argument('--from_pickle', '-p', action='store',
                        dest='from_pickle', required=False,
                        help="If we want to generate plots from a summary file")

    parser.add_argument('--debug', action='store_true', default=False,
                        dest='debug', required=False,
                        help="Will stop multiprocess in order to debug underlying code")

    parser.add_argument('--workers', action="store", default=1,
                        dest='workers', required=False, type=int,
                        help="Number of workers to process reads")

    args = parser.parse_args()
    return args


def get_summary_info_table(read_ids):
    row_names = read_ids
    d = pd.DataFrame(0.0, columns=COLUMN_NAMES, index=row_names)
    return d


def get_summary_info_row(read_id):
    d = pd.DataFrame(0.0, columns=COLUMN_NAMES, index=[read_id])
    return d


def get_alignment_summary_info_withdb(alignment_file, readdb, read_dirs, pass_threshold=7,
                                      gap_size=10, verbose=False, max_reads=100, number=0):
    """Filter fast5 files based on a quality threhsold and if there is an alignment"""
    assert alignment_file.endswith("bam"), "Alignment file must be in BAM format: {}".format(alignment_file)
    # grab aligned segment
    seen_counter = 0

    with closing(pysam.AlignmentFile(alignment_file, 'rb')) as bamfile:
        name_indexed = pysam.IndexedReads(bamfile)
        print("Indexing bam file by read name.")
        name_indexed.build()
        print("Finished.")
        print("Looping through readdb file.")
        for name, fast5 in parse_readdb(readdb, read_dirs):
            try:
                iterator = name_indexed.find(name)
                # create ability to only grab x number of reads
                if seen_counter >= max_reads:
                    break
                # need to start the data table with first row
                if seen_counter == 0:
                    pd_data = get_summary_info_row(name)
                    big_table = pd_data
                else:
                    big_table.append(pd_data)
                    pd_data = get_summary_info_row(name)

                # start tracking data
                pd_data["seen"] = 1
                seen_counter += 1

                cl_handle = CreateLabels(fast5, kmer_index=2)
                seq_start_time = cl_handle.raw_attributes['start_time']
                pd_data["seq_start_time"] = seq_start_time

                for aligned_segment in iterator:
                    if aligned_segment.is_secondary or aligned_segment.is_unmapped \
                            or aligned_segment.is_supplementary or aligned_segment.has_tag("SA"):
                        if aligned_segment.is_secondary:
                            pd_data["num_secondary_mappings"] += 1
                        if aligned_segment.is_unmapped:
                            pd_data["no_mapping"] = 1
                        if aligned_segment.is_supplementary or aligned_segment.has_tag("SA"):
                            pd_data["chimera_mapping"] += 1
                    else:
                        pd_data["map_q"] = aligned_segment.mapq
                        soft_clipped_percentage = \
                            1 - float(len(aligned_segment.query_alignment_sequence)) / len(aligned_segment.query_sequence)
                        pd_data["soft_clipped_percentage"] = soft_clipped_percentage

                        handle = AlignmentSegmentWrapper(aligned_segment)
                        handle.initialize()

                        accuracy = handle.alignment_accuracy()
                        pd_data["basecalled_accuracy"] = accuracy
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
                                pd_data["num_flagged_gaps"] = counter
                                pd_data["avg_flagged_gap_size"] = float(total_distance) / counter

                            q_score_average = 0
                            if aligned_segment.query_qualities is None:
                                print("Alignment done with fasta instead of fastq so read qualities will not be reported")
                            else:
                                q_score_average = np.mean(aligned_segment.query_qualities)

                            pd_data["q_score_average"] = q_score_average
                            print("pd_data['q_score_average']", pd_data["q_score_average"][0])
                            if pd_data["q_score_average"][0] > pass_threshold:
                                pd_data["pass"] = 1

                        except Exception as e:
                            pd_data["other_errors"] = 1
                            print("ERROR {}: {}".format(fast5, e), file=sys.stderr)

            except KeyError:
                pd_data["other_errors"] = 1
                print("Found no alignments for {}".format(fast5))

        return big_table


def get_summary_info_service(work_queue, done_queue, service_name="summary_info"):
    # prep
    total_handled = 0
    failure_count = 0
    mem_usages = list()

    # catch overall exceptions
    try:
        for f in iter(work_queue.get, 'STOP'):
            # catch exceptions on each element
            try:
                pd_data = create_summary_pd(**f)
                done_queue.put(pd_data)
            except Exception as e:
                # get error and log it
                message = "{}:{}".format(type(e), str(e))
                error = "{} '{}' failed with: {}".format(service_name, multithread.current_process().name, message)
                print("[{}] ".format(service_name) + error)
                done_queue.put(error)
                failure_count += 1

            # increment total handling
            total_handled += 1

    except Exception as e:
        # get error and log it
        message = "{}:{}".format(type(e), str(e))
        error = "{} '{}' critically failed with: {}".format(service_name, multithread.current_process().name, message)
        print("[{}] ".format(service_name) + error)
        done_queue.put(error)

    finally:
        # logging and final reporting
        print("[%s] '%s' completed %d calls with %d failures"
              % (service_name, multithread.current_process().name, total_handled, failure_count))
        done_queue.put("{}:{}".format(multithread.TOTAL_KEY, total_handled))
        done_queue.put("{}:{}".format(multithread.FAILURE_KEY, failure_count))
        if len(mem_usages) > 0:
            done_queue.put("{}:{}".format(multithread.MEM_USAGE_KEY, ",".join(map(str, mem_usages))))


def multiprocess_get_summary_info(alignment_file, readdb, read_dirs, get_summary_args, worker_count=1, debug=False):
    """Multiprocess get summary info"""
    assert alignment_file.endswith("bam"), "Alignment file must be in BAM format: {}".format(alignment_file)
    # grab aligned segment
    data = pd.DataFrame([])

    with closing(pysam.AlignmentFile(alignment_file, 'rb')) as bamfile:
        name_indexed = pysam.IndexedReads(bamfile)
        print("Indexing bam file by read name.")
        name_indexed.build()
        print("Finished.")
        print("Looping through readdb file.")
        if debug:
            for name, fast5, iterator in parse_readdb_wrapper(parse_readdb(readdb, read_dirs), name_indexed):
                pd_line = create_summary_pd(iterator, fast5, name, **get_summary_args)
                data = data.append(pd_line)
        else:
            total, failure, messages, output = multithread.run_service2(
                get_summary_info_service, parse_readdb_wrapper(parse_readdb(readdb, read_dirs), name_indexed),
                get_summary_args, ['name', "fast5", "sam_lines"], worker_count)
            for pd_line in output:
                if isinstance(pd_line, pd.DataFrame):
                    data = data.append(pd_line)
    return data


def parse_readdb_wrapper(parse_readdb_generator, name_indexed):
    """Wrap the generator to return the name, fast5 and iterator"""
    for name, fast5 in parse_readdb_generator:
        try:
            iterator = name_indexed.find(name)
            sam_lines = []
            for aligned_segment in iterator:
                sam_lines.append(aligned_segment.to_string())
            yield name, fast5, sam_lines
        except KeyError:
            print("Found no alignments for {}".format(fast5))


def create_summary_pd(sam_lines, fast5, name, number=0, pass_threshold=7, gap_size=10, verbose=False):
    """Create summary pd from an iterator"""
    pd_data = get_summary_info_row(name)

    cl_handle = CreateLabels(fast5, kmer_index=2)
    seq_start_time = cl_handle.raw_attributes['start_time']
    pd_data["seq_start_time"] = seq_start_time

    for sam_string in sam_lines:
        aligned_segment = sam_string_to_aligned_segment(sam_string)
        if aligned_segment.is_secondary or aligned_segment.is_unmapped \
                or aligned_segment.is_supplementary or aligned_segment.has_tag("SA"):
            if aligned_segment.is_secondary:
                pd_data["num_secondary_mappings"] += 1
            if aligned_segment.is_unmapped:
                pd_data["no_mapping"] = 1
            if aligned_segment.is_supplementary or aligned_segment.has_tag("SA"):
                pd_data["chimera_mapping"] += 1
        else:
            pd_data["map_q"] = aligned_segment.mapq
            soft_clipped_percentage = \
                1 - float(len(aligned_segment.query_alignment_sequence)) / len(aligned_segment.query_sequence)
            pd_data["soft_clipped_percentage"] = soft_clipped_percentage

            handle = AlignmentSegmentWrapper(aligned_segment)
            handle.initialize()

            accuracy = handle.alignment_accuracy()
            pd_data["basecalled_accuracy"] = accuracy
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
                    pd_data["num_flagged_gaps"] = counter
                    pd_data["avg_flagged_gap_size"] = float(total_distance) / counter

                q_score_average = 0
                if aligned_segment.query_qualities is None:
                    print("Alignment done with fasta instead of fastq so read qualities will not be reported")
                else:
                    q_score_average = np.mean(aligned_segment.query_qualities)

                pd_data["q_score_average"] = q_score_average
                print("pd_data['q_score_average']", pd_data["q_score_average"][0])
                if pd_data["q_score_average"][0] > pass_threshold:
                    pd_data["pass"] = 1

            except Exception as e:
                pd_data["other_errors"] = 1
                print("ERROR {}: {}".format(fast5, e), file=sys.stderr)
    return pd_data

def get_alignment_summary_info(fast5s, alignment_file, pass_threshold=7, gap_size=10, verbose=False,
                               max_reads=100, number=0):
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
    seen_counter = 0
    reads_seen = set()
    print("first_len reads_seen: {}".format(len(reads_seen)), file=sys.stderr)

    with closing(pysam.AlignmentFile(alignment_file, 'rb' if alignment_file.endswith("bam") else 'r')) as aln:
        for aligned_segment in aln.fetch(until_eof=True):
            if seen_counter > max_reads:
                break
            try:
                print("reads_seen: {}".format(len(reads_seen)), file=sys.stderr)
                read_name = aligned_segment.qname.split("_")[0]
                fast5_path = fast5_dict[read_name]
                if read_name not in reads_seen:
                    reads_seen |= {read_name}
                    seen_counter += 1
                    mapped_reads["seen"][read_name] = 1
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
                        or aligned_segment.is_supplementary or aligned_segment.has_tag("SA") \
                        or q_score_average < pass_threshold:

                    if aligned_segment.is_secondary:
                        mapped_reads["num_secondary_mappings"][read_name] += 1
                    if aligned_segment.is_unmapped:
                        mapped_reads["no_mapping"][read_name] = 1
                    if aligned_segment.is_supplementary or aligned_segment.has_tag("SA"):
                        mapped_reads["chimera_mapping"][read_name] += 1
                else:
                    mapped_reads["map_q"][read_name] = aligned_segment.mapq

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
                print(e, file=sys.stderr)

        return mapped_reads[mapped_reads["seen"] == 1]


def print_summary_information(summary_pd, pass_threshold=7):
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

        other_errors = failed_reads[failed_reads["q_score_average"] < pass_threshold]
        print("Fraction of failed reads which had q_score_average < {}: {}".format(pass_threshold, len(other_errors)/len(failed_reads)))

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
        output_file = os.path.join(output_dir, "mapping_errors_vs_qscore.jpg")
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

        plt.figure(figsize=(6, 8))
        panel1 = plt.axes([0.1, 0.1, .9, .9])

        if len(unmapped_reads) > 0:
            panel1.hist(unmapped_reads, bins=bins, color='blue', label="unmapped_reads")
        if len(chimera_reads) > 0:
            panel1.hist(q_scores, bins=bins, color='red', label='chimeric_reads')
        if len(multiple_mapping) > 0:
            panel1.scatter(q_scores, bins=bins, color='green', label='multiple_mapping_locations')

        panel1.set_xlabel("Average Quality Score (phred)")
        panel1.set_ylabel("Basecalling Accuracy: Matches/ Alignment Length")
        panel1.legend()

        if save_fig_path:
            plt.savefig(save_fig_path)
        else:
            plt.show()

    return True


def plot_accuracy_vs_qscore(summary_pd, save_fig_path=None):
    """Plot how basecalling accuracy maps with quality scores"""
    accuracy = summary_pd['basecalled_accuracy']
    q_scores = summary_pd['q_score_average']

    plt.figure(figsize=(6, 8))
    panel1 = plt.axes([0.1, 0.1, .9, .9])

    panel1.scatter(q_scores, accuracy, color='blue', label='accuracy vs avg. quality score')
    panel1.set_xlabel("Average Quality Score (phred)")
    panel1.set_ylabel("Basecalling Accuracy: Matches/ Alignment Length")
    panel1.legend()

    if save_fig_path:
        plt.savefig(save_fig_path)
    else:
        plt.show()

    return True


def main():
    start = timer()

    args = parse_args()
    assert os.path.isdir(args.output_dir), "{} is not a directory".format(args.output_dir)

    if args.from_pickle is not None:
        print("Loading sequencing summary info from pickle")
        summary_pd = pd.read_pickle(args.from_pickle)
    else:
        assert args.fast5_dir is not None, "Must select fast5_dir if not loading from pickle file"
        fast5s = list_dir(args.fast5_dir, ext='fast5')
        assert len(fast5s) > 0, "Check fast5_dir. No files with fast5 extension found: {}".format(args.fast5_dir)
        print("Creating alignment summary info")

        get_summary_args = dict(number=args.number, pass_threshold=args.quality_threshold,
                                gap_size=10, verbose=False)
        summary_pd = multiprocess_get_summary_info(args.alignment_file, args.readdb, [args.fast5_dir],
                                                   get_summary_args, worker_count=args.workers, debug=args.debug)

        summary_pd.to_pickle(os.path.join(args.output_dir, "summary_info.pkl"))

    print("Printing summary stats and creating plots")
    print_summary_information(summary_pd, pass_threshold=args.quality_threshold)
    plot_summary_information(summary_pd, output_dir=args.output_dir)
    stop = timer()
    print("Running Time = {} seconds".format(stop - start), file=sys.stderr)


if __name__ == "__main__":
    main()
    raise SystemExit
