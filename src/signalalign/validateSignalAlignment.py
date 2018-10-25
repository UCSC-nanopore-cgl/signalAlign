#!/usr/bin/env python3
"""Run signal-to-reference alignments
"""
from __future__ import print_function

import glob
import shutil
from signalalign.utils.fileHandlers import FolderHandler
from signalalign.motif import getDegenerateEnum
import sys
import math
import tempfile

from timeit import default_timer as timer
from signalalign.alignedsignal import CreateLabels
from signalalign.signalAlignment import multithread_signal_alignment, create_signalAlignment_args
from signalalign.visualization.plot_labelled_read import *


EVENTS="events"
EVENT_COUNT="event_count"
PEAK_DISTANCE="peak_distance"
MEA_PEAK_DISTANCE="mea_peak_distance"
CENTER_EVENT_ID= "ff_event_id"
CENTER_EVENT_START="ff_event_raw_start"

def parse_args(args=None):
    parser = ArgumentParser(description=__doc__)

    # input files
    parser.add_argument('--file_directory', '-i', action='store',
                        dest='files_dir', required=False, type=str, default=None,
                        help="directory with fast5 reads to align")
    parser.add_argument('--fast5_glob', '-g', action='store',
                        dest='fast5_glob', required=False, type=str, default=None,
                        help="glob matching fast5 reads to align")

    # output files
    parser.add_argument('--output_location', '-o', action='store', dest='out',
                        required=False, type=str, default=None,
                        help="directory to put the alignments. tempdir will be used if not specified")
    # parser.add_argument('--save_plots', '-p', action='store_true', dest='plots',
    #                     default=None, help='save plots in output directory')
    # parser.add_argument('--no_plots', '-P', action='store_false', dest='plots',
    #                     help='do not generate plots')

    # validation arguments
    parser.add_argument('--alignment_distance_threshold', '-d', action='store', dest='aln_dist_threshold',
                        required=False, type=int, default=10,
                        help="alignment distance threshold after which aligned events will be flagged")
    parser.add_argument('--force_kmer_event_alignment', '-k', action='store_true', dest="perform_kmer_event_alignment",
                        default=None, help="Forces SignalAlign to infer kmer-to-event alignment "
                                           "(by default, it will perform only if missing event table). ")
    parser.add_argument('--prevent_kmer_event_alignment', '-K', action='store_false', dest="perform_kmer_event_alignment",
                        help="Prevents SignalAlign from infer kmer-to-event alignment "
                             "(event if the fast5 is missing the event table). ")

    # signal align arguments
    parser.add_argument('--ref', '-r', action='store',
                        dest='ref', required=True, type=str,
                        help="reference sequence to align to, in FASTA")
    parser.add_argument('--in_template_hmm', '-T', action='store', dest='in_T_Hmm',
                        required=False, type=str, default=None,
                        help="input HMM for template events, if you don't want the default")
    parser.add_argument('--in_complement_hmm', '-C', action='store', dest='in_C_Hmm',
                        required=False, type=str, default=None,
                        help="input HMM for complement events, if you don't want the default")
    parser.add_argument('--templateHDP', '-tH', action='store', dest='templateHDP', default=None,
                        help="template serialized HDP file")
    parser.add_argument('--complementHDP', '-cH', action='store', dest='complementHDP', default=None,
                        help="complement serialized HDP file")
    parser.add_argument('--degenerate', '-x', action='store', dest='degenerate', default="variant",
                        help="Specify degenerate nucleotide options: "
                             "variant -> {ACGT}, cytosine2 -> {CE} cytosine3 -> {CEO} adenosine -> {AI}")
    parser.add_argument('--stateMachineType', '-smt', action='store', dest='stateMachineType', type=str,
                        default="threeState", help="decide which model to use, threeState by default")
    parser.add_argument('--threshold', '-t', action='store', dest='threshold', type=float, required=False,
                        default=None, help="posterior match probability threshold, Default: 0.01")
    parser.add_argument('--diagonalExpansion', '-e', action='store', dest='diag_expansion', type=int,
                        required=False, default=None, help="number of diagonals to expand around each anchor")
    parser.add_argument('--constraintTrim', '-m', action='store', dest='constraint_trim', type=int,
                        required=False, default=None, help='amount to remove from an anchor constraint')
    parser.add_argument("--kmer_size", action='store', dest="kmer_size", default=5, required=False,
                        help="size of kmers in fast5 file")

    # misc execution parameters
    parser.add_argument('--verbose', '-v', action='store_true', dest='verbose', required=False,
                        default=False, help="prints extra information")
    parser.add_argument('--dont_run_sa', action='store_true', dest='dont_run_sa', required=False,
                        default=False, help="Does not run signal align")

    args = parser.parse_args(args)

    # either: a) exactly one of [files_dir, fast5_glob], b) validation
    if not ((args.files_dir is None) ^ (args.fast5_glob is None)):
        raise Exception("Exactly one of --file_directory and --fast5_glob must be set")

    return args


def get_all_event_summaries(fast5s, alignment_args, aln_dist_threshold=10, generate_plot=True, verbose=False,
                            run_sa=True):
    start = timer()
    if run_sa:
        all_event_summaries = dict()
        # argments required for this to work
        alignment_args['output_format'] = 'full'
        alignment_args['embed'] = True
        alignment_args['check_for_temp_file_existance'] = False

        # run signal align
        print("\n[validateSignalAlignment]: running signalAlign in preparation for validation", file=sys.stdout)
        multithread_signal_alignment(alignment_args, fast5s, debug=True)

    print("\n[validateSignalAlignment]: performing validation", file=sys.stdout)
    for f5_path in fast5s:
        print("[validateSignalAlignment]: validating {}".format(f5_path), file=sys.stdout)
        # get data
        cl_handle = CreateLabels(f5_path)
        mea = cl_handle.add_mea_labels()
        sa_full = cl_handle.add_signal_align_predictions()
        matches, mismatches = cl_handle.add_basecall_alignment_prediction()
        basecall = list()
        basecall.extend(matches)
        basecall.extend(mismatches)
        basecall.sort(key=lambda x: x['raw_start'])

        event_summaries = analyze_event_skips(mea, sa_full, basecall, generate_plot=generate_plot)
        all_event_summaries[f5_path] = event_summaries

        # print stats on all event summaries
        total_events = len(event_summaries)
        bucket_size = 5.0
        max_event_aln_bucket = int(math.floor(max(list(map(lambda x: x[ABS_SA_ALIGNMENT_DIFF], event_summaries))) / bucket_size))
        aln_distance_buckets = {x:0 for x in range(0, max_event_aln_bucket + 1)}
        for summary in event_summaries:
            aln_distance_buckets[abs(math.floor(summary[ABS_SA_ALIGNMENT_DIFF] / bucket_size))] += 1
        max_bucket_count = math.log2(max(aln_distance_buckets.values()))
        print("Aligned Events: Absolute distance counts between SA aligned position and guide alignment position (log2)")
        running_total = 1.0 * total_events
        for x in range(0, max_event_aln_bucket + 1):
            hash_count = int(32.0 * math.log2(aln_distance_buckets[x]) / max_bucket_count) if aln_distance_buckets[x] != 0 else 0
            print("\t{:3d} to {:3d}: {} {:2.3f} {} \t raw_count:{:5d}  percentile: {:.5f}".format(
                int(x * bucket_size), int((x + 1) * bucket_size - 1), "#" * hash_count, math.log2(aln_distance_buckets[x]),
                " " * (32 - hash_count), aln_distance_buckets[x], running_total / total_events))
            running_total -= aln_distance_buckets[x]

        # gather stats on consecutive events
        if verbose: print("Consecutive events flagged by distance threshold")
        total_failed_events = 0
        all_flagged_event_sets = list()

        current_flagged_events = None
        for summary in event_summaries:
            # starting or in a flagged section
            if summary[ABS_SA_ALIGNMENT_DIFF] > aln_dist_threshold:
                if current_flagged_events is None:
                    current_flagged_events = list()
                current_flagged_events.append(summary)
                total_failed_events += 1

            # finishing a flagged section
            elif current_flagged_events is not None:
                current_flagged_mea_events = list(filter(lambda x: x[MEA], current_flagged_events))
                flagged_event_summary = {
                    EVENTS: current_flagged_events,
                    EVENT_COUNT: len(current_flagged_events),
                    PEAK_DISTANCE: max(list(map(lambda x: x[ABS_SA_ALIGNMENT_DIFF], current_flagged_events))),
                    MEA_PEAK_DISTANCE: 0 if len(current_flagged_mea_events) == 0 else max(list(map(
                        lambda x: x[ABS_SA_ALIGNMENT_DIFF], current_flagged_mea_events))),
                    CENTER_EVENT_ID: int(np.mean(list(map(lambda x: x[EVENT_INDEX], current_flagged_events)))),
                    CENTER_EVENT_START: int(np.mean(list(map(lambda x: x[RAW_START], current_flagged_events)))),
                }
                all_flagged_event_sets.append(flagged_event_summary)
                current_flagged_events = None
                if verbose:
                    print("\tFlaggedEventSet:\tcount:{:3d}\tpeak_dist:{:3d}\tpeak_mea_dist:{:3d}\tcenter_event:{:5d}\tcenter_event_start:{:7d}\t".format(
                        flagged_event_summary[EVENT_COUNT], flagged_event_summary[PEAK_DISTANCE],
                        flagged_event_summary[MEA_PEAK_DISTANCE], flagged_event_summary[CENTER_EVENT_ID],
                        flagged_event_summary[CENTER_EVENT_START]
                    ))

            # not in or ending a flagged section
            else:
                pass

        # summarize all flaged events
        print("Found {} flagged event sets".format(len(all_flagged_event_sets)))
        print("Of {} total events, {} were flagged ({:2.5f}%)".format(
            total_events, total_failed_events, 100.0 * total_failed_events / total_events))
        if len(all_flagged_event_sets) > 0:
            # maybe do some analysis of all these?
            pass

    stop = timer()
    print("Running Time = {} seconds".format(stop - start), file=sys.stderr)
    return all_event_summaries


def main():

    args = parse_args()

    # general sanity checks
    if not os.path.isfile(args.ref):
        print("[validateSignalAlignment] Did not find valid reference file", file=sys.stderr)
        sys.exit(1)

    # make directory to put temporary files and output location
    output_location = os.path.abspath(args.out)
    if output_location is None:
        output_root = tempfile.TemporaryDirectory()
    else:
        if not os.path.isdir(output_location):
            os.mkdir(output_location)
        output_root = output_location
    temp_signal_align_dir = os.path.join(output_root, "temp_signalAlign")
    if os.path.isdir(temp_signal_align_dir):
        shutil.rmtree(temp_signal_align_dir)
        assert not os.path.isdir(temp_signal_align_dir)

    # get input files
    input_glob = args.fast5_glob if args.fast5_glob is not None else os.path.join(args.files_dir, "*.fast5")
    orig_fast5s = glob.glob(input_glob)
    if len(orig_fast5s) == 0:
        print("[validateSignalAlignment] Did not find any files matching {}".format(input_glob), file=sys.stderr)
        sys.exit(1)
    print("[validateSignalAlignment] Found {} files".format(len(orig_fast5s)), file=sys.stdout)
    #
    alignment_args = create_signalAlignment_args(
        # signal align args
        destination=temp_signal_align,
        stateMachineType=args.stateMachineType,
        bwa_reference=args.ref,
        forward_reference=args.ref,
        in_templateHmm=args.in_T_Hmm,
        in_complementHmm=args.in_C_Hmm,
        in_templateHdp=args.templateHDP,
        in_complementHdp=args.complementHDP,
        threshold=args.threshold,
        diagonal_expansion=args.diag_expansion,
        constraint_trim=args.constraint_trim,
        degenerate=getDegenerateEnum(args.degenerate),
        perform_kmer_event_alignment=args.perform_kmer_event_alignment,
        alignment_file=args.alignment_file
    )
    get_all_event_summaries(orig_fast5s, alignment_args, aln_dist_threshold=args.aln_dist_threshold,
                            verbose=args.verbose, run_sa=not args.dont_run_sa)


if __name__ == "__main__":
    sys.exit(main())
