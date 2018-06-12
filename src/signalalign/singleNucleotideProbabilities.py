#!/usr/bin/env python3
"""Run signal-to-reference alignments
"""
from __future__ import print_function

import numpy as np
import glob
import shutil
import subprocess
import os
import sys
import pysam
from argparse import ArgumentParser
from random import shuffle
from contextlib import closing
from signalalign.nanoporeRead import NanoporeRead
from signalalign.signalAlignment import multithread_signal_alignment
from signalalign.scripts.alignmentAnalysisLib import CallMethylation
from signalalign.utils.fileHandlers import FolderHandler
from signalalign.utils.sequenceTools import get_full_nucleotide_read_from_alignment, replace_periodic_reference_positions
from signalalign.utils.multithread import *
from signalalign.motif import getDegenerateEnum
from signalalign.event_detection import generate_events_and_alignment


READ_NAME_KEY = "read_name"
FAST5_NAME_KEY = "fast5_name"
LENGTH_KEY = "length"
CONSENSUS_IDENTITY_KEY = "consensus_identity"
POSTERIOR_IDENTITY_KEY = "posterior_identity"
POSITIVE_STRAND_KEY = "positive_strand"
ALIGNED_IDENTITY_KEY="aligned_identity"
NO_GAP_ALIGNED_IDENTITY_KEY="no_gap_aligned_identity"
MEM_USAGES="mem_usages"

VALID_IDENTITY_RATIO = .85

def parse_args(args=None):
    parser = ArgumentParser(description=__doc__)

    parser.add_argument('--file_directory', '-d', action='store',
                        dest='files_dir', required=False, type=str, default=None,
                        help="directory with fast5 reads to align")
    parser.add_argument('--fast5_glob', '-g', action='store',
                        dest='fast5_glob', required=False, type=str, default=None,
                        help="glob matching fast5 reads to align")
    parser.add_argument('--ref', '-r', action='store',
                        dest='ref', required=True, type=str,
                        help="reference sequence to align to, in FASTA")
    parser.add_argument('--output_location', '-o', action='store', dest='out',
                        required=True, type=str, default=None,
                        help="directory to put the alignments")
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
    parser.add_argument('--stateMachineType', '-smt', action='store', dest='stateMachineType', type=str,
                        default="threeState", help="decide which model to use, threeState by default")
    parser.add_argument('--threshold', '-t', action='store', dest='threshold', type=float, required=False,
                        default=None, help="posterior match probability threshold, Default: 0.01")
    parser.add_argument('--diagonalExpansion', '-e', action='store', dest='diag_expansion', type=int,
                        required=False, default=None, help="number of diagonals to expand around each anchor")
    parser.add_argument('--constraintTrim', '-m', action='store', dest='constraint_trim', type=int,
                        required=False, default=None, help='amount to remove from an anchor constraint')
    parser.add_argument('--target_regions', '-q', action='store', dest='target_regions', type=str,
                        required=False, default=None, help="tab separated table with regions to align to")
    parser.add_argument('---un-banded', '-ub', action='store_false', dest='banded',
                        default=True, help='flag, turn off banding')
    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False,
                        default=1, type=int, help="number of jobs to run concurrently")
    parser.add_argument('--nb_files', '-n', action='store', dest='nb_files', required=False,
                        default=None, type=int, help="maximum number of reads to align")
    parser.add_argument("--kmer_size", action='store', dest="kmer_size", default=5, required=False,
                        help="size of kmers in fast5 file")
    parser.add_argument("--step_size", action='store', dest="step_size", default=10, required=False,
                        help="distance between positions of uncertainty")
    parser.add_argument("--alignment_file", action='store', dest="alignment_file", default=None, required=False,
                        help="a SAM/BAM with alignments of reads.  if set, cigar strings will be used only from this file")

    parser.add_argument("--validate", action='store', dest='validation_file', default=None, required=False,
                        help="validate an output file or directory (signalAlign will not be run)")

    args = parser.parse_args(args)

    # either: a) exactly one of [files_dir, fast5_glob], b) validation
    if not ((args.files_dir is None) ^ (args.fast5_glob is None)) or args.validation_file is not None:
        raise Exception("Unless validating, exactly one of --file_directory and --fast5_glob must be set")

    return args



def resolvePath(p):
    if p is None:
        return None
    elif p.startswith("/"):
        return p
    else:
        return os.path.abspath(p)


def organize_fast5s(fast5_locations):
    # gathered data
    fast5_to_read_id = dict()
    requires_event_calling = list()

    # examine each fast5
    for fast5 in fast5_locations:
        npr = NanoporeRead(fast5)
        success = npr.Initialize()
        read_id = npr.read_label
        fast5_id = os.path.basename(fast5)[:-6]
        fast5_to_read_id[fast5_id] = read_id
        if not success:
            requires_event_calling.append((fast5, read_id))
        npr.close()

    return fast5_to_read_id, requires_event_calling


def get_reference_sequence(ref_location, contig, start_pos, end_pos):
    # ensure faidx
    if not os.path.isfile("{}.fai".format(ref_location)): subprocess.check_call(['samtools', 'faidx', ref_location])
    if not os.path.isfile("{}.fai".format(ref_location)): pysam.faidx(ref_location)

    # use pysam
    with closing(pysam.FastaFile(ref_location)) as ref:
        return ref.fetch(reference=contig, start=start_pos, end=end_pos)



def validate_snp_directory(snp_directory, reference_sequence_path, alignment_file_location=None,
                           print_summary=False, move_files=False, make_plots=False):
    # prep
    all_consensus_identities = list()
    all_posterior_identities = list()
    all_lengths = list()
    all_summaries = list()
    # full_reference_sequence = get_first_sequence(reference_sequence_path).upper()
    files = glob.glob(os.path.join(snp_directory, "*.tsv"))
    if len(files) == 0:
        print("\n[singleNucleotideProbabilities] No files in {}\n".format(len(files), snp_directory))
        return

    print("\n[singleNucleotideProbabilities] Validating {} files in {}\n".format(len(files), snp_directory))

    # for rewriting
    if move_files:
        new_file_destination = snp_directory[:-1] if snp_directory.endswith("/") else snp_directory
        new_file_destination += ".filtered"
        if not os.path.isdir(new_file_destination): os.mkdir(new_file_destination)
        print("[singleNucleotideProbabilities] Copying files into {}\n".format(new_file_destination))
    else:
        new_file_destination = None

    # look at all files
    orig_cnt = 0
    prob_cnt = 0
    rewr_cnt = 0
    norw_cnt = 0

    for file in files:
        summary, problem = validate_snp_file(file, reference_sequence_path, print_summary=print_summary)
        consensus_identity = summary[CONSENSUS_IDENTITY_KEY]
        posterior_identity = summary[POSTERIOR_IDENTITY_KEY]
        length = summary[LENGTH_KEY]
        consensus_ratio = 1.0 * consensus_identity / length
        posterior_ratio = 1.0 * posterior_identity / length
        all_consensus_identities.append(consensus_identity)
        all_posterior_identities.append(posterior_identity)
        all_lengths.append(length)
        all_summaries.append(summary)

        if move_files:
            if problem:
                print("{}: not rewriting problem file".format(os.path.basename(file)))
                prob_cnt += 1
                continue
            # make a new file
            new_file = os.path.join(new_file_destination, os.path.basename(file))
            if consensus_ratio < VALID_IDENTITY_RATIO:
                print("{}:\tfile identity ratio is below threshold {}".format(
                    os.path.basename(file), VALID_IDENTITY_RATIO))
                norw_cnt += 1
            else:

                shutil.copyfile(file, new_file)
                orig_cnt += 1

    # printing results
    print("\n[singleNucleotideProbabilities] Summary of {} files:".format(len(files)))
    print("\tAVG Identity:             {}".format(np.mean(all_consensus_identities)))
    print("\tAVG Length:               {}".format(np.mean(all_lengths)))
    print("\tPosterior Identity Ratio: {}".format(1.0 * sum(all_posterior_identities) / sum(all_lengths)))
    print("\tConsensus Identity Ratio: {}".format(1.0 * sum(all_consensus_identities) / sum(all_lengths)))

    # forward and backward differences
    forwards = list(filter(lambda x: x[POSITIVE_STRAND_KEY] is not None and x[POSITIVE_STRAND_KEY], all_summaries))
    backwards = list(filter(lambda x: x[POSITIVE_STRAND_KEY] is not None and not x[POSITIVE_STRAND_KEY], all_summaries))
    if len(forwards) != 0:
        forward_consensus_iden = np.sum(list(map(lambda x: float(x[CONSENSUS_IDENTITY_KEY]), forwards))) / \
                                 np.sum(list(map(lambda x: x[LENGTH_KEY], forwards)))
        print("\tForward Identity Ratio:   {} ({})".format(forward_consensus_iden, len(forwards)))
    if len(backwards) != 0:
        backward_consensus_iden = np.sum(list(map(lambda x: float(x[CONSENSUS_IDENTITY_KEY]), backwards))) / \
                                 np.sum(list(map(lambda x: x[LENGTH_KEY], backwards)))
        print("\tBackward Identity Ratio:  {} ({})".format(backward_consensus_iden, len(backwards)))

    # analyze new ones if filtered
    if new_file_destination is not None:
        total_cnt = orig_cnt + prob_cnt + rewr_cnt + norw_cnt
        print("[singleNucleotideProbabilities] Summary of rewriting:")
        print("\tIncluded:               {} ({}%)".format(orig_cnt, int(100 * orig_cnt / total_cnt)))
        print("\tRem  - Below Threshold: {} ({}%)".format(norw_cnt, int(100 * norw_cnt / total_cnt)))
        print("\tRem  - SNP Prob:        {} ({}%)".format(prob_cnt, int(100 * prob_cnt / total_cnt)))
        print("\n[singleNucleotideProbabilities] rewrote {} ({}%) files.  Rerunning validation."
              .format(orig_cnt, int(100 * orig_cnt / len(files))))
        validate_snp_directory(new_file_destination, reference_sequence_path, print_summary=False, move_files=False,
                               make_plots=False)

    if make_plots:
        if alignment_file_location is None:
            print("\n[singleNucleotideProbabilities] Cannot make plots without alignment_file_location")
            return
        else:
            print("\n[singleNucleotideProbabilities] Generating plots for {} forward- and {} backwards-aligned reads"
                  .format(len(forwards), len(backwards)))

        #get the data we need
        missing_reads = 0
        min_read = 1
        min_signal = 1
        read_identities = get_read_identities_from_sam(alignment_file_location, reference_sequence_path)
        for read in all_summaries:
            if read[READ_NAME_KEY] in read_identities:
                read[NO_GAP_ALIGNED_IDENTITY_KEY] = read_identities[read[READ_NAME_KEY]][NO_GAP_ALIGNED_IDENTITY_KEY]
                min_read = min(min_read, read[NO_GAP_ALIGNED_IDENTITY_KEY])
                min_signal = min(min_signal, 1.0 * read[POSTERIOR_IDENTITY_KEY] / read[LENGTH_KEY])
            else:
                missing_reads += 1
                read[NO_GAP_ALIGNED_IDENTITY_KEY] = -1 # doesn't appear on the plot
                (forwards if read[POSITIVE_STRAND_KEY] else backwards).remove(read)
        if missing_reads > 0:
            print("[singleNucleotideProbabilities] Missing {} / {} reads in {}"
                  .format(missing_reads, len(all_summaries), alignment_file_location))

        # plot it
        import matplotlib.pyplot as plt

        # points
        fwd_read = list()
        fwd_signal = list()
        bwd_read = list()
        bwd_signal = list()
        for read in forwards:
            fwd_read.append(read[NO_GAP_ALIGNED_IDENTITY_KEY])
            fwd_signal.append(1.0 * read[POSTERIOR_IDENTITY_KEY] / read[LENGTH_KEY])
        for read in backwards:
            bwd_read.append(read[NO_GAP_ALIGNED_IDENTITY_KEY])
            bwd_signal.append(1.0 * read[POSTERIOR_IDENTITY_KEY] / read[LENGTH_KEY])

        f, = plt.plot(fwd_read, fwd_signal, 'b.', alpha=.25, label="Forward Strand")
        b, = plt.plot(bwd_read, bwd_signal, 'r.', alpha=.25, label="Backward Strand")
        plt.axis([max(0,min_read-(1-min_read)*.1),1,max(0,min_signal-(1-min_signal)*.1),1])
        plt.title("Identities")
        plt.legend(handles=[f, b])
        plt.xlabel("Read Alignment Identity")
        plt.ylabel("Signal Alignment Identity")
        plt.show()


def validate_snp_file(snp_file, reference_sequence_path, print_sequences=False, print_summary=False):
    identifier = os.path.basename(snp_file)
    consensus_sequence = list()
    all_probabilities = list()
    header_positions = list()
    header_characters = list()
    first_pos = None
    last_pos = None
    problem = False
    duplicated_positions = 0
    unspecified_positions = 0

    forward = None
    read_name = None
    fast5_name =  None
    contig = None
    with open(snp_file, 'r') as snp:
        for line in snp:
            if line.startswith("##"):
                if "strand" in line:
                    forward = "template" in line
                elif "read_id" in line:
                    read_name = line.split(":")[1].strip()
                elif "fast5_input" in line:
                    fast5_name = line.split(":")[1].strip()
                elif "contig" in line:
                    contig = line.split(":")[1].strip()
                continue
            elif line.startswith("#"):
                identifier = "{}/{}".format(read_name, fast5_name)
                line = line.split("\t")
                i = 0
                for l in line:
                    if l.startswith("p"):
                        header_positions.append(i)
                        header_characters.append(l[1].upper())
                    i += 1
            else:
                line = line.split("\t")
                # positions
                pos = int(line[1])
                # set first_position (for reference matching)
                if first_pos is None: first_pos = pos
                # for cases where positions are duplicated or missing
                if last_pos is not None:
                    if last_pos >= pos:
                        duplicated_positions += 1
                        continue
                    while last_pos + 1 < pos:
                        probabilities = dict()
                        for char in header_characters:
                            probabilities[char] = 0.0
                        all_probabilities.append(probabilities)
                        consensus_sequence.append("-")
                        unspecified_positions += 1
                        last_pos += 1
                last_pos = pos
                #probabilities
                probabilities = dict()
                #consensus
                max_prob = -1.0
                max_prob_idx = None
                idx = 0
                for pos in header_positions:
                    prob = float(line[pos].strip())
                    probabilities[header_characters[idx]] = prob
                    if prob > max_prob:
                        max_prob = prob
                        max_prob_idx = idx
                    idx += 1
                consensus_sequence.append(header_characters[max_prob_idx])
                all_probabilities.append(probabilities)

    # get sequences
    consensus_sequence = "".join(consensus_sequence)
    reference_sequence = get_reference_sequence(reference_sequence_path, contig,
                                                min(first_pos, last_pos),max(first_pos, last_pos)+1).upper()

    # this is our quality metric
    length = 0
    identity = 0
    posterior_identity = 0.0
    for c, r, p in zip(consensus_sequence, reference_sequence, all_probabilities):
        length += 1
        if c == r: identity += 1
        posterior_identity += p[r]

    # for separating results
    if print_sequences or print_summary: print("")

    # sanity check
    if duplicated_positions > 0:
        print("{}: Found {} duplicated positions!"
              .format(identifier, duplicated_positions))
    if unspecified_positions * 100 > length:
        print("{}: Found {} unspecified positions ({}% of total length)"
              .format(identifier, unspecified_positions, int(100.0 * unspecified_positions / length)))
        problem = True

    # printing full sequences
    if print_sequences:
        print("{}: Whole Sequences:".format(identifier))
        print("\treference:  {}".format(reference_sequence))
        print("\tconsensus:  {}".format(consensus_sequence))
        # print("\tcomplement: {}".format(reverse_complement(consensus_sequence, complement=True, reverse=False)))

    if print_summary:
        strand_char = " " if forward is None else ("+" if forward else "-")
        print("%s:\tdirection: %s\tlength:%8d\t\tconsensus identity:%8d (%8f)\t\tposterior_identity:%8d (%8f)" %
              (identifier, strand_char, length, identity, 1.0*identity/length,
               int(posterior_identity), posterior_identity/length))

    summary = {
        READ_NAME_KEY: read_name,
        FAST5_NAME_KEY: fast5_name,
        LENGTH_KEY: length,
        CONSENSUS_IDENTITY_KEY: identity,
        POSTERIOR_IDENTITY_KEY: posterior_identity,
        POSITIVE_STRAND_KEY: forward
    }
    return summary, problem


def get_read_identities_from_sam(alignment_location, reference_location):

    # what we track
    read_to_identity = dict()

    # ref mgmt
    reference_string = None
    ref_contig = None

    with closing(pysam.AlignmentFile(alignment_location, 'rb' if alignment_location.endswith("bam") else 'r')) as aln:
        for aligned_segment in aln:
            if not aligned_segment.is_secondary and not aligned_segment.is_unmapped:
                contig = aln.getrname(aligned_segment.rname)
                if ref_contig is None or contig != ref_contig:
                    reference_string = get_reference_sequence(reference_location, contig, None, None)
                query_name = aligned_segment.qname
                read = aligned_segment.query_sequence
                aligned_tuples = aligned_segment.get_aligned_pairs()

                # counts
                length = 0
                no_gap_len = 0
                iden = 0.0
                no_gap_iden = 0.0

                # iterate
                for tuple in aligned_tuples:
                    # get characters
                    if tuple[0] is None: read_char = '-'
                    elif len(read) <= tuple[0]:
                        print("Alignment pos {} exceeds read length {}: {}".format(tuple[0], len(read), query_name))
                        read_char = '-'
                    else: read_char = read[tuple[0]]
                    if tuple[1] is None: ref_char = '-'
                    else: ref_char = reference_string[tuple[1]]

                    # upper
                    read_char = read_char.upper()
                    ref_char = ref_char.upper()

                    #counts
                    length += 1
                    if read_char == ref_char: iden += 1
                    if read_char != '-' and ref_char != '-':
                        no_gap_len += 1
                        if read_char == ref_char: no_gap_iden += 1

                # save
                read_to_identity[query_name] = {
                    READ_NAME_KEY: query_name,
                    ALIGNED_IDENTITY_KEY: iden / length,
                    NO_GAP_ALIGNED_IDENTITY_KEY: no_gap_iden / no_gap_len
                }

    return read_to_identity


def variant_caller(work_queue, done_queue, service_name="variant_caller"):
    # prep
    total_handled = 0
    failure_count = 0

    #catch overall exceptions
    try:
        for f in iter(work_queue.get, 'STOP'):
            # catch exceptions on each element
            try:
                work = CallMethylation(**f)
                work.write()

            except Exception as e:
                # get error and log it
                message = "{}:{}".format(type(e), str(e))
                error = "{} '{}' failed with: {}".format(service_name, current_process().name, message)
                print("[{}] ".format(service_name) + error)
                done_queue.put(error)
                failure_count += 1

            # increment total handling
            total_handled += 1

    except Exception as e:
        # get error and log it
        message = "{}:{}".format(type(e), str(e))
        error = "{} '{}' critically failed with: {}".format(service_name, current_process().name, message)
        print("[{}] ".format(service_name) + error)
        done_queue.put(error)

    finally:
        # logging and final reporting
        print("[%s] '%s' completed %d calls with %d failures"
              % (service_name, current_process().name, total_handled, failure_count))
        done_queue.put("{}:{}".format(TOTAL_KEY, total_handled))
        done_queue.put("{}:{}".format(FAILURE_KEY, failure_count))


def event_detection(work_queue, done_queue, alignment_file, event_detection_strategy=None, event_detection_params=None,
                    service_name="event_detection"):
    # prep
    total_handled = 0
    failure_count = 0

    #catch overall exceptions
    try:
        for fast5, read_id in iter(work_queue.get, 'STOP'):
            # catch exceptions on each element
            try:
                nucleotide, qualities, hardcode_front, hardcode_back = get_full_nucleotide_read_from_alignment(
                    alignment_file, read_id)
                success, _, _ = generate_events_and_alignment(fast5, nucleotide, nucleotide_qualities=qualities,
                                                              save_to_fast5=False, overwrite=False)
                if not success:
                    failure_count += 1

            except Exception as e:
                # get error and log it
                message = "{}:{}".format(type(e), str(e))
                error = "{} '{}' failed with: {}".format(service_name, current_process().name, message)
                print("[{}] ".format(service_name) + error)
                done_queue.put(error)
                failure_count += 1

            # increment total handling
            total_handled += 1

    except Exception as e:
        # get error and log it
        message = "{}:{}".format(type(e), str(e))
        error = "{} '{}' critically failed with: {}".format(service_name, current_process().name, message)
        print("[{}] ".format(service_name) + error)
        done_queue.put(error)

    finally:
        # logging and final reporting
        print("[%s] '%s' completed %d calls with %d failures"
              % (service_name, current_process().name, total_handled, failure_count))
        done_queue.put("{}:{}".format(TOTAL_KEY, total_handled))
        done_queue.put("{}:{}".format(FAILURE_KEY, failure_count))


def discover_single_nucleotide_probabilities(args, working_folder, kmer_length, reference_location,
                                             list_of_fast5s, alignment_args, workers, step_size,
                                             output_directory=None, use_saved_alignments=True, save_alignments=True):

    # read fast5s and extract read ids
    fast5_to_read, requires_event_calling = organize_fast5s(list_of_fast5s)
    print("[info] built map of fast5 identifiers to read ids with {} elements".format(len(fast5_to_read)))

    # promethION fast5s do not come with events, need to do event calling
    if len(requires_event_calling) > 0:
        # log and error checking
        print("[info] {}/{} fast5s need event detection and read alignment".format(
            len(requires_event_calling), len(fast5_to_read)))
        if args.alignment_file is None:
            print("[error] Cannot perform event detection without alignmentfile", file=sys.stderr)
            sys.exit(1)

        # prep for event detection
        event_detection_service_args = {
            'alignment_file': args.alignment_file,
            'event_detection_strategy': None,
            'event_detection_params': None
        }

        # do event detection
        total, failure, messages = run_service(event_detection, requires_event_calling, {}, 'fast5', workers,
                                               service_arguments=event_detection_service_args)

        # loggit and continue
        print("[info] {}/{} fast5s successfully had events detected".format(total - failure, total))





    # do alignment and calling for each step
    for s in range(step_size):
        # log start
        print("\n[info] starting step %d / %d" % (s + 1, step_size))

        # prep (this is where alignments will live if they have been saved)
        saved_step_dir = os.path.join(working_folder.path, "step_{}".format(s))

        # do or get alignments
        if use_saved_alignments and os.path.isdir(saved_step_dir):
            # getting alignments
            alignments = [x for x in glob.glob(os.path.join(saved_step_dir, "*.tsv")) if os.stat(x).st_size != 0]
            print("[info] using {} saved alignments (of {} fast5s) in {}".format(
                len(alignments), len(list_of_fast5s), saved_step_dir))
        else:
            # file locations
            sub_fasta_path = working_folder.add_file_path(
                "ref_ambig.s{}.o{}.{}".format(s, step_size, os.path.basename(reference_location)))

            # build reference
            substitution_ref = replace_periodic_reference_positions(reference_location, sub_fasta_path, step_size, s)
            alignment_args['forward_reference'] = substitution_ref
            # run alignment
            multithread_signal_alignment(alignment_args, list_of_fast5s, workers)

            # get alignments
            alignments = [x for x in glob.glob(os.path.join(working_folder.path, "*.tsv")) if os.stat(x).st_size != 0]

        # sanity check
        alignment_count = len(alignments)
        if alignment_count == 0:
            print("[error] Didn't find any alignment files here {}".format(working_folder.path))
            sys.exit(1)
        else:
            print("[info] Found %d alignment files (%d input fast5s) here %s" %
                  (alignment_count, len(list_of_fast5s), working_folder.path))

        # prep to get marginals
        marginal_probability_prefix = working_folder.add_file_path("marginals.{step}".format(step=s))
        proposal_args = {
            "sequence": None,
            "out_file_prefix": marginal_probability_prefix,
            "step_size": step_size,
            "step_offset": s,
            "degenerate_type": alignment_args["degenerate"],
            "kmer_length": kmer_length
        }

        # get marginals
        print("[info] running variant_caller on %d alignments files with %d workers" % (alignment_count, workers))
        run_service(variant_caller, alignments, proposal_args, "alignment_file", workers)

        # remove or save old alignments
        files = glob.glob(working_folder.path + "*.tsv")
        if save_alignments:
            if not os.path.isdir(saved_step_dir): os.mkdir(saved_step_dir)
            print("[info] saving {} alignment files into {}".format(len(files), saved_step_dir))
            for f in files:
                os.rename(f, os.path.join(saved_step_dir, os.path.basename(f)))
        else:
            print("[info] deleting {} alignment files".format(len(files)))
            for f in files:
                os.remove(f)
        print("[info] step %d completed\n" % s)

    # per fast5, we want to coalesce all step calling into one file (named by read)
    if output_directory is None:
        output_directory = os.path.join(working_folder.path, "reads")
    if not os.path.isdir(output_directory): os.mkdir(output_directory)
    print("[info] writing output to {}".format(output_directory))
    output_files = list()

    # iterate over input fast5s
    for fast5_id, read_id in fast5_to_read.items():
        # get files
        files = glob.glob(os.path.join(working_folder.path,
                                       "marginals*{}*{}".format(read_id, CallMethylation.FILE_EXTENSION)))
        if len(files) != step_size:
            print("[error] input fast5 '{}' for read {} yielded {} output files, expected {}".format(
                fast5_id, read_id, len(files), step_size))
            if len(files) == 0:
                continue

        # read all lines in all files
        output_lines = list()
        template_count = 0
        complement_count = 0
        contigs = set()
        for file in files:
            if "forward" in file.split("."):
                template_count += 1
            if "backward" in file.split("."):
                complement_count += 1
            with open(file, 'r') as input:
                for line in input:
                    line = line.split("\t")
                    contigs.add(line[0])
                    output_lines.append(line)
        # sort based on position
        output_lines.sort(key=lambda x: int(x[1]))

        # template/complement and sanity checks
        reverse = complement_count > template_count
        strand_identifier = "complement" if reverse else "template"
        if template_count == 0 and complement_count == 0:
            print("[warn] {}: could not infer template or complement".format(fast5_id))
            strand_identifier = "unknown"
        if template_count != 0 and complement_count != 0:
            print("[warn] {}: got {} template and {} complement calls after variant calling"
                  .format(fast5_id, template_count, complement_count))
            strand_identifier += " (by majority)"

        # write output
        output_filename = "{}.tsv".format(read_id)
        output_file = os.path.join(output_directory, output_filename)
        with open(output_file, 'w') as output:
            output.write("## fast5_input: {}.fast5\n".format(fast5_id))
            output.write("## read_id: {}\n".format(read_id))
            output.write("## contig: {}\n".format(",".join(contigs)))
            output.write("## strand: {}\n".format(strand_identifier))
            output.write("#CHROM\tPOS\tpA\tpC\tpG\tpT\n")
            for line in output_lines:
                if not reverse:
                    line = [line[0], line[1], line[3], line[4], line[5], line[6]]
                else:
                    line = [line[0], line[1], line[6], line[5], line[4], line[3]]
                output.write("\t".join(line) + "\n")
        #save
        output_files.append(output_file)

    # document and return
    print("[info] wrote {} output files ({} input fast5s) in {}"
          .format(len(output_files), len(fast5_to_read), output_directory))
    return output_files


def main(args):
    # parse args
    args = parse_args(args)

    command_line = " ".join(sys.argv[:])
    print("[singleNucleotideProbabilities] Command Line: {cmdLine}\n".format(cmdLine=command_line), file=sys.stderr)

    # first: see if we want to validate and return
    if args.validation_file is not None:
        if os.path.isfile(args.validation_file):
            validate_snp_file(args.validation_file, args.ref, print_sequences=True, print_summary=True)
        elif os.path.isdir(args.validation_file):
            validate_snp_directory(args.validation_file, args.ref, print_summary=False, move_files=False,
                                   make_plots=True, alignment_file_location=args.alignment_file)
        else:
            print("[error] got invalid validation location: {}".format(args.validation_file))
        return 0

    # get absolute paths to inputs
    args.files_dir           = resolvePath(args.files_dir)
    args.fast5_glob          = resolvePath(args.fast5_glob)
    args.ref                 = resolvePath(args.ref)
    args.out                 = resolvePath(args.out)
    args.in_T_Hmm            = resolvePath(args.in_T_Hmm)
    args.in_C_Hmm            = resolvePath(args.in_C_Hmm)
    args.templateHDP         = resolvePath(args.templateHDP)
    args.complementHDP       = resolvePath(args.complementHDP)
    args.target_regions      = resolvePath(args.target_regions)
    args.alignment_file      = resolvePath(args.alignment_file)

    # assert integers
    args.step_size = int(args.step_size)
    args.kmer_size = int(args.kmer_size)

    # get input glob
    input_glob = args.fast5_glob if args.fast5_glob is not None else os.path.join(args.files_dir, "*.fast5")

    start_message = """
#   Single Nucleotide Probabilities
#
#   Aligning files matching: {inputGlob}
#   Aligning to reference: {reference}
#   Aligning maximum of {nbFiles} files
#   Using model: {model}
#   Using banding: {banding}
#   Aligning to regions in: {regions}
#   Non-default template HMM: {inThmm}
#   Non-default complement HMM: {inChmm}
#   Template HDP: {tHdp}
#   Complement HDP: {cHdp}
#   Kmer size: {kmerSize}
#   Step size: {stepSize}
#   Alignment File: {alignmentFile}
    """.format(inputGlob=input_glob, reference=args.ref, banding=args.banded, nbFiles=args.nb_files,
               inThmm=args.in_T_Hmm, inChmm=args.in_C_Hmm, model=args.stateMachineType, regions=args.target_regions,
               tHdp=args.templateHDP, cHdp=args.complementHDP, kmerSize=args.kmer_size, stepSize=args.step_size,
               alignmentFile=args.alignment_file)
    print(start_message, file=sys.stdout)

    # prep
    if not os.path.isdir(args.out): os.mkdir(args.out)

    # get fast5 locations and prune
    fast5s = glob.glob(input_glob)
    if args.nb_files is not None and args.nb_files < len(fast5s):
        print("[singleNucleotideProbabilities] pruning {} fast5 files down to configured max {}".format(len(fast5s), args.nb_files))
        shuffle(fast5s)
        fast5s = fast5s[:args.nb_files]

    # get the (input) reference sequence
    if not os.path.isfile(args.ref):
        print("[singleNucleotideProbabilities] Did not find valid reference file", file=sys.stderr)
        sys.exit(1)

    # make a working folder in the specified directory
    temp_folder = FolderHandler()
    temp_dir_path = temp_folder.open_folder(os.path.join(args.out, "tempFiles_errorCorrection"))

    # alignment args are the parameters to the HMM/HDP model, and don't change
    alignment_args = {
        # "path_to_EC_refs": None,
        "destination": temp_dir_path,
        "stateMachineType": args.stateMachineType,
        "bwa_reference": args.ref,
        "in_templateHmm": args.in_T_Hmm,
        "in_complementHmm": args.in_C_Hmm,
        "in_templateHdp": args.templateHDP,
        "in_complementHdp": args.complementHDP,
        "threshold": args.threshold,
        "diagonal_expansion": args.diag_expansion,
        "constraint_trim": args.constraint_trim,
        "target_regions": None,
        "degenerate": getDegenerateEnum("variant"),
        "alignment_file": args.alignment_file,
        'track_memory_usage': False,
        'get_expectations': False
    }

    # get the sites that have proposed edits
    print("\n\n[singleNucleotideProbabilities] scanning for proposals with %d fast5s" % len(fast5s))
    output_files = discover_single_nucleotide_probabilities(args, temp_folder, args.kmer_size, args.ref, fast5s, alignment_args,
                                                            args.nb_jobs, args.step_size, output_directory=args.out)
    print("\n[singleNucleotideProbabilities] got {} output files:".format(len(output_files)))
    i = 0
    for output_file in output_files:
        print("\t{}".format(output_file))
        i += 1
        if i > 10 and len(output_files) > 10:
            print("\t...")
            break

    #validation
    if len(output_files) != 0:
        validate_snp_directory(os.path.dirname(output_files[0]), args.ref, print_summary=True)

    print("\n\n[singleNucleotideProbabilities] fin\n")

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
