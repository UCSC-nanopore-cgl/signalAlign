#!/usr/bin/env python
"""Train HMMs for alignment of signal data from the MinION
"""

import sys
import os
import glob
import pandas as pd
import numpy as np

from timeit import default_timer as timer
from argparse import ArgumentParser
from shutil import copyfile
from subprocess import check_call
import shutil
import tempfile
from py3helpers.utils import create_dot_dict, merge_lists, all_string_permutations, save_json, load_json, \
    count_lines_in_file, merge_dicts, list_dir
from py3helpers.multiprocess import *

from signalalign.signalAlignment import multithread_signal_alignment_samples, create_signalAlignment_args, \
    SignalAlignSample
from signalalign.hiddenMarkovModel import HmmModel, parse_assignment_file, parse_alignment_file, read_in_alignment_file
from signalalign.utils.fileHandlers import FolderHandler
from signalalign.utils.parsers import read_fasta
from signalalign.utils.sequenceTools import get_motif_kmers, get_sequence_kmers, CustomAmbiguityPositions
from signalalign.build_alignments import generate_top_n_kmers_from_sa_output

def make_master_assignment_table(list_of_assignment_paths, min_probability=0.0, full=False):
    """Create a master assignment table from a list of assignment paths

    :param list_of_assignment_paths: list of all paths to assignment.tsv files to concat
    :param min_probability: minimum probabilty to keep
    :return: pandas DataFrame of all assignments
    """
    assignment_dfs = []
    for f in list_of_assignment_paths:
        assignment_dfs.append(get_assignment_table(f, min_probability, full))
    return pd.concat(assignment_dfs, ignore_index=True)


def multiprocess_make_master_assignment_table(list_of_assignment_paths, min_probability=0.0, full=False,
                                              worker_count=1):
    """Create a master assignment table from a list of assignment paths

    :param list_of_assignment_paths: list of all paths to assignment.tsv files to concat
    :param min_probability: minimum probabilty to keep
    :return: pandas DataFrame of all assignments
    """
    extra_args = {"min_probability": min_probability,
                  "full": full}
    service = BasicService(get_assignment_table, service_name="multiprocess_make_master_assignment_table")
    total, failure, messages, output = run_service(service.run, list_of_assignment_paths,
                                                   extra_args, ["file_path"], worker_count)
    return pd.concat(output, ignore_index=True)


def get_assignment_table(file_path, min_probability, full):
    if full:
        data = parse_alignment_file(file_path)
    else:
        data = parse_assignment_file(file_path)
    return data.loc[data['prob'] >= min_probability]


def multiprocess_make_kmer_assignment_tables(list_of_assignment_paths, kmers, strands, min_probability=0.0,
                                             verbose=True, full=False, max_assignments=10,
                                             worker_count=1):
    """Create a master assignment tables from a list of assignment paths

    :param kmers:
    :param strands:
    :param verbose:
    :param full:
    :param max_assignments:
    :param worker_count:
    :param list_of_assignment_paths: list of all paths to assignment.tsv files to concat
    :param min_probability: minimum probabilty to keep
    :return: pandas DataFrame of all assignments
    """
    # just in case we get a set
    kmers = list(kmers)
    # Multiprocess reading in files
    extra_args = {"min_probability": min_probability,
                  "full": full,
                  "kmers": kmers}
    service = BasicService(get_assignment_kmer_tables, service_name="get_assignment_kmer_tables")
    total, failure, messages, output = run_service(service.run, list_of_assignment_paths,
                                                   extra_args, ["file_path"], worker_count)
    kmer_tables = [(pd.concat(x, ignore_index=True), kmers[i]) for i, x in enumerate(zip(*output))]
    # Multiprocess sorting each kmer table
    extra_args = {"strands": strands,
                  "verbose": verbose,
                  "max_assignments": max_assignments}
    service = BasicService(sort_dataframe_wrapper, service_name="sort_dataframe_wrapper")
    total, failure, messages, output = run_service(service.run, kmer_tables,
                                                   extra_args, ["data_table", "kmer"], worker_count)

    return pd.concat(output, ignore_index=True)


def get_assignment_kmer_tables(file_path, kmers, min_probability, full):
    if full:
        data = parse_alignment_file(file_path)
    else:
        data = parse_assignment_file(file_path)
    data = data.loc[data['prob'] >= min_probability]
    all_kmers = []
    for k in kmers:
        all_kmers.append(data.loc[data.kmer == k])
    return all_kmers


def sort_dataframe_wrapper(data_table, kmer, max_assignments=10, verbose=False, strands=('t', 'c')):
    assert len(strands) > 0, \
        "strands must be a list and not be empty. strands: {}".format(strands)
    final_output = []
    if data_table.empty and verbose:
        print("missing kmer {}, continuing".format(kmer))
        final_output = data_table
    else:
        for strand in strands:
            by_strand = data_table.loc[data_table['strand'] == strand]
            kmer_assignments = by_strand.sort_values(['prob'], ascending=0)[:max_assignments]

            if len(kmer_assignments) < max_assignments and verbose:
                print("WARNING didn't find {max} requested assignments for {kmer} only found {found}"
                      "".format(max=max_assignments, kmer=kmer, found=len(kmer_assignments)))
            final_output.append(kmer_assignments)
        final_output = pd.concat(final_output, ignore_index=True)
    return final_output


def get_kmers(kmer_len, alphabet="ATGC", motifs=None, reference=None):
    """Get all kmers from the reference sequence, all possible kmers from an alphabet and/or
        all kmers that cover a modified nucelotide

    :param kmer_len: length of kmer
    :param alphabet: alphabet
    :param motifs: if you want specific motifs to be looked at
    :param reference: reference sequence to get kmers from
    :return: set of desired kmers
    """
    kmers = set()
    # if motifs is present, process for all motifs with modified base
    if motifs is not None:
        for motif in motifs:
            kmers |= get_motif_kmers(motif, kmer_len, alphabet=alphabet)
    # if we want to limit kmers which were seen in reference sequence
    if reference is not None:
        for _, _, sequence in read_fasta(reference):
            kmers |= get_sequence_kmers(sequence, k=kmer_len, rev_comp=True)
    else:
        kmers |= {x for x in all_string_permutations(alphabet, length=kmer_len)}

    return kmers


def make_master_full_alignments_table(files, min_probability=0.0):
    """Get the forward mapped reads from the full output of signalAlign
    :param files: path to directory of output from signalalign
    :param min_probability: minimum probability
    """
    assignment_dfs = []
    for f in files:
        data = read_in_alignment_file(f)
        assignment_dfs.append(data.loc[data['posterior_probability'] >= min_probability])
    return pd.concat(assignment_dfs, ignore_index=True)


def write_and_generate_build_alignments_positions(full_alignments_dirs, positions_file, outpath, min_probability=0.0,
                                                  max_assignments=10, verbose=True):
    """Writes the built alignments from 'full' assignments given that the kmers cover the positions file
    :param full_alignments_dirs: list of assignment directories
    :param positions_file: positions file
    :param min_probability: minimum probability for kmer assignments
    :param max_assignments: maximum number of kmer assignments
    :param verbose: boolean option for print statements
    :param outpath: where to write the output file
    :return: all kmer assignemnts with at least the min probability and no more than max assignments
    """
    data = generate_build_alignments_positions(full_alignments_dirs, positions_file, min_probability=min_probability,
                                               max_assignments=max_assignments, verbose=verbose)
    data.rename(index=str, columns={"path_kmer": "kmer", "descaled_event_mean": "level_mean",
                                    "posterior_probability": "prob"},
                inplace=True)
    data = data[["kmer", "strand", "level_mean", "prob"]]
    data.to_csv(outpath, sep='\t', header=False, index=False)
    return data


def generate_build_alignments_positions(full_alignments_dirs, positions_file, min_probability=0.0,
                                        max_assignments=10, verbose=True):
    """Generate built alignments from 'full' assignments given that the kmers cover the positions file
    :param full_alignments_dirs: list of assignment directories
    :param positions_file: positions file
    :param min_probability: minimum probability for kmer assignments
    :param max_assignments: maximum number of kmer assignments
    :param verbose: boolean option for print statements
    :return: all kmer assignemnts with at least the min probability and no more than max assignments
    """
    all_data = []
    positions_data = CustomAmbiguityPositions.parseAmbiguityFile(positions_file)

    for full_alignments_dir in full_alignments_dirs:
        # Get forward data and backward data from directory
        forward_data = make_master_full_alignments_table(list_dir(full_alignments_dir, ext="forward.tsv"),
                                                         min_probability=min_probability)

        backward_data = make_master_full_alignments_table(list_dir(full_alignments_dir, ext="backward.tsv"),
                                                          min_probability=min_probability)

        forward_data = build_alignments_positions(forward_data, positions_data, forward=True)
        backward_data = build_alignments_positions(backward_data, positions_data, forward=False)
        all_data.extend([forward_data, backward_data])

    final_data = filter_top_n_kmers(pd.concat(all_data, ignore_index=True), max_n=max_assignments, verbose=verbose)
    return final_data


def get_kmers_covering_positions(full_sa_output_pd, positions, kmer_len):
    """Get all kmers that cover a list of positions. It is assumed that the contig/chromosome is the same
    :param full_sa_output_pd: data that we want to filter for kmers that cover certain positions
    :param positions: list of positions we want to include
    :param kmer_len: kmer_len
    :return: kmer alignment data that covers the positions passed in
    """

    all_positions_to_keep = set(
        merge_lists([list(range(x - (kmer_len - 1), x + 1)) for x in positions]))

    keepers = full_sa_output_pd[[x in all_positions_to_keep for x in full_sa_output_pd["reference_index"]]]
    return keepers


def build_alignments_positions(full_sa_output_pd, positions_data, forward):
    """Convert assignments to alignment line format for HDP training given the kmer covers certain reference positions

    Filter assignments on a minimum probability, read strand, and a max number of kmer assignments

    :param full_sa_output_pd: all data with "contig" "strand" and "reference_index"
    :param positions_data: positions data with "contig", "strand" and "position" fields
    :param forward: boolean option for dealing with forward positions or backward strand positions
    """
    # loop through for each strand in the assignments
    final_output = []

    kmer_len = len(full_sa_output_pd["path_kmer"][0])
    if forward is True:
        strand_pairs = [["+", "t"], ["-", "c"]]
    else:
        strand_pairs = [["-", "t"], ["+", "c"]]

    for contig in set(full_sa_output_pd["contig"]):
        data_by_contig = full_sa_output_pd[full_sa_output_pd['contig'] == contig]
        contig_positions = positions_data[positions_data["contig"] == contig]

        for strand_pair in strand_pairs:
            stranded_positions = contig_positions[contig_positions["strand"] == strand_pair[0]]

            data_by_strand = data_by_contig[data_by_contig['strand'] == strand_pair[1]]
            final_output.append(get_kmers_covering_positions(data_by_strand,
                                                             stranded_positions["position"],
                                                             kmer_len=kmer_len))
    return pd.concat(final_output, ignore_index=True)


def filter_top_n_kmers(kmer_data, max_n=10, verbose=True):
    """Filter out only the top N most probable kmers
    :param verbose: boolean option to print warnings when there are not enough events for a given kmer
    :param kmer_data: pd dataframe with "path_kmer" and "posterior_probability"
    :param max_n: max number of kmers to keep
    """
    output = []
    kmers = set(kmer_data["path_kmer"])
    for kmer in kmers:
        specific_kmer_data = kmer_data[kmer_data["path_kmer"] == kmer]
        top_n_kmers = specific_kmer_data.sort_values(["posterior_probability"], ascending=False)[:max_n]
        output.append(top_n_kmers)
        if len(top_n_kmers) < max_n and verbose:
            print("WARNING didn't find {max} requested assignments for {kmer} only found {found}"
                  "".format(max=max_n, kmer=kmer, found=len(top_n_kmers)))

    return pd.concat(output, ignore_index=True)


def generate_buildAlignments(assignments_pd, kmer_list, max_assignments=10, strands=('t', 'c'), verbose=False):
    """Convert assignments to alignment line format for HDP training.

    Filter assignments on a minimum probability, read strand, and a max number of kmer assignments

    :param assignments_pd: giant assignments pandas data table
    :param verbose: option to print update statements
    :param strands: 't' or 'c' representing template or complement strand of read
    :param kmer_list: list of kmers to write to alignment file
    :param max_assignments: max number of assignments to process for each kmer
    :param min_probability: the minimum probability to use for assigning kmers
    """
    # loop through for each strand in the assignments
    assert len(strands) > 0, \
        "strands must be a list and not be empty. strands: {}".format(strands)
    final_output = []
    for strand in strands:
        by_strand = assignments_pd.loc[assignments_pd['strand'] == strand]
        by_strand.sort_values(by='kmer', inplace=True)
        by_strand.set_index(keys=['kmer'], drop=False, inplace=True)
        for k in kmer_list:
            final_output.append(kmer_assignments_wrapper(by_strand, k, max_assignments, verbose))
    return pd.concat(final_output, ignore_index=True)


def multiprocess_generate_buildAlignments(assignments_pd, kmer_list, max_assignments=10, strands=('t', 'c'),
                                          verbose=False, worker_count=1):
    """Convert assignments to alignment line format for HDP training.

    Filter assignments on a minimum probability, read strand, and a max number of kmer assignments

    :param assignments_pd: giant assignments pandas data table
    :param verbose: option to print update statements
    :param strands: 't' or 'c' representing template or complement strand of read
    :param kmer_list: list of kmers to write to alignment file
    :param max_assignments: max number of assignments to process for each kmer
    :param min_probability: the minimum probability to use for assigning kmers
    """
    # loop through for each strand in the assignments
    assert len(strands) > 0, \
        "strands must be a list and not be empty. strands: {}".format(strands)
    final_output = []
    for strand in strands:
        by_strand = assignments_pd.loc[assignments_pd['strand'] == strand]
        by_strand.set_index(keys=['kmer'], drop=False, inplace=True)
        extra_args = {"by_strand": by_strand,
                      "max_assignments": max_assignments,
                      "verbose": verbose}
        service = BasicService(kmer_assignments_wrapper, service_name="multiprocess_generate_buildAlignments")
        total, failure, messages, output = run_service(service.run, kmer_list,
                                                       extra_args, ["k"], worker_count)

        final_output.extend(output)
    return pd.concat(final_output, ignore_index=True)


def kmer_assignments_wrapper(by_strand, k, max_assignments=10, verbose=False):
    kmer_assignments = by_strand.loc[by_strand.kmer == k]
    if kmer_assignments.empty and verbose:
        print("missing kmer {}, continuing".format(k))
    kmer_assignments = kmer_assignments.sort_values(['prob'], ascending=0)[:max_assignments]
    if len(kmer_assignments) < max_assignments and verbose:
        print("WARNING didn't find {max} requested assignments for {kmer} only found {found}"
              "".format(max=max_assignments, kmer=k, found=len(kmer_assignments)))
    return kmer_assignments


def generate_buildAlignments2(assignments_pd, kmer_list, max_assignments=10, strands=('t', 'c'), verbose=False):
    """Convert assignments to alignment line format for HDP training.

    Filter assignments on a minimum probability, read strand, and a max number of kmer assignments

    :param assignments_pd: giant assignments pandas data table
    :param verbose: option to print update statements
    :param strands: 't' or 'c' representing template or complement strand of read
    :param kmer_list: list of kmers to write to alignment file
    :param max_assignments: max number of assignments to process for each kmer
    """
    # loop through for each strand in the assignments
    assert len(strands) > 0, \
        "strands must be a list and not be empty. strands: {}".format(strands)
    final_output = []
    for strand in strands:
        by_strand = assignments_pd.loc[assignments_pd['strand'] == strand]

        for k in kmer_list:
            kmer_assignments = by_strand.loc[by_strand['kmer'] == k]
            if kmer_assignments.empty and verbose:
                print("missing kmer {}, continuing".format(k))
                continue
            kmer_assignments = kmer_assignments.sort_values(['prob'], ascending=0)[:max_assignments]
            final_output.append(kmer_assignments)
            if len(kmer_assignments) < max_assignments and verbose:
                print("WARNING didn't find {max} requested assignments for {kmer} only found {found}"
                      "".format(max=max_assignments, kmer=k, found=len(kmer_assignments)))
    return pd.concat(final_output, ignore_index=True)


class CreateHdpTrainingData(object):
    """Process the assignment files created from SignalAlign for the HDP distribution estimation"""

    def __init__(self, samples, out_file_path, template=True, complement=False, verbose=True, alphabet="ATGC",
                 jobs=1):
        """
        Control how each kmer/event assignment is processed given a set of samples and the parameters associated with
        each sample

        :param out_file_path: path to ouptut file
        :param jobs: number of jobs to multiprocess with
        :param samples: SignalAlignSamples
        :param template: generate kmers for template read strand: default: True
        :param complement: generate kmers for complement read strand: default: True
        :param verbose: option to print update statements
        :param alphabet: alphabet of sequencing experiment
        """
        self.jobs = jobs
        self.strands = []
        if template:
            self.strands.append('t')
        if complement:
            self.strands.append('c')
        assert self.strands != [], 'template or complement need to be set to True. ' \
                                   'complement: {}, template: {}'.format(complement, template)

        for sample in samples:
            assert isinstance(sample, SignalAlignSample)

        self.alphabet = alphabet
        self.samples = samples
        self.out_file_path = out_file_path
        self.template = template
        self.complement = complement
        self.verbose = verbose
        self.k = 0
        self.n_assignments = 0

    def write_hdp_training_file(self, verbose=False):
        """Write a hdp training file to a specified location"""
        # final_output = []
        for sample in self.samples:
            print("Filtering and gathering {} assignment.tsv files".format(sample.name))
            if len(sample.analysis_files) == 0:
                assert sample.assignments_dir is not None, \
                    "Received no assignments_dir or analysis files in sample {}".format(sample.name)
                full = False
                assignment_files = list_dir(sample.assignments_dir, ext="assignments.tsv")
                # get assignments
                if len(assignment_files) == 0:
                    print("[CreateHdpTrainingData] filtering 'full' output files")
                    assignment_files = list_dir(sample.assignments_dir, ext="ard.tsv")
                    full = True
            else:
                full = False

                if sample.assignments_dir is not None:
                    print("[CreateHdpTrainingData] WARNING: Using sample analysis files when "
                          "assignments_dir is also set: {}".format(sample.name))
                assignment_files = [x for x in sample.analysis_files if x.endswith("assignments.tsv")]
            # TODO need to get this to infer type before reading in
            # infer kmer length and get kmers
            sample_assignment_table = get_assignment_table(assignment_files[0], 0, full)
            self.set_kmer_len(len(sample_assignment_table.iloc[0]['kmer']))
            # kmers = self.get_sample_kmers(sample)
            base_dir = os.path.dirname(self.out_file_path)
            temp_dir = os.path.join(base_dir, "top_n_kmers")
            os.mkdir(temp_dir)
            generate_top_n_kmers_from_sa_output(assignment_files, temp_dir, self.out_file_path,
                                                sample.number_of_kmer_assignments, alphabet=self.alphabet,
                                                kmer_len=self.k, min_prob=sample.probability_threshold,
                                                worker_count=self.jobs, random=False,
                                                complement=True,
                                                remove=False, alignment=full)
            shutil.rmtree(temp_dir)
            # final_output.append(
            #     multiprocess_make_kmer_assignment_tables(assignment_files, kmers, self.strands,
            #                                              min_probability=sample.probability_threshold,
            #                                              verbose=verbose, full=full,
            #                                              max_assignments=sample.number_of_kmer_assignments,
            #                                              worker_count=self.jobs))

        # master_assignment_table = pd.concat(final_output, ignore_index=True)
        # self.n_assignments = len(master_assignment_table)
        # master_assignment_table.to_csv(self.out_file_path, sep='\t', header=False, index=False)
        return self.out_file_path

    # def write_hdp_training_file2(self, verbose=False):
    #     """Write a hdp training file to a specified location"""
    #     final_output = []
    #     for sample in self.samples:
    #         print("Filtering and gathering {} assignment.tsv files".format(sample.name))
    #         if len(sample.analysis_files) == 0:
    #             assert sample.assignments_dir is not None, \
    #                 "Received no assignments_dir or analysis files in sample {}".format(sample.name)
    #             assignment_files = list_dir(sample.assignments_dir, ext="assignments.tsv")
    #             if len(assignment_files) > 0:
    #                 sample_assignment_table = \
    #                     multiprocess_make_master_assignment_table(assignment_files,
    #                                                               min_probability=sample.probability_threshold,
    #                                                               worker_count=self.jobs)
    #             else:
    #                 print("[CreateHdpTrainingData] filtering 'full' output files")
    #                 assignment_files = list_dir(sample.assignments_dir, ext="ard.tsv")
    #                 sample_assignment_table = \
    #                     multiprocess_make_master_assignment_table(assignment_files,
    #                                                               min_probability=sample.probability_threshold,
    #                                                               full=True, worker_count=self.jobs)
    #         else:
    #             if sample.assignments_dir is not None:
    #                 print("[CreateHdpTrainingData] WARNING: Using sample analysis files when "
    #                       "assignments_dir is also set: {}".format(sample.name))
    #             assignment_files = [x for x in sample.analysis_files if x.endswith("assignments.tsv")]
    #             sample_assignment_table = \
    #                 multiprocess_make_master_assignment_table(assignment_files,
    #                                                           min_probability=sample.probability_threshold,
    #                                                           worker_count=self.jobs)
    #
    #         self.set_kmer_len(len(sample_assignment_table.iloc[0]['kmer']))
    #         # get kmers associated with each sample
    #         kmers = self.get_sample_kmers(sample)
    #         # write correctly formated output
    #         final_output.append(generate_buildAlignments(sample_assignment_table, kmers,
    #                                                      max_assignments=sample.number_of_kmer_assignments,
    #                                                      strands=self.strands, verbose=verbose))
    #     master_assignment_table = pd.concat(final_output, ignore_index=True)
    #     self.n_assignments = len(master_assignment_table)
    #     master_assignment_table.to_csv(self.out_file_path, sep='\t', header=False, index=False)
    #     return self.out_file_path

    def set_kmer_len(self, k):
        """Set the kmer length of the samples"""
        if self.k == 0:
            self.k = k
        else:
            assert self.k == k, "self.k was already set and does not match new k. self.k {} != k {}".format(self.k, k)

    def get_sample_kmers(self, sample, kmer_len=None):
        """Get all kmers from a sample, either from the reference sequence, all possible kmers from an alphabet or
            all kmers that cover a modified nucelotide

        :param sample: AbstractSamples object
        :return: set of desired kmers
        """
        reference = None
        if self.k == 0:
            assert kmer_len is not None, "Kmer length was not set. Must set kmer length in order to get sample kmers"
            self.set_kmer_len(kmer_len)
        if sample.kmers_from_reference:
            reference = sample.bwa_reference

        return get_kmers(self.k, alphabet=self.alphabet, motifs=sample.motifs, reference=reference)


def get_hdp_type(requested_type):
    """Get the integer representing the model type for the buildHdpUtils.c program

    :param requested_type: string associated with each type of hdp model
    :return: integer associated with that
    """
    hdp_types = {
        "singleLevelFixed": 0,
        "singleLevelPrior": 1,
        "multisetFixed": 2,
        "multisetPrior": 3,
        "compFixed": 4,
        "compPrior": 5,
        "middleNtsFixed": 6,
        "middleNtsPrior": 7,
        "groupMultisetFixed": 8,
        "groupMultisetPrior": 9,
        "singleLevelPrior2": 10,
        "multisetPrior2": 11,
        "multisetPriorEcoli": 12,
        "singleLevelPriorEcoli": 13,
        "singleLevelFixedCanonical": 14,
        "singleLevelFixedM6A": 15
    }
    assert (requested_type in list(hdp_types.keys())), "Requested HDP type is invalid, got {}".format(requested_type)
    return hdp_types[requested_type]


class TrainSignalAlign(object):
    """A single class which takes in the only config file used for training and allows for users to train
    either the transitions or emissions of the HMM model
    """
    # global hdp types for specific alphabets and 1D read options
    HDP_TYPES_ACEGOT = [
        ("singleLevelFixed", 0),
        ("singleLevelPrior", 1),
        ("multisetFixed", 2),
        ("multisetPrior", 3),
        ("compFixed", 4),
        ("compPrior", 5),
        ("middleNtsFixed", 6),
        ("middleNtsPrior", 7),
        ("groupMultisetFixed", 8),
        ("groupMultisetPrior", 9),
    ]

    HDP_TYPES_1D = [
        ("singleLevelPrior2", 10),
        ("multisetPrior2", 11),
        ("singleLevelFixedCanonical", 14),
        ("singleLevelFixedM6A", 15),
    ]

    HDP_TYPES_ACEGT = [
        ("singleLevelPrior2", 10),
        ("multisetPrior2", 11),
    ]

    HDP_TYPES_ACGT = [
        ("singleLevelFixedCanonical", 14)
    ]

    HDP_TYPES_ACEGIT = [
        ("multisetPriorEcoli", 12),
        ("singleLevelPriorEcoli", 13),
    ]
    HDP_TYPES_ACFGT = [
        ("singleLevelFixedM6A", 15),
    ]

    def __init__(self, args):
        # TODO Need to create docs here
        """Initialize all objects the training routine may need"""
        # executable
        self.buildHdpUtil = None
        # HDP type
        self.int_hdp_type = None
        # load json and create dot dictionary of all the parameters
        self.args = args
        # check output directory
        self.args.output_dir = os.path.abspath(self.args.output_dir)
        assert os.path.exists(self.args.output_dir), "Output directory does not exist. " \
                                                     "output_dir: {}".format(self.args.output_dir)
        self.working_folder = FolderHandler()
        self.working_path = self.working_folder.open_folder(os.path.join(self.args.output_dir, "tempFiles_trainModels"))

        # create samples from self.args.samples
        self.samples = self._create_samples()

        # Current model paths
        self.template_hmm_model_path = self.args.template_hmm_model
        self.template_hdp_model_path = self.args.template_hdp_model
        self.complement_hmm_model_path = self.args.complement_hmm_model
        self.complement_hdp_model_path = self.args.complement_hdp_model

        # Current SignalHmm model objects
        self.complement_model = None
        self.template_model = None

        # globals for experiments
        self.path_to_bin = self.args.path_to_bin
        self.debug = self.args.debug
        self.two_d = self.args.two_d
        self.job_count = self.args.job_count
        # state machine type changes for SignalAlignment so it can expect an HDP or not
        self.state_machine_type = "threeState"
        self.kmer_length = None
        self.alphabet = None
        # check config file
        self._check_config()

    def _create_samples(self):
        """Create SignalAlignSample for each sample"""
        sa_args = [merge_dicts([s,
                                {"quality_threshold": self.args.filter_reads, "workers": self.args.job_count}])
                   for s in self.args.samples]
        return [SignalAlignSample(working_folder=self.working_folder, **s) for s in sa_args]

    def new_working_folder(self, append):
        """Create new working folder in order to keep track of each new run of analysis"""
        self.working_folder = FolderHandler()
        self.working_path = self.working_folder.open_folder(os.path.join(self.args.output_dir,
                                                                         "tempFiles_trainModels_" + str(append)))

    def train_hdp(self, iteration=""):
        """Train hdp.... duh?
        :param outpath: output file path
        :param build_alignment: path to alignment file
        :param num_alignments: number of alignments in alignment file
        :param threshold:
        :param verbose:
        :param path_to_bin
        :param twoD:
        :param hdp_type: Build Hdp, specify type, options: "Prior, Fixed, twoWay. twoWay is a Prior-type model (recommended)"
        # initial HDP
        :param template_model: Input template lookup table
        :param complement_model: Input complement lookup table
        # fixed concentration models
        :param base_gamma:
        :param middle_gamma:
        :param leaf_gamma:
        # gamma prior models
        :param base_alpha:
        :param base_beta:
        :param middle_alpha:
        :param middle_beta:
        :param leaf_alpha:
        :param leaf_beta:
        # gibbs sampling
        :param gibbs_samples: number of gibbs samples
        :param thinning: how many thinning draws?
        # sample grid
        :param grid_start:
        :param grid_end:
        :param grid_length:
        :param kmer_length: length of kmer
        :return: dictionary of hdp training options
        """
        if iteration:
            iteration = iteration + "."

        if self.args.hdp_args.built_alignments:
            assert os.path.isfile(self.args.hdp_args.built_alignments), \
                "Build alignment file does not exist. {}".format(self.args.hdp_args.built_alignments)
            build_alignment_path = self.args.hdp_args.built_alignments
            num_alignments = count_lines_in_file(build_alignment_path)
        else:
            # set strands which will built
            template = True
            complement = False
            if self.two_d:
                complement = True
            # create instance
            hdp_data = CreateHdpTrainingData(self.samples, os.path.join(self.working_path,
                                                                        "buildAlignment" + iteration + ".tsv"),
                                             template=template,
                                             complement=complement,
                                             verbose=self.debug,
                                             alphabet=self.template_model.alphabet,
                                             jobs=self.job_count)
            # write an hdp training file to path
            build_alignment_path = hdp_data.write_hdp_training_file(verbose=True)
            num_alignments = hdp_data.n_assignments

        verbose_flag = "--verbose "
        # create the output paths for the models
        template_hdp_location = os.path.join(self.working_path,
                                             "template." + iteration + self.args.hdp_args.hdp_type + ".nhdp")
        complement_hdp_location = None
        if self.two_d:
            one_d = None
            complement_hdp_location = os.path.join(self.working_path,
                                                   "complement." + iteration + self.args.hdp_args.hdp_type + ".nhdp")
        else:
            one_d = '--oneD'

        # if we're making a HDP with fixed concentration parameters
        build_initial_hdp_command = "{buildHdpUtil} {verbose}-p {hdpType} -v {tHdpLoc} -w {cHdpLoc} -l {buildAln} " \
                                    "-a {kmerLength} -n {gibbs_samples} -I {burnIn} -t {thin} -s {start} -e {end} " \
                                    "-k {len} {oneD} -C {cL} -T {tL} " \
                                    "-g {Ba} -r {Bb} -j {Ma} -y {Mb} -i {La} -u {Lb} -B {base} -M {middle} -L {leaf} " \
                                    "".format(buildHdpUtil=self.buildHdpUtil,
                                              hdpType=self.int_hdp_type,
                                              tHdpLoc=template_hdp_location,
                                              cHdpLoc=complement_hdp_location,
                                              buildAln=build_alignment_path,
                                              gibbs_samples=self.args.hdp_args.gibbs_samples,
                                              burnIn=min(30000000,
                                                         int(self.args.hdp_args.burnin_multiplier * num_alignments)),
                                              thin=self.args.hdp_args.thinning,
                                              start=self.args.hdp_args.grid_start,
                                              end=self.args.hdp_args.grid_end,
                                              len=self.args.hdp_args.grid_length,
                                              verbose=verbose_flag,
                                              tL=self.template_hmm_model_path,
                                              cL=self.complement_hmm_model_path,
                                              kmerLength=self.kmer_length,
                                              oneD=one_d,
                                              Ba=self.args.hdp_args.base_alpha,
                                              Bb=self.args.hdp_args.base_beta,
                                              Ma=self.args.hdp_args.middle_alpha,
                                              Mb=self.args.hdp_args.middle_beta,
                                              La=self.args.hdp_args.leaf_alpha,
                                              Lb=self.args.hdp_args.leaf_beta,
                                              base=self.args.hdp_args.base_gamma,
                                              middle=self.args.hdp_args.middle_gamma,
                                              leaf=self.args.hdp_args.leaf_gamma)

        print("[[trainModels_buildHdpUtil] Command: {}\n".format(build_initial_hdp_command))
        check_call(build_initial_hdp_command.split())

        print("[trainModels_buildHdpUtil] - finished training HDP emissions routine")

        # check if the HDP created models
        assert os.path.exists(template_hdp_location), "HDP training did not create template hdp model. {}".format(
            template_hdp_location)
        if complement_hdp_location:
            assert os.path.exists(
                complement_hdp_location), "HDP training did not create complement hdp model. {}".format(
                complement_hdp_location)
        # set class parameters
        self.template_hdp_model_path = template_hdp_location
        self.complement_hdp_model_path = complement_hdp_location
        self.state_machine_type = "threeStateHdp"
        return self.template_hdp_model_path, self.complement_hdp_model_path

    def train_normal_hmm(self, transitions=True, emissions=False, iteration=""):
        """Train model transitions"""
        i = 0
        if iteration:
            iteration = "_" + iteration
        # start iterating
        while i < self.args.transitions_args.iterations:
            # align all the samples
            self.run_signal_align(get_expectations=True, trim=self.args.transitions_args.training_bases)
            all_sample_files = merge_lists([sample.analysis_files for sample in self.samples])
            assert len(all_sample_files) > 0, "Something failed in multithread signal alignment. We got no sample files"
            # load then normalize the expectations
            template_expectations_files = [x for x in all_sample_files
                                           if x.endswith(".template.expectations.tsv")]

            if len(template_expectations_files) > 0:
                new_template_hmm = self.working_folder.add_file_path("template_trained_{}{}.hmm".format(i, iteration))
                copyfile(self.template_hmm_model_path, new_template_hmm)
                self.template_hmm_model_path = new_template_hmm
                self.template_model.add_and_normalize_expectations(files=template_expectations_files,
                                                                   hmm_file=self.template_hmm_model_path,
                                                                   update_transitions=transitions,
                                                                   update_emissions=emissions)
            if self.two_d:
                new_complement_hmm = self.working_folder.add_file_path(
                    "complement_trained_{}{}.hmm".format(i, iteration))
                copyfile(self.complement_hmm_model_path, new_complement_hmm)
                self.complement_hmm_model_path = new_complement_hmm

                complement_expectations_files = [x for x in all_sample_files
                                                 if x.endswith(".complement.expectations.tsv")]
                if len(complement_expectations_files) > 0:
                    self.complement_model.add_and_normalize_expectations(files=complement_expectations_files,
                                                                         hmm_file=self.complement_hmm_model_path,
                                                                         update_transitions=transitions,
                                                                         update_emissions=emissions)

            # log the running likelihood
            if len(self.template_model.running_likelihoods) > 0 and \
                    (self.two_d and len(self.complement_model.running_likelihoods)) > 0:
                print("[trainModels_transitions] {i}| {t_likelihood}\t{c_likelihood}".format(
                    t_likelihood=self.template_model.running_likelihoods[-1],
                    c_likelihood=self.complement_model.running_likelihoods[-1],
                    i=i))
                if self.args.transitions_args.test and (len(self.template_model.running_likelihoods) >= 2) and \
                        (self.two_d and len(self.complement_model.running_likelihoods) >= 2):
                    assert (self.template_model.running_likelihoods[-2] < self.template_model.running_likelihoods[
                        -1]) and \
                           (self.complement_model.running_likelihoods[-2] < self.complement_model.running_likelihoods[
                               -1]), "Testing: Likelihood error, went up"
            elif len(self.template_model.running_likelihoods) > 0:
                print(
                    "[trainModels_transitions] {i}| {t_likelihood}".format(
                        t_likelihood=self.template_model.running_likelihoods[-1],
                        i=i))
                if self.args.transitions_args.test and (len(self.template_model.running_likelihoods) >= 2):
                    assert (self.template_model.running_likelihoods[-2] <
                            self.template_model.running_likelihoods[-1]), "Testing: Likelihood error, went up"

            i += 1

        print("[trainModels_transitions] - finished training transitions routine")
        return self.template_hmm_model_path, self.complement_hmm_model_path

    def expectation_maximization_training(self):
        """Complete the entire pipeline of training a new HMM-HDP model

        Note: If expectation_maximization is set to true, both the transitions and hdp/hmm_emissions will be trained
        """
        start = timer()

        if self.args.training.normal_emissions:
            print("[trainModels] Training HMM emission distributions is not currently available.")
        if self.args.training.expectation_maximization:
            for i in range(1, self.args.training.em_iterations + 1):
                print("[trainModels] Training HMM transition distributions. iteration: {}".format(i))
                # first train the model transitions
                self.train_normal_hmm(iteration=str(i))
                print("[trainModels] Running Assignment with new HMM transition distributions. "
                      "iteration: {}".format(i))
                # next get assignments
                self.run_signal_align()
                print("[trainModels] Training HDP emission distributions. iteration: {}".format(i))
                # make new hdp
                self.train_hdp(iteration=str(i))
                print([sample.analysis_files for sample in self.samples])
                print(self.template_hdp_model_path)
                print(self.template_hmm_model_path)
                print(self.complement_hmm_model_path)
                print(self.complement_hdp_model_path)
                # self.new_working_folder(append=str(i))
        elif self.args.training.transitions or self.args.training.hdp_emissions or self.args.training.normal_emissions:
            if self.args.training.transitions or self.args.training.normal_emissions:
                print("[trainModels] Training HMM transition distributions.")
                # self.train_transitions()
                self.train_normal_hmm(iteration="",
                                      transitions=self.args.training.transitions,
                                      emissions=self.args.training.normal_emissions)
            if self.args.training.hdp_emissions:
                print("[trainModels] Training HDP emission distributions.")
                if not self.args.hdp_args.built_alignments:
                    self.run_signal_align(check_samples=True)
                self.train_hdp(iteration="")
        else:
            raise AssertionError("Must set one of the following to True. "
                                 "training.transitions: {}, training.hdp_emissions: {}, "
                                 "training.expectation_maximization: "
                                 "{}, training.normal_emissions: {}".format(self.args.training.transitions,
                                              self.args.training.hdp_emissions,
                                              self.args.training.expectation_maximization,
                                              self.args.training.normal_emissions))

        stop = timer()
        print("[trainModels] Complete")
        print("Training Time = {} seconds".format(stop - start))
        print(self.template_hmm_model_path, self.complement_hmm_model_path,
              self.template_hdp_model_path, self.complement_hdp_model_path)

        return self.template_hmm_model_path, self.complement_hmm_model_path, \
               self.template_hdp_model_path, self.complement_hdp_model_path

    def load_hmm_models(self):
        """Load in the correct models depending on what is going to be trained. """
        # load template model
        assert self.template_hmm_model_path, "Missing template model %s" % (self.template_hmm_model_path)
        self.template_hmm_model_path = os.path.abspath(self.template_hmm_model_path)
        self.template_model = HmmModel(self.template_hmm_model_path)
        new_template_hmm = self.working_folder.add_file_path("template_trained.hmm")
        copyfile(self.template_hmm_model_path, new_template_hmm)
        assert os.path.exists(new_template_hmm), "Problem copying default model to {}".format(new_template_hmm)
        self.template_hmm_model_path = new_template_hmm
        # set alphabet and kmer_length
        self.kmer_length = self.template_model.kmer_length
        self.alphabet = self.template_model.alphabet
        # load complement model if 2D
        if self.two_d:
            assert self.complement_hmm_model_path, "Missing complement model: {}".format(self.complement_hmm_model_path)
            self.complement_hmm_model_path = os.path.abspath(self.complement_hmm_model_path)
            self.complement_model = HmmModel(self.complement_hmm_model_path)
            new_complement_hmm = self.working_folder.add_file_path("complement_trained.hmm")
            copyfile(self.complement_hmm_model_path, new_complement_hmm)
            assert os.path.exists(new_complement_hmm), "Problem copying default model to {}".format(new_complement_hmm)
            self.complement_hmm_model_path = new_complement_hmm
            # make sure models match
            assert self.complement_model.kmer_length == self.template_model.kmer_length, \
                "Template model and complement model kmer lengths do not match." \
                " template: {} != complement: {}".format(self.complement_model.kmer_length,
                                                         self.template_model.kmer_length)
            assert self.complement_model.alphabet == self.template_model.alphabet, \
                "Template model and complement model alphabets do not match." \
                " template: {} != complement: {}".format(self.complement_model.alphabet,
                                                         self.template_model.alphabet)
        # get the input HDP models, if they can be found
        if self.template_hdp_model_path:
            self.state_machine_type = "threeStateHdp"
            assert os.path.exists(self.template_hdp_model_path), \
                "Template HDP path not found {}".format(self.template_hdp_model_path)
            self.template_hdp_model_path = os.path.abspath(self.template_hdp_model_path)
            new_template_hdp = self.working_folder.add_file_path(
                "{}".format(os.path.basename(self.template_hdp_model_path)))
            copyfile(self.template_hdp_model_path, new_template_hdp)
            self.complement_hdp_model_path = new_template_hdp
        # same for complement hdp
        if self.complement_hdp_model_path and self.two_d:
            assert os.path.exists(self.complement_hdp_model_path), \
                "Complement HDP path not found {}".format(self.complement_hdp_model_path)
            self.complement_hdp_model_path = os.path.abspath(self.complement_hdp_model_path)
            new_complement_hdp = \
                self.working_folder.add_file_path("{}".format(os.path.basename(self.complement_hdp_model_path)))
            copyfile(self.complement_hdp_model_path, new_complement_hdp)
            self.complement_hdp_model_path = new_complement_hdp

    def _check_train_transitions_config(self):
        assert isinstance(self.args.transitions_args.iterations, int), \
            "args.transitions_args.iterations must be an integer. {}".format(self.args.transitions_args.iterations)
        assert isinstance(self.args.job_count, int), \
            "args.job_count must be an integer. {}".format(self.args.job_count)

    def _check_train_hdp_config(self):
        """Check if the input parameters will for training the HDP."""
        # make sure hdp type works with alphabet and 1D
        self.int_hdp_type = get_hdp_type(self.args.hdp_args.hdp_type)
        if not self.args.two_d:
            assert (self.args.hdp_args.hdp_type, self.int_hdp_type) in set(self.HDP_TYPES_1D), \
                "HDP type is not compatible with 1D. {}: 1D types {}".format(self.args.hdp_type,
                                                                             self.HDP_TYPES_1D)
        if self.alphabet == "ACEGOT":
            assert (self.args.hdp_args.hdp_type, self.int_hdp_type) in set(self.HDP_TYPES_ACEGOT), \
                "HDP type is not compatible with alphabet=ACEGOT." \
                "Hdp_type: {}, ACEGOT HDP types:  {}".format(self.args.hdp_type, self.HDP_TYPES_ACEGOT)

        elif self.alphabet == "ACEGIT":
            assert (self.args.hdp_args.hdp_type, self.int_hdp_type) in set(self.HDP_TYPES_ACEGIT), \
                "HDP type is not compatible with alphabet=ACEGIT." \
                "Hdp_type: {}, ACEGIT HDP types:  {}".format(self.args.hdp_type, self.HDP_TYPES_ACEGIT)

        elif self.alphabet == "ACEGT":
            assert (self.args.hdp_args.hdp_type, self.int_hdp_type) in set(self.HDP_TYPES_ACEGT), \
                "HDP type is not compatible with alphabet=ACEGT." \
                "Hdp_type: {}, ACEGT HDP types:  {}".format(self.args.hdp_type, self.HDP_TYPES_ACEGT)

        elif self.alphabet == "ACGT":
            assert (self.args.hdp_args.hdp_type, self.int_hdp_type) in set(self.HDP_TYPES_ACGT), \
                "HDP type is not compatible with alphabet=ACGT." \
                "Hdp_type: {}, ACGT HDP types:  {}".format(self.args.hdp_type, self.HDP_TYPES_ACGT)
        elif self.alphabet == "ACFGT":
            assert (self.args.hdp_args.hdp_type, self.int_hdp_type) in set(self.HDP_TYPES_ACFGT), \
                "HDP type is not compatible with alphabet=ACFGT." \
                "Hdp_type: {}, ACFGT HDP types:  {}".format(self.args.hdp_type, self.HDP_TYPES_ACFGT)
        else:
            raise AssertionError("Cannot create a HDP with proved alphabet")

        # check buildHdpUtil executable
        self.buildHdpUtil = os.path.join(self.args.path_to_bin, "./buildHdpUtil")
        assert (os.path.exists(self.buildHdpUtil)), "ERROR: Didn't find buildHdpUtil. {}".format(self.buildHdpUtil)
        # check other parameter inconsistencies
        if self.args.hdp_args.built_alignments:
            assert self.args.training.expectation_maximization is not True, "Cannot use 'built_alignments' file for " \
                                                                            "EM training. Either set " \
                                                                            "training.expectation_maximization to " \
                                                                            "false or change " \
                                                                            "hdp_args.built_alignments to null"
            assert os.path.isfile(self.args.hdp_args.built_alignments), \
                "Build alignment file does not exist. {}".format(self.args.hdp_args.built_alignments)

    def _check_config(self):
        """Make sure training configuration file is correctly filled out"""
        # check model files and load HMM models into memory for training transitions
        self.load_hmm_models()
        # check path to bin
        assert os.path.isdir(self.path_to_bin), "path_to_bin does not exist. " \
                                                "path_to_bin: {}".format(self.path_to_bin)

        # check if signalMachine is found
        assert os.path.exists(os.path.join(self.args.path_to_bin, "./signalMachine")), \
            "ERROR: Didn't find signalMachine executable. {}".format(os.path.join(self.args.path_to_bin,
                                                                                  "./signalMachine"))

        if self.args.training.transitions or self.args.training.expectation_maximization:
            self._check_train_transitions_config()

        if self.args.training.hdp_emissions or self.args.training.expectation_maximization:
            self._check_train_hdp_config()

        return self.args

    def run_signal_align(self, output_format="assignments", get_expectations=False, trim=False, check_samples=True):
        """Run signal align with specified arguments"""
        alignment_args = create_signalAlignment_args(
            destination=self.working_path,
            stateMachineType=self.state_machine_type,
            in_templateHmm=self.template_hmm_model_path,
            in_complementHmm=self.complement_hmm_model_path,
            in_templateHdp=self.template_hdp_model_path,
            in_complementHdp=self.complement_hdp_model_path,
            diagonal_expansion=self.args.diagonal_expansion,
            constraint_trim=self.args.constraint_trim,
            traceBackDiagonals=self.args.traceBackDiagonals,
            twoD_chemistry=self.two_d,
            get_expectations=get_expectations,
            path_to_bin=self.path_to_bin,
            check_for_temp_file_existance=True,
            threshold=self.args.signal_alignment_args.threshold,
            track_memory_usage=self.args.signal_alignment_args.track_memory_usage,
            embed=self.args.signal_alignment_args.embed,
            event_table=self.args.signal_alignment_args.event_table,
            output_format=output_format,
            filter_reads=self.args.filter_reads,
            delete_tmp=self.args.signal_alignment_args.delete_tmp,
            rna=self.args.rna)

        dont_run_sa_samples = []
        run_sa_samples = []
        if check_samples:
            for sample in self.samples:
                if sample.assignments_dir is not None:
                    dont_run_sa_samples.append(sample)
                else:
                    run_sa_samples.append(sample)
            if len(run_sa_samples) > 0:
                run_sa_samples = multithread_signal_alignment_samples(run_sa_samples, alignment_args, self.job_count,
                                                                      trim=trim, debug=self.debug)
            self.samples = merge_lists([run_sa_samples, dont_run_sa_samples])
        else:
            self.samples = multithread_signal_alignment_samples(self.samples, alignment_args, self.job_count,
                                                                trim=trim, debug=self.debug)
        return self.samples


def make_alignment_line(strand, kmer, prob, event):
    """Convert strand, kmer, probability and event to a correctly formatted alignment file line
    :param strand: 't' or 'c' representing template or complement strand of read
    :param kmer: nucleotide kmer
    :param prob: probability of kmer coming from certain event
    :param event: mean of the corresponding event
    :return: correctly formatted alignment line
    """
    assert strand in ['c', 't'], "Strand must be either 'c' or 't'. strand: {}".format(strand)
    assert isinstance(prob, float), "probability must be a float: prob {}".format(prob)
    assert isinstance(kmer, str), "kmer must be a string: kmer {}".format(kmer)
    assert isinstance(event, float), "event must be a float: event {}".format(event)
    entry_line = "blank\t0\tblank\tblank\t{strand}\t0\t0.0\t0.0\t0.0\t{kmer}\t0.0\t0.0\t{prob}\t{event}\t0.0\n"
    return entry_line.format(strand=strand, kmer=kmer, prob=prob, event=event)


def main():
    def parse_args():
        parser = ArgumentParser()
        subparsers = parser.add_subparsers(dest="command")

        # parsers for running the full pipeline
        run_parser = subparsers.add_parser("run", help="runs full workflow ")
        run_parser.add_argument('--config', '-c', type=str,
                                help='Path to the (filled in) config file, generated with "generate".')
        return parser.parse_args()

    args = parse_args()
    if args.command == "run":
        if not os.path.exists(args.config):
            print("{config} not found".format(config=args.config))
            exit(1)
        # run training
        config_args = create_dot_dict(load_json(args.config))
        copyfile(args.config, os.path.join(config_args.output_dir, os.path.basename(args.config)))
        TrainSignalAlign(config_args).expectation_maximization_training()
    else:
        print("Error, try: `trainModels run --config path/to/config.json`")


if __name__ == "__main__":
    sys.exit(main())
