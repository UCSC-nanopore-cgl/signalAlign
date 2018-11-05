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

from py3helpers.utils import create_dot_dict, merge_lists, all_string_permutations, save_json, load_json, \
    count_lines_in_file

from signalalign.signalAlignment import multithread_signal_alignment_samples, create_signalAlignment_args, \
    SignalAlignSample
from signalalign.hiddenMarkovModel import HmmModel
from signalalign.utils.fileHandlers import FolderHandler
from signalalign.utils.parsers import read_fasta
from signalalign.utils.sequenceTools import get_motif_kmers, get_sequence_kmers


def parse_assignment_file(file_path):
    """Parse the .assignments.tsv output file from signalAlign:

    :param file_path: path to assignments file
    :return: panda DataFrame with column names "kmer", "strand", "level_mean", "prob"
    """
    data = pd.read_table(file_path,
                         usecols=(0, 1, 2, 3),
                         names=["kmer", "strand", "level_mean", "prob"],
                         dtype={"kmer": np.str, "strand": np.str, "level_mean": np.float64, "prob": np.float64},
                         header=None
                         )
    return data


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


def make_master_assignment_table(list_of_assignment_paths):
    """Create a master assignment table from a list of assignment paths

    :param list_of_assignment_paths: list of all paths to assignment.tsv files to concat
    :return: pandas DataFrame of all assignments
    """
    assignment_dfs = []
    for f in list_of_assignment_paths:
        assignment_dfs.append(parse_assignment_file(f))
    return pd.concat(assignment_dfs)


class CreateHdpTrainingData(object):
    """Process the assignment files created from SignalAlign for the HDP distribution estimation"""

    def __init__(self, samples, out_file_path, template=True, complement=False, verbose=True):
        """
        Control how each kmer/event assignment is processed given a set of samples and the parameters associated with
        each sample

        :param samples:
        :param out_file:
        :param template: generate kmers for template read strand: default: True
        :param complement: generate kmers for complement read strand: default: True
        :param min_probability: the minimum probability to use for assigning kmers
        :param verbose: option to print update statements
        """
        self.strands = []
        if template:
            self.strands.append('t')
        if complement:
            self.strands.append('c')
        assert self.strands != [], 'template or complement need to be set to True. ' \
                                   'complement: {}, template: {}'.format(complement, template)

        for sample in samples:
            assert isinstance(sample, SignalAlignSample)

        self.canonical = "ATGC"
        self.samples = samples
        self.out_file_path = out_file_path
        self.template = template
        self.complement = complement
        self.verbose = verbose
        self.master_assignment_table = \
            make_master_assignment_table(sorted(merge_lists([sample.analysis_files for sample in self.samples])))
        self.k = len(self.master_assignment_table.iloc[0]['kmer'])
        self.n_assignments = len(self.master_assignment_table)

    def generate_hdp_training_lines_wrapper(self, kmer_list, max_assignments=100, min_probability=0.8):
        """Convert assignments to alignment line format for HDP training.

        Filter assignments on a minimum probability, read strand, and a max number of kmer assignments

        :param kmer_list: list of kmers to write to alignment file
        :param max_assignments: max number of assignments to process for each kmer
        :param min_probability: the minimum probability to use for assigning kmers
        """
        # loop through for each strand in the assignments
        return self._generate_hdp_training_lines(self.master_assignment_table, kmer_list,
                                                 max_assignments=max_assignments,
                                                 strands=self.strands, min_probability=min_probability,
                                                 verbose=self.verbose)

    def write_hdp_training_file(self):
        """Write a hdp training file to a specified location"""
        with open(self.out_file_path, 'w') as out_file:
            for sample in self.samples:
                # get kmers associated with each sample
                kmers = self.get_sample_kmers(sample)
                # write correctly formated output
                for line in self.generate_hdp_training_lines_wrapper(kmers,
                                                                     max_assignments=sample.number_of_kmer_assignments,
                                                                     min_probability=sample.probability_threshold):
                    out_file.write(line)
        return self.out_file_path

    def get_sample_kmers(self, sample):
        """Get all kmers from a sample, either from the reference sequence, all possible kmers from an alphabet or
            all kmers that cover a modified nucelotide

        :param sample: AbstractSamples object
        :return: set of desired kmers
        """
        kmers = set()
        # if motifs is present, process for all motifs with modified base
        if sample.motifs:
            for motif in sample.motifs:
                kmers |= get_motif_kmers(motif, self.k, alphabet=self.canonical)
        # if we want to limit kmers which were seen in reference sequence
        if sample.kmers_from_reference:
            for _, _, sequence in read_fasta(sample.bwa_reference):
                kmers |= get_sequence_kmers(sequence, k=self.k, rev_comp=True)
        else:
            kmers |= {x for x in all_string_permutations(self.canonical, length=self.k)}

        return kmers

    @staticmethod
    def _generate_hdp_training_lines(assignments, kmer_list, max_assignments=10,
                                     strands=('t', 'c'), min_probability=0.8, verbose=False):
        """Convert assignments to alignment line format for HDP training.

        Filter assignments on a minimum probability, read strand, and a max number of kmer assignments

        :param assignments: pandas array of assignments to search
        :param template: generate kmers for template read strand: default: True
        :param complement: generate kmers for complement read strand: default: True
        :param verbose: option to print update statements
        :param kmer_list: list of kmers to write to alignment file
        :param max_assignments: max number of assignments to process for each kmer
        :param min_probability: the minimum probability to use for assigning kmers
        """
        # loop through for each strand in the assignments
        assert isinstance(strands, list) and len(strands) > 0, "strands must be a list and not be empty. strands: {}" \
                                                               "".format(strands)
        for strand in strands:
            by_strand = assignments.loc[(assignments['strand'] == strand)
                                        & (assignments['prob'] >= min_probability)]

            for k in kmer_list:
                kmer_assignments = by_strand.loc[by_strand['kmer'] == k]
                if kmer_assignments.empty and verbose:
                    print("missing kmer {}, continuing".format(k))
                    continue
                kmer_assignments = kmer_assignments.sort_values(['prob'], ascending=0)
                n = 0
                for _, r in kmer_assignments.iterrows():
                    yield make_alignment_line(strand=r['strand'], kmer=r['kmer'], event=r['level_mean'], prob=r['prob'])
                    n += 1
                    if n >= max_assignments:
                        break
                if n < max_assignments and verbose:
                    print("WARNING didn't find {max} requested assignments for {kmer} only found {found}"
                          "".format(max=max_assignments, kmer=k, found=n))


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
        "singleLevelFixedCanonical": 14
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
        ("singleLevelFixedCanonical", 14)
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
        return [SignalAlignSample(working_folder=self.working_folder, **s) for s in self.args.samples]

    def new_working_folder(self, append):
        """Create new working folder in order to keep track of each new run of analysis"""
        self.working_folder = FolderHandler()
        self.working_path = self.working_folder.open_folder(os.path.join(self.args.output_dir,
                                                                         "tempFiles_trainModels_" + str(append)))

    def train_hdp(self):
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
            hdp_data = CreateHdpTrainingData(self.samples, os.path.join(self.working_path, "buildAlignment.tsv"),
                                             template=template,
                                             complement=complement,
                                             verbose=self.debug)
            # write an hdp training file to path
            build_alignment_path = hdp_data.write_hdp_training_file()
            num_alignments = hdp_data.n_assignments

        verbose_flag = "--verbose " if self.debug is True else ""
        # create the output paths for the models
        template_hdp_location = os.path.join(self.working_path, "template." + self.args.hdp_args.hdp_type + ".nhdp")
        complement_hdp_location = None
        if self.two_d:
            one_d = None
            complement_hdp_location = os.path.join(self.working_path,
                                                   "complement." + self.args.hdp_args.hdp_type + ".nhdp")
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
                                              burnIn=int(self.args.hdp_args.burnin_multiplier * num_alignments),
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

    def train_normal_hmm(self, transitions=True, emissions=False):
        """Train model transitions"""
        i = 0
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
                self.template_model.add_and_normalize_expectations(files=template_expectations_files,
                                                                   hmm_file=self.template_hmm_model_path,
                                                                   update_transitions=transitions,
                                                                   update_emissions=emissions)
            if self.two_d:
                complement_expectations_files = [x for x in all_sample_files
                                                 if x.endswith(".complement.expectations.tsv")]
                if len(complement_expectations_files) > 0:
                    self.complement_model.add_and_normalize_expectations(files=complement_expectations_files,
                                                                         hmm_file=self.complement_model_path,
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
                self.train_normal_hmm()
                print("[trainModels] Running Assignment with new HMM transition distributions. "
                      "iteration: {}".format(i))
                # next get assignments
                self.run_signal_align()
                print("[trainModels] Training HDP emission distributions. iteration: {}".format(i))
                # make new hdp
                self.train_hdp()
                print([sample.analysis_files for sample in self.samples])
                print(self.template_hdp_model_path)
                print(self.template_hmm_model_path)
                print(self.complement_hmm_model_path)
                print(self.complement_hdp_model_path)
                # self.new_working_folder(append=str(i))
        elif self.args.training.transitions or self.args.training.hdp_emissions:
            if self.args.training.transitions:
                print("[trainModels] Training HMM transition distributions.")
                # self.train_transitions()
                self.train_normal_hmm()
            if self.args.training.hdp_emissions:
                print("[trainModels] Training HDP emission distributions.")
                if not self.args.hdp_args.built_alignments:
                    self.run_signal_align()
                self.train_hdp()
        else:
            raise AssertionError("Must set one of the following to True. "
                                 "training.transitions: {}, training.hdp_emissions: {}, "
                                 "training.expectation_maximization: "
                                 "{}.".format(self.args.training.transitions,
                                              self.args.training.hdp_emissions,
                                              self.args.training.expectation_maximization))

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

    def run_signal_align(self, output_format="assignments", get_expectations=False, trim=False):
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
            delete_tmp=self.args.signal_alignment_args.delete_tmp)

        self.samples = multithread_signal_alignment_samples(self.samples, alignment_args, self.job_count,
                                                            trim=trim, debug=self.debug)
        return self.samples


def main():
    def parse_args():
        parser = ArgumentParser()
        subparsers = parser.add_subparsers(dest="command")

        # parsers for running the full pipeline
        run_parser = subparsers.add_parser("run", help="runs full workflow ")
        run_parser.add_argument('--config', default='trainModels-config.yaml', type=str,
                                help='Path to the (filled in) config file, generated with "generate".')
        return parser.parse_args()

    args = parse_args()
    if args.command == "run":
        if not os.path.exists(args.config):
            print("{config} not found".format(config=args.config))
            exit(1)
        # run training
        config_args = create_dot_dict(load_json(args.config))
        TrainSignalAlign(config_args).expectation_maximization_training()
    else:
        print("Error, try: `trainModels run --config path/to/config.json`")


if __name__ == "__main__":
    sys.exit(main())
