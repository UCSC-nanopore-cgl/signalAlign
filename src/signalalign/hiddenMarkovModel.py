#!/usr/bin/env python
"""hiddenMarkovModel.py contains objects for handling HMMs for SignalAlign"""
########################################################################
# File: hiddenMarkovModel.py
#  executable: hiddenMarkovModel.py
#
# Author: Andrew Bailey
# History: 08/10/18 Created
########################################################################


from __future__ import print_function
import sys
import os
import numpy as np
import pandas as pd

from itertools import product
from scipy.stats import norm, invgauss, entropy
from scipy.spatial.distance import euclidean
from sklearn.neighbors import KernelDensity
from py3helpers.utils import all_string_permutations
from py3helpers.seq_tools import is_non_canonical_iupac_base

import matplotlib as mpl
if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
mpl.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Globals
NORM_DIST_PARAMS = 2
NB_MODEL_PARAMS = 5
_SQRT2 = np.sqrt(2)


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


class HmmModel(object):
    def __init__(self, ont_model_file, hdp_model_file=None, rna=False):
        # TODO Need to create docs here
        assert os.path.exists(ont_model_file)
        self.rna = rna
        self.ont_model_file = ont_model_file
        self.match_model_params = 5  # level_mean, level_sd, noise_mean, noise_sd, noise_lambda
        self.state_number = 3
        self.transitions = np.zeros(self.state_number**2)
        self.transitions_expectations = np.zeros(self.state_number**2)
        self.likelihood = 0.0
        self.running_likelihoods = []
        self.alphabet_size = 0
        self.alphabet = ""
        self.kmer_length = 0
        self.has_ont_model = False
        self.normalized = False
        self.sorted_kmer_tuple = tuple()
        self.num_kmers = 0
        # HDP stuff here
        self.kmer_assignments = []
        self.event_assignments = []
        self.assignments_record = []
        self.symbol_set_size = 0
        # event model for describing normal distributions for each kmer
        self.event_model = {"means": np.zeros(self.symbol_set_size),
                            "SDs": np.zeros(self.symbol_set_size),
                            "noise_means": np.zeros(self.symbol_set_size),
                            "noise_SDs": np.zeros(self.symbol_set_size),
                            "noise_lambdas": np.zeros(self.symbol_set_size)}
        self.set_default_transitions()

        # bins for expectations
        self.load_model(self.ont_model_file)
        self.mean_expectations = np.zeros(self.symbol_set_size)
        self.sd_expectations = np.zeros(self.symbol_set_size)
        self.posteriors = np.zeros(self.symbol_set_size)
        self.observed = np.zeros(self.symbol_set_size, dtype=bool)

        self.has_hdp_model = False
        # Load HDP model if passed
        if hdp_model_file:
            assert os.path.exists(hdp_model_file)
            self.hdp_path = hdp_model_file
            self.splines_finalized = False
            self.has_data = False
            self.sample_gamma = False
            self.num_dps = 0
            self.mu = 0
            self.nu = 0
            self.alpha = 0
            self.beta = 0
            self.grid_start = 0
            self.grid_stop = 0
            self.grid_length = 0
            self.gamma_alpha = 0
            self.gamma_beta = 0
            self.w = 0
            self.s = 0
            self.all_posterior_pred = []
            self.all_spline_slopes = []

            self._initialize_hdp_model()

    def normalize_transitions_expectations(self):
        """Normalize transitions from each state to the other states

        eg: MATCH_CONTINUE = MATCH_CONTINUE / (GAP_OPEN_Y + GAP_OPEN_X + MATCH_CONTINUE)
        """
        for from_state in range(self.state_number):
            i = self.state_number * from_state
            j = sum(self.transitions_expectations[i:i + self.state_number])
            for to_state in range(self.state_number):
                self.transitions_expectations[i + to_state] = self.transitions_expectations[i + to_state] / j

    def set_default_transitions(self):
        MATCH_CONTINUE = np.exp(-0.23552123624314988)     # stride
        GAP_OPEN_X = np.exp(-1.6269694202638481)          # skip
        GAP_OPEN_Y = np.exp(-4.3187242127300092)          # 1 - (skip + stride)

        MATCH_FROM_GAP_X = np.exp(-0.21880828092192281)   # 1 - skip'
        GAP_EXTEND_X = np.exp(-1.6269694202638481)        # skip'
        GAP_SWITCH_TO_Y = 0.0

        GAP_EXTEND_Y = np.exp(-4.3187242127239411)        # stay (1 - (skip + stay))
        MATCH_FROM_GAP_Y = np.exp(-0.013406326748077823)  # 1 - (skip + stay)
        GAP_SWITCH_TO_X = 0.000000001
        self.transitions = [
            MATCH_CONTINUE, GAP_OPEN_X, GAP_OPEN_Y,
            MATCH_FROM_GAP_X, GAP_EXTEND_X, GAP_SWITCH_TO_Y,
            MATCH_FROM_GAP_Y, GAP_SWITCH_TO_X, GAP_EXTEND_Y
        ]
        return

    def check_header_line(self, line, expectations_file):
        """Make sure that the header line of an expectations file matches the model we are training

        :param line: split header line. eg: ['3', '4', "ACGT", '5']
        :param expectations_file: path to expectations file for error reporting
        :return: True if assert statements pass
        """
        assert len(line) == 4, "signalHmm.check_header_line - incorrect header (param line): {}".format(expectations_file)
        assert int(line[0]) == self.state_number, "signalHmm.check_header_line - state number error should be {exp} " \
                                                  "got {obs}".format(exp=self.state_number, obs=line[0])
        assert int(line[1]) == self.alphabet_size, "signalHmm.check_header_line - alphabet size error incorrect " \
                                                   "parameters: {file}, line {line}".format(file=expectations_file,
                                                                                            line=''.join(line))
        assert line[2] == self.alphabet, "signalHmm.check_header_line - incorrect parameters: {file}, line {line}" \
                                         "".format(file=expectations_file, line=''.join(line))
        assert int(line[3]) == self.kmer_length, "signalHmm.check_header_line - incorrect parameters: {file}, " \
                                                 "line {line}".format(file=expectations_file, line=''.join(line))
        return True

    def load_model(self, model_file):
        """Load HMM model from model file

        the model file has the format:
        line 0: stateNumber \t alphabetSize \t alphabet \t kmerLength
        line 1: match->match \t match->gapX \t match->gapY \t
                gapX->match \t gapX->gapX \t gapX->gapY \t
                gapY->match \t gapY->gapX \t gapY->gapY \n
        line 2: [level_mean] [level_sd] [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n

        :param model_file: path to model file
        """
        assert os.path.exists(model_file), "signalHmm.load_model - didn't find model here: {}".format(model_file)

        with open(model_file, 'r') as fH:

            line = fH.readline().split()
            # check for correct header length
            assert len(line) == 4, "signalHmm.load_model - incorrect line length line:{}".format(''.join(line))
            # check stateNumber
            assert int(line[0]) == self.state_number, "signalHmm.load_model - incorrect stateNumber got {got} should be {exp}" \
                                                      "".format(got=int(line[0]), exp=self.state_number)
            # load model parameters
            self.alphabet_size = int(line[1])
            self.alphabet = line[2]
            self.kmer_length = int(line[3])
            self.symbol_set_size = self.alphabet_size**self.kmer_length
            assert self.symbol_set_size > 0, "signalHmm.load_model - Got 0 for symbol_set_size"
            assert self.symbol_set_size <= 6**6, "signalHmm.load_model - Got more than 6^6 for symbol_set_size got {}" \
                                                 "".format(self.symbol_set_size)

            line = list(map(float, fH.readline().split()))
            assert len(line) == len(self.transitions) + 1, "signalHmm.load_model incorrect transitions line"
            self.transitions = line[:-1]
            self.likelihood = line[-1]

            line = list(map(float, fH.readline().split()))
            assert len(line) == self.symbol_set_size * NB_MODEL_PARAMS, \
                "signalHmm.load_model incorrect event model line"
            self.event_model["means"] = line[::NB_MODEL_PARAMS]
            self.event_model["SDs"] = line[1::NB_MODEL_PARAMS]
            self.event_model["noise_means"] = line[2::NB_MODEL_PARAMS]
            self.event_model["noise_SDs"] = line[3::NB_MODEL_PARAMS]
            self.event_model["noise_lambdas"] = line[4::NB_MODEL_PARAMS]

            assert not np.any(self.event_model["means"] == 0.0), "signalHmm.load_model, this model has 0 E_means"
            assert not np.any(self.event_model["SDs"] == 0.0), "signalHmm.load_model, this model has 0 E_means"
            assert not np.any(self.event_model["noise_means"] == 0.0), "signalHmm.load_model, this model has 0 E_noise_means"
            assert not np.any(self.event_model["noise_SDs"] == 0.0), "signalHmm.load_model, this model has 0 E_noise_SDs"
            self._create_kmer_index_map()
            self.has_ont_model = True

    def write(self, out_file):
        """Write out model file to out_file path
        :param out_file: path to write hmm model file
        """
        # the model file has the format:
        # line 0: stateNumber \t alphabetSize \t alphabet \t kmerLength
        # line 1: match->match \t match->gapX \t match->gapY \t
        #         gapX->match \t gapX->gapX \t gapX->gapY \t
        #         gapY->match \t gapY->gapX \t gapY->gapY \n
        # line 2: [level_mean] [level_sd] [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n
        assert self.has_ont_model, "Shouldn't be writing down a Hmm that has no Model"
        assert self.normalized, "Shouldn't be writing down a not normalized HMM"

        with open(out_file, 'w') as f:

            # line 0
            f.write("{stateNumber}\t{alphabetSize}\t{alphabet}\t{kmerLength}\n"
                    "".format(stateNumber=self.state_number, alphabetSize=self.alphabet_size,
                              alphabet=self.alphabet, kmerLength=self.kmer_length))
            # line 1 transitions
            for i in range(self.state_number * self.state_number):
                f.write("{transition}\t".format(transition=str(self.transitions[i])))
            # likelihood
            f.write("{}\n".format(str(self.likelihood)))

            # line 2 Event Model
            for k in range(self.symbol_set_size):
                f.write("{level_mean}\t{level_sd}\t{noise_mean}\t{noise_sd}\t{noise_lambda}\t"
                        "".format(level_mean=self.event_model["means"][k], level_sd=self.event_model["SDs"][k],
                                  noise_mean=self.event_model["noise_means"][k], noise_sd=self.event_model["noise_SDs"][k],
                                  noise_lambda=self.event_model["noise_lambdas"][k]))
            f.write("\n")

    @staticmethod
    def _get_kmer_index(kmer, alphabet, kmer_length, alphabet_size):
        """Get the model index for a given kmer

        ex: get_kmer_index(AAAAA) = 0
        :param kmer: nucleotide sequence
        """
        assert set(kmer).issubset(set(alphabet)) is True, "Nucleotide not found in model alphabet: kmer={}, " \
                                                          "alphabet={}".format(kmer, alphabet)
        assert len(kmer) == kmer_length, "Kmer length does not match model kmer length"

        alphabet_dict = {base: index for index, base in enumerate(sorted(alphabet))}
        kmer_index = 0
        for index, nuc in enumerate(kmer):
            kmer_index += alphabet_dict[nuc]*(alphabet_size**(kmer_length-index-1))
        return kmer_index

    def get_kmer_index(self, kmer):
        """Get the model index for a given kmer

        ex: get_kmer_index(AAAAA) = 0
        :param kmer: nucleotide sequence
        """
        return self._get_kmer_index(kmer, self.alphabet, self.kmer_length, self.alphabet_size)

    def _create_kmer_index_map(self):
        """Create the kmer_to_index_map and index_to_kmer_map"""
        sorted_kmer_list = []
        for i, kmer_list in enumerate(product(self.alphabet, repeat=self.kmer_length)):
            sorted_kmer_list.append(''.join(kmer_list))
        self.sorted_kmer_tuple = tuple(sorted_kmer_list)
        self.num_kmers = len(self.sorted_kmer_tuple)
        return self.sorted_kmer_tuple

    def index_to_kmer(self, index):
        """Get kmer from a given index

        ex: index_to_kmer(0) = "AAAAA"
        :param index: number representing kmer
        """
        assert index < self.num_kmers, \
            "The kmer index is out of bounds given the alphabet and kmer length. {} > {}".format(index, self.num_kmers)
        return self.sorted_kmer_tuple[index]

    def get_event_mean_gaussian_parameters(self, kmer):
        """Get the model's Normal distribution parameters to model the mean of a specific kmer

        :param kmer: kmer that can fit in model
        """
        kmer_index = self.get_kmer_index(kmer)
        normal_mean = self.event_model["means"][kmer_index]
        normal_sd = self.event_model["SDs"][kmer_index]

        return normal_mean, normal_sd

    def get_event_sd_inv_gaussian_parameters(self, kmer):
        """Get the model's inverse gaussian distribution parameters to model the mean of a specific kmer

        :param kmer: kmer that can fit in model
        """
        kmer_index = self.get_kmer_index(kmer)
        inv_gauss_mean = self.event_model["noise_means"][kmer_index]
        inv_gauss_lambda = self.event_model["noise_lambdas"][kmer_index]
        return inv_gauss_mean, inv_gauss_lambda

    def log_event_mean_gaussian_probability_match(self, event_mean, kmer):
        """Get the probability of the event_mean coming from the model's kmer gaussian/normal distribution
        :param event_mean: mean of event
        :param kmer: nucleotide sequence to check
        """
        normal_mean, normal_sd = self.get_event_mean_gaussian_parameters(kmer)
        return norm.logpdf(event_mean, normal_mean, normal_sd)

    def log_event_sd_inv_gaussian_probability_match(self, event_sd, kmer):
        """Get the probability of the event_sd coming from the model's kmer inv-gaussian distribution
        :param event_sd: sd of event
        :param kmer: kmer for model distribution selection
        """
        inv_gauss_mean, inv_gauss_lambda = self.get_event_sd_inv_gaussian_parameters(kmer)
        return invgauss(inv_gauss_mean/inv_gauss_lambda, scale=inv_gauss_lambda).logpdf(event_sd)

    def add_expectations_file(self, expectations_file):
        """Add expectations file to the HMM. This is used for generating expectations of transition probabilities or
        emission probabilities


                expectations files have the format:
        line 0: stateNumber \t alphabetSize \t alphabet \t kmerLength
        line 1: match->match \t match->gapX \t match->gapY \t
                gapX->match \t gapX->gapX \t gapX->gapY \t
                gapY->match \t gapY->gapX \t gapY->gapY \n
        line 2: [level_mean] [level_sd] [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n
        line 3: event expectations [mean] [sd] / kmer \n
        line 4: posteriors 1 per kmer \n
        line 5: observed 1 per kmer \n


        :param expectations_file: path to signalAlign expectations file
        :return: True if expectations file was in correct format
        """
        if not os.path.exists(expectations_file) or os.stat(expectations_file).st_size == 0:
            print("Empty or missing file {}".format(expectations_file), file=sys.stderr)
            return False

        with open(expectations_file, 'r') as fH:
            # line 0
            line = fH.readline().split()
            self.check_header_line(line=line, expectations_file=expectations_file)
            # line 1: transitions, likelihood
            # check if valid
            line = list(map(float, fH.readline().split()))
            assert len(line) == (len(self.transitions) + 1), \
                "HMM.add_expectations_file - problem with file {f} " \
                "transitions line {l}, incorrect length".format(f=expectations_file, l=''.join(line))

            self.likelihood += line[-1]
            self.transitions_expectations = [sum(x) for x in zip(self.transitions_expectations, line[0:-1])]

            # # line 2: event model
            # line = list(map(float, fH.readline().split()))
            # assert len(line) == self.symbol_set_size * NB_MODEL_PARAMS, "HMM.add_expectations_file - problem with " \
            #                                                             "event model in file {ef}".format(ef=expectations_file)
            #
            # # line 3 event expectations [E_mean, E_sd]
            # line = list(map(float, fH.readline().split()))
            # assert len(line) == self.symbol_set_size * NORM_DIST_PARAMS, \
            #     'HMM: check_file - bad file (event expectations): {}'.format(expectations_file)
            #
            # self.event_assignments += line
            # self.mean_expectations = [i + j for i, j in zip(self.mean_expectations, line[::NORM_DIST_PARAMS])]
            # self.sd_expectations = [i + j for i, j in zip(self.sd_expectations, line[1::NORM_DIST_PARAMS])]
            #
            # # line 4, posteriors
            # line = list(map(float, fH.readline().split()))
            # assert len(line) == self.symbol_set_size, "HMM: check_file - bad file (posteriors): {}".format(expectations_file)
            #
            # self.kmer_assignments += line
            #
            # # line 5, probabilities
            # self.posteriors = [sum(x) for x in zip(self.posteriors, line)]
            # line = list(map(bool, fH.readline().split()))
            # assert len(line) == self.symbol_set_size, "HMM: check_file - bad file (observations): {}".format(expectations_file)
            # self.observed = [any(b) for b in zip(self.observed, line)]
            return True

    def normalize(self, update_transitions, update_emissions):
        """Normalize the transitions and emission probabilities

        :param update_transitions: boolean option to update transitions
        :param update_emissions: boolean option to update emissions
        """
        # update
        if update_transitions is True:
            # normalize transitions expectations
            self.normalize_transitions_expectations()
            for i in range(self.state_number**2):
                    self.transitions[i] = self.transitions_expectations[i]

        # calculate the new expected mean and standard deviation for the kmer normal distributions
        if update_emissions:
            # print(self.observed)
            for k in range(self.symbol_set_size):  # TODO implement learning rate
                # print(k)
                if self.observed[k] is True:
                    u_k = self.mean_expectations[k] / self.posteriors[k]
                    o_k = np.sqrt(self.sd_expectations[k] / self.posteriors[k])
                    if u_k > 0:
                        self.event_model["means"][k] = u_k
                        self.event_model["SDs"][k] = o_k
                else:
                    continue
        self.normalized = True

    def reset_assignments(self):
        """Keep track of number of event assignments processed and reset event and kmer assignments"""
        self.assignments_record.append(len(self.event_assignments))
        self.event_assignments = []
        self.kmer_assignments = []

    def add_and_normalize_expectations(self, files, hmm_file, update_transitions=True, update_emissions=False):
        """Add expectations file to HMM model and update transitions. Emissions are currently unable to be updated
        :param files: list of 'expectation' files to add to model
        :param hmm_file: path to HMM file to write new model
        :param update_transitions: boolean option to update transitions
        :param update_emissions: boolean option to update emissions
        """
        if update_emissions is False and update_transitions is False:
            print("[trainModels] NOTICE: Training transitions by default\n", file=sys.stderr)
            update_transitions = True

        # reset model likelihood and keep track of passing and failing files
        self.likelihood = 0
        files_added_successfully = 0
        files_with_problems = 0
        for f in files:
            try:
                success = self.add_expectations_file(f)
                if success:
                    files_added_successfully += 1
                    os.remove(f)

                else:
                    files_with_problems += 1
            except Exception as e:
                files_with_problems += 1
                print("Problem adding expectations file {file} got error {e}".format(file=f, e=e),
                      file=sys.stderr)

        # normalize, write and keep track of running likelihood
        self.normalize(update_transitions=update_transitions, update_emissions=update_emissions)
        self.write(hmm_file)
        self.running_likelihoods.append(self.likelihood)
        self.reset_assignments()
        print("[trainModels] NOTICE: Added {success} expectations files successfully, {problem} files had problems\n"
              "".format(success=files_added_successfully, problem=files_with_problems), file=sys.stderr)

    def _initialize_hdp_model(self):
        """Read in HDP model and make sure parameters match the ONT model"""
        with open(self.hdp_path, 'r') as hdp_fh:
            hdp_alphabet_size = int(hdp_fh.readline())
            assert self.alphabet_size == hdp_alphabet_size, \
                "ONT Alphabet size does not match HDP model ({} != {})".format(self.alphabet_size, hdp_alphabet_size)
            hdp_alphabet = hdp_fh.readline().rstrip()
            assert self.alphabet == hdp_alphabet, \
                "ONT Alphabet size does not match HDP model ({} != {})".format(self.alphabet, hdp_alphabet)
            hdp_kmer_length = int(hdp_fh.readline())
            assert self.kmer_length == hdp_kmer_length, \
                "ONT Kmer length size does not match HDP model ({} != {})".format(self.kmer_length, hdp_kmer_length)

            self.splines_finalized = bool(int(hdp_fh.readline()))
            self.has_data = bool(int(hdp_fh.readline()))
            self.sample_gamma = bool(int(hdp_fh.readline()))
            self.num_dps = int(hdp_fh.readline())
            self.data = [float(x) for x in hdp_fh.readline().split()]
            self.dp_ids = [int(x) for x in hdp_fh.readline().split()]
            unpack_line = hdp_fh.readline().split()
            self.mu = float(unpack_line[0])
            self.nu = float(unpack_line[1])
            self.alpha = float(unpack_line[2])
            self.beta = float(unpack_line[3])
            unpack_line = hdp_fh.readline().split()
            self.grid_start = int(unpack_line[0])
            self.grid_stop = int(unpack_line[1])
            self.grid_length = int(unpack_line[2])
            self.linspace = np.linspace(self.grid_start, self.grid_stop, num=self.grid_length)
            self.gamma_params = [float(x) for x in hdp_fh.readline().split()]
            if self.sample_gamma:
                self.gamma_alpha = [float(x) for x in hdp_fh.readline().split()]
                self.gamma_beta = [float(x) for x in hdp_fh.readline().split()]
                self.w = [float(x) for x in hdp_fh.readline().split()]
                self.s = [bool(int(x)) for x in hdp_fh.readline().split()]

            for i in range(self.num_dps):
                line = hdp_fh.readline().split()
                parent_id = line[0]
                num_factor_children = line[1]
                if parent_id == '-':
                    pass
                    # print(num_factor_children)

            for id in range(self.num_dps):
                post_pred = [float(x) for x in hdp_fh.readline().split()]
                self.all_posterior_pred.append(post_pred)

            for id in range(self.num_dps):
                spline_slopes = [float(x) for x in hdp_fh.readline().split()]
                self.all_spline_slopes.append(spline_slopes)

            line = hdp_fh.readline()
            factor_list = []
            while line:
                items = line.split()
                if int(items[0]) == 0:
                    fctr = "SOMETHING"
                    param_array = items[2].split(';')
                if int(items[0]) == 1:
                    # new_middle_factor
                    fctr = items[2]
                if int(items[0]) == 2:
                    # new_data_pt_factor
                    fctr = items[2]
                factor_list.append(fctr)
                if items[1] != '-':
                    pass
                line = hdp_fh.readline()
            self.has_hdp_model = True

    @staticmethod
    def grid_spline_interp(query_x, x, y, slope, length):
        # if event mean is below start of grid
        if query_x <= x[0]:
            return y[0] - slope[0] * (x[0] - query_x)
        # if event mean is above end grid
        elif query_x >= x[length - 1]:
            n = length - 1
            return y[n] + slope[n] * (query_x - x[n])
        else:
            dx = x[1] - x[0]
            idx_left = int((query_x - x[0]) // dx)
            idx_right = idx_left + 1

            dy = y[idx_right] - y[idx_left]

            a = slope[idx_left] * dx - dy
            b = dy - slope[idx_right] * dx

            t_left = (query_x - x[idx_left]) / dx
            t_right = 1.0 - t_left

            return t_right * y[idx_left] + t_left * y[idx_right] + t_left * t_right * (a * t_right + b * t_left)

    def plot_kmer_distribution(self, kmer, alignment_file=None, alignment_file_data=None, savefig_dir=None):
        """Plot the distribution of a kmer with ONT and/or HDP distributions"""
        assert self.has_ont_model, "Must have ONT model loaded"
        if savefig_dir:
            assert os.path.exists(savefig_dir), "Save figure directory does not exist: {}".format(savefig_dir)
        # keep track of handles and text depending on which models are loaded
        handles1 = []
        legend_text1 = []
        handles2 = []
        legend_text2 = []

        normal_mean, normal_sd = self.get_event_mean_gaussian_parameters(kmer)

        plt.figure(figsize=(12, 8))
        panel1 = plt.axes([0.1, 0.1, .6, .8])
        panel1.set_xlabel('pA')
        panel1.set_ylabel('Density')
        panel1.grid(color='black', linestyle='-', linewidth=1, alpha=0.5)
        panel1.xaxis.set_major_locator(ticker.AutoLocator())
        panel1.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        min_x = normal_mean-(5*normal_sd)
        max_x = normal_mean+(5*normal_sd)
        panel1.set_xlim(min_x, max_x)
        panel1.set_title(label=kmer)
        # plot ont normal distribution
        x = np.linspace(normal_mean - 4*normal_sd, normal_mean + 4*normal_sd, 200)
        ont_handle, = panel1.plot(x, norm.pdf(x, normal_mean, normal_sd))
        # panel1.plot([normal_mean, normal_mean], [0, norm.pdf(normal_mean, normal_mean, normal_sd)], lw=2)
        ont_model_name = os.path.basename(self.ont_model_file)
        txt_handle1, = panel1.plot([], [], ' ')
        txt_handle2, = panel1.plot([], [], ' ')
        txt_handle3, = panel1.plot([], [], ' ')

        handles1.append(ont_handle)
        legend_text1.append("ONT Normal Distribution")

        handles2.extend([txt_handle1, txt_handle2, txt_handle3])
        legend_text2.extend(["ONT Model: \n  {}".format(ont_model_name), "ONT Event Mean: {}".format(normal_mean), "ONT Event SD: {}".format(normal_sd)])

        if self.has_hdp_model:
            # plot HDP predicted distribution
            kmer_id = self.get_kmer_index(kmer)
            x = self.linspace
            panel1.set_xlim(min(x), max(x))
            hdp_y = self.all_posterior_pred[kmer_id]
            hdp_handle, = panel1.plot(x, hdp_y, '-')
            # compute entropy and hellinger distance
            ont_normal_dist = norm.pdf(self.linspace, normal_mean, normal_sd)
            kl_distance = entropy(pk=hdp_y, qk=ont_normal_dist, base=2)
            h_distance = hellinger2(p=hdp_y, q=ont_normal_dist)
            # deal with some extra text
            txt_handle4, = panel1.plot([], [], ' ')
            txt_handle5, = panel1.plot([], [], ' ')
            txt_handle6, = panel1.plot([], [], ' ')

            hdp_model_name = os.path.basename(self.hdp_path)

            handles1.append(hdp_handle)
            legend_text1.append("HDP Distribution")

            handles2.extend([txt_handle4, txt_handle5, txt_handle6])
            legend_text2.extend(["HDP Model: \n  {}".format(hdp_model_name), "Kullback–Leibler divergence: {}".format(np.round(kl_distance, 4)),
                                "Hellinger distance: {}".format(np.round(h_distance, 4))])

        if alignment_file is not None or alignment_file_data is not None:
            # option to parse file or not
            if alignment_file is not None:
                data = parse_assignment_file(alignment_file)
            else:
                data = alignment_file_data

            kmer_assignments = data.loc[data['kmer'] == kmer]
            kmer_data = kmer_assignments["level_mean"]
            # get event means and linspace in correct format
            x = np.asarray(kmer_data).reshape(len(kmer_data), 1)
            x_plot = self.linspace[:, np.newaxis]
            # get estimate for data
            if len(kmer_data) > 0:

                kde = KernelDensity(kernel="gaussian", bandwidth=0.5).fit(x)
                # estimate across the linspace
                log_dens = kde.score_samples(x_plot)
                kde_handle, = panel1.plot(x_plot[:, 0], np.exp(log_dens), '-')
                raw_data_handle, = panel1.plot(x[:, 0], -0.005 - 0.01 * np.random.random(x.shape[0]), '+k')
                # add to legend
                handles1.extend([kde_handle, raw_data_handle])
                legend_text1.extend(["Gaussian KDE Estimate", "Event Means: {} points".format(len(kmer_data))])
                txt_handle7, = panel1.plot([], [], ' ')
                if alignment_file:
                    alignment_file_name = os.path.basename(alignment_file)
                    handles2.append(txt_handle7)
                    legend_text2.append("RAW event data file: \n  {}".format(alignment_file_name))
            else:
                print("{} not found in alignment file".format(kmer))
        # create legend
        first_legend = panel1.legend(handles1, legend_text1, fancybox=True, shadow=True,
                                     loc='lower left', bbox_to_anchor=(1, .8))
        ax = plt.gca().add_artist(first_legend)

        panel1.legend(handles2, legend_text2, loc='upper left', bbox_to_anchor=(1, 0.2))

        # option to save figure or just show it
        if savefig_dir:
            base_name = "DNA_"
            if self.rna:
                base_name = "RNA_"
            name = "{}{}.png".format(base_name, kmer)
            out_path = os.path.join(savefig_dir, name)
            plt.savefig(out_path)
        else:
            plt.show()

    def get_kl_divergence(self, kmer):
        """Get Kullback–Leibler divergence between the HDP and ONT models for a specific kmer"""
        kmer_id = self.get_kmer_index(kmer)
        hdp_y = self.all_posterior_pred[kmer_id]
        if len(hdp_y) == 0:
            # print("[Kullback–Leibler divergence] No HDP data for {}".format(kmer))
            return None

        normal_mean, normal_sd = self.get_event_mean_gaussian_parameters(kmer)
        ont_normal_dist = norm.pdf(self.linspace, normal_mean, normal_sd)
        kl_divergence = entropy(pk=hdp_y, qk=ont_normal_dist, base=2)
        if kl_divergence == np.inf:
            # print("[Kullback–Leibler divergence] Zero probability for {}".format(kmer))
            return None

        return kl_divergence

    def get_hellinger_distance(self, kmer):
        """Get Hellinger distance between the HDP and ONT models for a specific kmer"""
        kmer_id = self.get_kmer_index(kmer)
        hdp_y = self.all_posterior_pred[kmer_id]
        if len(hdp_y) == 0:
            # print("[Hellinger Distance] No HDP data for {}".format(kmer))
            return None
        normal_mean, normal_sd = self.get_event_mean_gaussian_parameters(kmer)
        ont_normal_dist = norm.pdf(self.linspace, normal_mean, normal_sd)
        h_distance = hellinger2(p=hdp_y, q=ont_normal_dist)
        return h_distance

    def get_median_delta(self, kmer):
        """Calculate the difference between the max value of HDP and ONT kmer distributions"""
        kmer_id = self.get_kmer_index(kmer)
        hdp_y = self.all_posterior_pred[kmer_id]
        if len(hdp_y) == 0:
            # print("[Median Delta] No HDP data for {}".format(kmer))
            return None
        normal_mean, normal_sd = self.get_event_mean_gaussian_parameters(kmer)
        delta = self.linspace[hdp_y.index(max(hdp_y))]-normal_mean
        return abs(delta)

    def compare_distributions(self):
        """Calculate hellinger divergence and kl divergence between the HDP and ONT models for each kmer"""
        hellinger_distances = []
        kl_divergences = []
        median_deltas = []
        for kmer in self.sorted_kmer_tuple:
            # if statements used if the HDP model does not have information on the kmer distribution
            h_dist = self.get_hellinger_distance(kmer)
            if h_dist:
                hellinger_distances.append(h_dist)

            kl_divergence = self.get_kl_divergence(kmer)
            if kl_divergence:
                kl_divergences.append(kl_divergence)
                # print(kmer, kl_divergence)
            median_delta = self.get_median_delta(kmer)
            if median_delta:
                if len(median_deltas) > 0 and median_delta > max(median_deltas):
                    pass
                    # print(kmer, median_delta)
                median_deltas.append(median_delta)

        return hellinger_distances, kl_divergences, median_deltas

    def write_new_model(self, out_path, alphabet, replacement_base):
        """Write a correctly formatted new model file with a new alphabet.
        :param out_path: path to output hmm
        :param alphabet: new alphabet
        :param replacement_base: base to replace new character

        note: will retain same kmer size and assumes only one new character
        """
        # the model file has the format:
        # line 0: stateNumber \t alphabetSize \t alphabet \t kmerLength
        # line 1: match->match \t match->gapX \t match->gapY \t
        #         gapX->match \t gapX->gapX \t gapX->gapY \t
        #         gapY->match \t gapY->gapX \t gapY->gapY \n
        # line 2: [level_mean] [level_sd] [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n
        assert self.has_ont_model, "Shouldn't be writing down a Hmm that has no Model"
        if not self.normalized:
            self.normalize_transitions_expectations()

        alphabet = "".join(sorted(alphabet.upper()))
        for base in alphabet:
            assert not is_non_canonical_iupac_base(base), \
                "You cannot use IUPAC character to represent multiple bases. {}".format(base)

        replacement_base = replacement_base.upper()
        new_base = (set(alphabet) - set(self.alphabet)).pop()

        alphabet_size = len(alphabet)
        new_kmers = all_string_permutations(alphabet)
        with open(out_path, 'w') as f:

            # line 0
            f.write("{stateNumber}\t{alphabetSize}\t{alphabet}\t{kmerLength}\n"
                    "".format(stateNumber=self.state_number, alphabetSize=alphabet_size,
                              alphabet=alphabet, kmerLength=self.kmer_length))
            # line 1 transitions
            for i in range(self.state_number * self.state_number):
                f.write("{transition}\t".format(transition=str(self.transitions[i])))
            # likelihood
            f.write("{}\n".format(str(self.likelihood)))

            # line 2 Event Model
            for kmer in new_kmers:
                generic_kmer = kmer.replace(new_base, replacement_base)
                k = self.get_kmer_index(generic_kmer)
                f.write("{level_mean}\t{level_sd}\t{noise_mean}\t{noise_sd}\t{noise_lambda}\t"
                        "".format(level_mean=self.event_model["means"][k], level_sd=self.event_model["SDs"][k],
                                  noise_mean=self.event_model["noise_means"][k], noise_sd=self.event_model["noise_SDs"][k],
                                  noise_lambda=self.event_model["noise_lambdas"][k]))
            f.write("\n")


def hellinger2(p, q):
    return euclidean(np.sqrt(p), np.sqrt(q)) / _SQRT2


def parse_alignment_file(file_path):
    """Parse the buildAlignment.tsv output file from CreateHdpTrainingData

    :param file_path: path to alignment file
    :return: panda DataFrame with column names "kmer", "strand", "level_mean", "prob"
    """
    assert os.path.exists(file_path), "File path does not exist: {}".format(file_path)
    data = pd.read_table(file_path,
                         usecols=(4, 9, 12, 13),
                         names=["strand", "kmer", "prob", "level_mean"],
                         dtype={"kmer": np.str, "strand": np.str, "level_mean": np.float64, "prob": np.float64},
                         header=None)
    return data


def main():
    rna_ont_model = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/models/testModelR9p4_5mer_acgt_RNA.model"
    test_output_dir = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/test_plot_distributions"
    rna_hdp_model = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/models/template.singleLevelFixedCanonical.nhdp"
    savefig_dir = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/test_plot_distributions"
    dna_ont_model = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/models/testModelR9_5mer_acegit_template.model"
    dna_hdp_model = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/models/template.singleLevelPriorEcoli.nhdp"
    # plot_rna_dna_ont_hdp_comparison(ont_model, hdp_model, dna_ont_model, dna_hdp_model)#, savefig_dir=savefig_dir)
    alignment_file = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/test_plot_distributions/buildAlignment.tsv"
    # dna_hmm_handle = HmmModel(dna_ont_model, dna_hdp_model)
    #
    # dna_hmm_handle.plot_kmer_distribution("TAAAA", alignment_file=alignment_file) #savefig_dir=savefig_dir

    rna_hmm_handle = HmmModel(rna_ont_model, rna_hdp_model, rna=True)

    # rna_hmm_handle.plot_kmer_distribution("GGACT", alignment_file=alignment_file, savefig_dir=savefig_dir)

    # ﻿ RRACH
    #     R= A/G
    #     H = A/C/T




if __name__ == "__main__":
    main()
    raise SystemExit
