"""hiddenMarkovModel.py contains objects for handling HMMs for SignalAlign"""


from __future__ import print_function
import sys
import os
import numpy as np
from scipy.stats import norm, invgauss
# Globals
NORM_DIST_PARAMS = 2
NB_MODEL_PARAMS = 5


class HmmModel(object):
    def __init__(self, model_file):
        # TODO Need to create docs here
        self.match_model_params = 5  # level_mean, level_sd, noise_mean, noise_sd, noise_lambda
        self.state_number = 3
        self.transitions = np.zeros(self.state_number**2)
        self.transitions_expectations = np.zeros(self.state_number**2)
        self.likelihood = 0.0
        self.running_likelihoods = []
        self.alphabet_size = 0
        self.alphabet = ""
        self.kmer_length = 0
        self.kmer_index = dict()
        self.has_model = False
        self.normalized = False
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
        self.load_model(model_file)
        self.mean_expectations = np.zeros(self.symbol_set_size)
        self.sd_expectations = np.zeros(self.symbol_set_size)
        self.posteriors = np.zeros(self.symbol_set_size)
        self.observed = np.zeros(self.symbol_set_size, dtype=bool)


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
        MATCH_FROM_GAP_X = np.exp(-0.21880828092192281)   # 1 - skip'
        MATCH_FROM_GAP_Y = np.exp(-0.013406326748077823)  # 1 - (skip + stay)
        GAP_OPEN_X = np.exp(-1.6269694202638481)          # skip
        GAP_OPEN_Y = np.exp(-4.3187242127300092)          # 1 - (skip + stride)
        GAP_EXTEND_X = np.exp(-1.6269694202638481)        # skip'
        GAP_EXTEND_Y = np.exp(-4.3187242127239411)        # stay (1 - (skip + stay))
        GAP_SWITCH_TO_X = 0.000000001
        GAP_SWITCH_TO_Y = 0.0
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
            self.has_model = True

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
        assert self.has_model, "Shouldn't be writing down a Hmm that has no Model"
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

    def get_kmer_index(self, kmer):
        """Get the model index for a given kmer

        ex: get_kmer_index(AAAAA) = 0
        :param kmer: nucleotide sequence
        """
        assert set(kmer).issubset(set(self.alphabet)) is True, "Nucleotide not found in model alphabet: kmer={}, " \
                                                               "alphabet={}".format(kmer, self.alphabet)
        assert len(kmer) == self.kmer_length, "Kmer length does not match model kmer length"

        alphabet_dict = {base: index for index, base in enumerate(sorted(self.alphabet))}
        kmer_index = 0
        for index, nuc in enumerate(kmer):
            kmer_index += alphabet_dict[nuc]*(self.alphabet_size**(self.kmer_length-index-1))
        return kmer_index

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
