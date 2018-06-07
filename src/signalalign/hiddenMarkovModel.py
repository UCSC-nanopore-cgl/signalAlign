"""hiddenMarkovModel.py contains objects for handling HMMs for SignalAlign"""


from __future__ import print_function
import sys
import os
import numpy as np
from scipy.stats import norm, invgauss
# Globals
NORM_DIST_PARAMS = 2
NB_MODEL_PARAMS = 5

#emissions_signal_strawManGetKmerEventMatchProbWithDescaling

# int64_t kmerIndex = kmer_id(kmer_i, self->model.alphabet, self->model.alphabetSize, self->model.kmerLength);
#
# double *eventModel = match ? self->model.EMISSION_MATCH_MATRIX : self->model.EMISSION_GAP_Y_MATRIX;
# // get the µ and σ for the level and noise for the model
# double levelMean = emissions_signal_getModelLevelMean(eventModel, kmerIndex);
# double levelStdDev = emissions_signal_getModelLevelSd(eventModel, kmerIndex);
# eventMean = emissions_signal_descaleEventMean_JordanStyle(eventMean, levelMean, self->model.scale,
#                                                                                       self->model.shift, self->model.var);
#
# double noiseMean = emissions_signal_getModelFluctuationMean(eventModel, kmerIndex);
# //double noiseStdDev = emissions_signal_getModelFluctuationSd(eventModel, kmerIndex);
#
# double modelNoiseLambda = emissions_signal_getModelFluctuationLambda(eventModel, kmerIndex);
#
# double l_probEventMean = emissions_signal_logGaussPdf(eventMean, levelMean, levelStdDev);
#
# //double l_probEventNoise = emissions_signal_logGaussPdf(eventNoise, noiseMean, noiseStdDev);
# double l_probEventNoise = emissions_signal_logInvGaussPdf(eventNoise, noiseMean, modelNoiseLambda);
#
# // clean
# free(kmer_i);
#
# // debugging
# //double prob = l_probEventMean + l_probEventNoise;
# //st_uglyf("MATCHING--x_i:%s (index: %lld), e_j mean: %f, \n modelMean: %f, modelLsd: %f probEvent: %f probNoise: %f, combined: %f\n",
#   //         kmer_i, kmerIndex, eventMean, levelMean, levelStdDev, l_probEventMean, l_probEventNoise, prob);
#
# return l_probEventMean + l_probEventNoise;
# }

class SignalHmm(object):
    def __init__(self, model_type):
        self.match_model_params = 5  # level_mean, level_sd, noise_mean, noise_sd, noise_lambda
        self.model_type = model_type  # ID of model type
        self.state_number = {"threeState": 3, "threeStateHdp": 3}[model_type]
        self.symbol_set_size = 0
        self.transitions = np.zeros(self.state_number**2)
        self.transitions_expectations = np.zeros(self.state_number**2)
        self.likelihood = 0.0
        self.running_likelihoods = []
        self.alphabet_size = 0
        self.alphabet = ""
        self.kmer_length = 0
        self.has_model = False
        self.normalized = False
        self.kmer_index = dict()
        # event model for describing normal distributions for each kmer
        self.event_model = {"means": np.zeros(self.symbol_set_size),
                            "SDs": np.zeros(self.symbol_set_size),
                            "noise_means": np.zeros(self.symbol_set_size),
                            "noise_SDs": np.zeros(self.symbol_set_size),
                            "noise_lambdas": np.zeros(self.symbol_set_size)}

    def normalize_transitions_expectations(self):
        # normalize transitions
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
        if len(line) != 4:
            print("signalHmm.check_header_line - incorrect header (param line): {}".format(expectations_file), file=sys.stderr)
            return False
        if int(line[0]) != self.state_number:
            print("signalHmm.check_header_line - state number error should be {exp} got {obs}"
                  "".format(exp=self.state_number, obs=line[0]), file=sys.stderr)
            return False
        if int(line[1]) != self.alphabet_size:
            print("signalHmm.check_header_line - alphabet size error incorrect parameters: {file}, line {line}"
                  "".format(file=expectations_file, line=''.join(line)), file=sys.stderr)
            return False
        if line[2] != self.alphabet:
            print("signalHmm.check_header_line - incorrect parameters: {file}, line {line}"
                  "".format(file=expectations_file, line=''.join(line)), file=sys.stderr)
            return False
        if int(line[3]) != self.kmer_length:
            print("signalHmm.check_header_line - incorrect parameters: {file}, line {line}"
                  "".format(file=expectations_file, line=''.join(line)), file=sys.stderr)
            return False
        return True

    def load_model(self, model_file):
        # the model file has the format:
        # line 0: stateNumber \t alphabetSize \t alphabet \t kmerLength
        # line 1: match->match \t match->gapX \t match->gapY \t
        #         gapX->match \t gapX->gapX \t gapX->gapY \t
        #         gapY->match \t gapY->gapX \t gapY->gapY \n
        # line 2: [level_mean] [level_sd] [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n
        assert os.path.exists(model_file), "signalHmm.load_model - didn't find model here{}?".format(model_file)

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
        # the model file has the format:
        # line 0: stateNumber \t alphabetSize \t alphabet \t kmerLength
        # line 1: match->match \t match->gapX \t match->gapY \t
        #         gapX->match \t gapX->gapX \t gapX->gapY \t
        #         gapY->match \t gapY->gapX \t gapY->gapY \n
        # line 2: [level_mean] [level_sd] [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n
        assert self.has_model, "Shouldn't be writing down a Hmm that has no Model"
        assert self.normalized, "Shouldn't be writing down a not normalized HMM"

        f = open(out_file, 'w')

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

        f.close()

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

    # def get_event_probability(self, mean, stdev, kmer):
    #     """Get probability of an event aligning to a specific kmer"""



class ContinuousPairHmm(SignalHmm):
    def __init__(self, model_type):
        super(ContinuousPairHmm, self).__init__(model_type=model_type)
        self.set_default_transitions()

        # bins for expectations
        self.mean_expectations = np.zeros(self.symbol_set_size)
        self.sd_expectations = np.zeros(self.symbol_set_size)
        self.posteriors = np.zeros(self.symbol_set_size)
        self.observed = np.zeros(self.symbol_set_size, dtype=bool)
        self.has_model = False
        self.normalized = False

    def add_expectations_file(self, expectations_file):
        # expectations files have the format:
        # line 0: stateNumber \t alphabetSize \t alphabet \t kmerLength
        # line 1: match->match \t match->gapX \t match->gapY \t
        #         gapX->match \t gapX->gapX \t gapX->gapY \t
        #         gapY->match \t gapY->gapX \t gapY->gapY \n
        # line 2: [level_mean] [level_sd] [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n
        # line 3: event expectations [mean] [sd] / kmer \n
        # line 4: posteriors 1 per kmer \n
        # line 5: observed 1 per kmer \n
        if not os.path.exists(expectations_file) or os.stat(expectations_file).st_size == 0:
            print("Empty or missing file {}".format(expectations_file), file=sys.stderr)
            return False

        fH = open(expectations_file, 'r')

        # line 0
        line = fH.readline().split()
        header_line_check = self.check_header_line(line=line, expectations_file=expectations_file)

        if header_line_check is False:
            fH.close()
            return False

        # line 1: transitions, likelihood
        line = list(map(float, fH.readline().split()))
        # check if valid
        if len(line) != (len(self.transitions) + 1):
            print("cpHMM: check_file - bad file (transitions expectations): {}".format(expectations_file),
                  file=sys.stderr)
            fH.close()
            return False

        self.likelihood += line[-1]
        self.transitions_expectations = [sum(x) for x in zip(self.transitions_expectations, line[0:-1])]
        # line 2: event model
        line = list(map(float, fH.readline().split()))
        if len(line) != self.symbol_set_size * NB_MODEL_PARAMS:
            print("cpHMM: check_file - bad file (event model): {}".format(expectations_file), file=sys.stderr)
            fH.close()
            return False

        # line 3 event expectations [E_mean, E_sd]
        line = list(map(float, fH.readline().split()))
        if len(line) != self.symbol_set_size * NORM_DIST_PARAMS:
            print("cpHMM: check_file - bad file (event expectations): {}".format(expectations_file), file=sys.stderr)
            fH.close()
            return False

        self.mean_expectations = [i + j for i, j in zip(self.mean_expectations, line[::NORM_DIST_PARAMS])]
        self.sd_expectations = [i + j for i, j in zip(self.sd_expectations, line[1::NORM_DIST_PARAMS])]

        # line 4, posteriors
        line = list(map(float, fH.readline().split()))
        if len(line) != self.symbol_set_size:
            print("cpHMM: check_file - bad file (posteriors): {}".format(expectations_file), file=sys.stderr)
            fH.close()
            return False

        self.posteriors = [sum(x) for x in zip(self.posteriors, line)]
        line = list(map(bool, fH.readline().split()))
        if len(line) != self.symbol_set_size:
            print("cpHMM: check_file - bad file (observations): {}".format(expectations_file), file=sys.stderr)
            fH.close()
            return False

        self.observed = [any(b) for b in zip(self.observed, line)]

        fH.close()
        return True

    def normalize(self, update_transitions, update_emissions):
        # normalize transitions expectations
        self.normalize_transitions_expectations()

        # update
        if update_transitions is True:
            for i in range(self.state_number**2):
                self.transitions[i] = self.transitions_expectations[i]

        # calculate the new expected mean and standard deviation for the kmer normal distributions
        if update_emissions:
            for k in range(self.symbol_set_size):  # TODO implement learning rate
                if self.observed[k] is True:
                    u_k = self.mean_expectations[k] / self.posteriors[k]
                    o_k = np.sqrt(self.sd_expectations[k] / self.posteriors[k])
                    if u_k > 0:
                        self.event_model["means"][k] = u_k
                        self.event_model["SDs"][k] = o_k
                else:
                    continue
        self.normalized = True


class HdpSignalHmm(SignalHmm):
    def __init__(self, model_type):
        super(HdpSignalHmm, self).__init__(model_type=model_type)
        self.set_default_transitions()
        self.kmer_assignments = []
        self.event_assignments = []
        self.assignments_record = []

    def add_expectations_file(self, expectations_file):
        # expectations files have the format:
        # line 0: stateNumber \t alphabetSize \t alphabet \t kmerLength
        # line 1: match->match \t match->gapX \t match->gapY \t
        #         gapX->match \t gapX->gapX \t gapX->gapY \t
        #         gapY->match \t gapY->gapX \t gapY->gapY \n
        # line 2: [level_mean] [level_sd] [noise_mean] [noise_sd] [noise_lambda ](.../kmer) \n
        # line 3: event assignments
        # line 4: kmer assignments
        if not os.path.exists(expectations_file) or os.stat(expectations_file).st_size == 0:
            print("Empty or missing file {}".format(expectations_file))
            return

        fH = open(expectations_file, 'r')

        # line 0
        line = fH.readline().split()
        header_line_check = self.check_header_line(line=line, expectations_file=expectations_file)
        if header_line_check is False:
            fH.close()
            return False

        # line 1: transitions, likelihood
        line = list(map(float, fH.readline().split()))
        if len(line) != (len(self.transitions) + 1):
            print("hdpHmm.add_expectations_file - problem with file {f} transitions line {l}, incorrect length"
                  "".format(f=expectations_file, l=''.join(line)), file=sys.stdout)
            fH.close()
            return False

        self.likelihood += line[-1]
        self.transitions_expectations = [sum(x) for x in zip(self.transitions_expectations, line[0:-1])]

        # line 2: event model
        line = list(map(float, fH.readline().split()))
        if len(line) != self.symbol_set_size * NB_MODEL_PARAMS:
            print("hdpHmm.add_expectations_file - problem with event model in file {}"
                  "".format(expectations_file), file=sys.stderr)
            fH.close()
            return False

        # line 3: event assignments
        line = list(map(float, fH.readline().split()))
        self.event_assignments += line

        # line 4: kmer assignments
        line = list(map(str, fH.readline().split()))
        self.kmer_assignments += line

        fH.close()
        return True

    def reset_assignments(self):
        self.assignments_record.append(len(self.event_assignments))
        self.event_assignments = []
        self.kmer_assignments = []

    def normalize(self, update_transitions, update_emissions=None):
        self.normalize_transitions_expectations()
        if update_transitions is True:
            for i in range(self.state_number**2):
                self.transitions[i] = self.transitions_expectations[i]
        self.normalized = True
