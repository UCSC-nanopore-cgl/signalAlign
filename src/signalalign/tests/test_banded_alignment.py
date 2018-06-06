#!/usr/bin/env python
"""
    Place unit tests for banded_alignment.py
"""
########################################################################
# File: banded_alignment.py
#  executable: banded_alignment.py
# Purpose: banded alignment functions
#
# Author: Andrew Bailey
# History: 5/31/2018 Created
########################################################################
import unittest
import os
import numpy as np
import threading
import time
from signalalign.hiddenMarkovModel import *
from signalalign.banded_alignment import *
from py3helpers.utils import all_string_permutations, get_random_string


class HiddenMarkovTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(HiddenMarkovTests, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-2])
        # fast5_file = os.path.join(cls.HOME, "tests/minion_test_reads/canonical_ecoli_R9/miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_strand.fast5")
        cls.model_file = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/models/testModelR9p4_5mer_acgt_RNA.model"
        cls.model = get_model(model_type="threeState", model_file=cls.model_file)

    def test_get_kmer_index(self):
        all_kmers = [x for x in all_string_permutations("ATGC", 5)]
        for x in range(10):
            kmer = get_random_string(5, chars="ATGC")
            self.assertEqual(all_kmers.index(kmer), self.model.get_kmer_index(kmer))

    def test_log_event_mean_gaussian_probability_match(self):
        def emissions_signal_logGaussPdf(x, mu, sigma):
            log_inv_sqrt_2pi = -0.91893853320467267
            l_sigma = np.log(sigma)
            a = (x - mu) / sigma
            # // returns Log-space
            return log_inv_sqrt_2pi - l_sigma + (-0.5 * a * a)

        for x in range(10):
            kmer = get_random_string(5, chars="ATGC")
            mu, sigma = self.model.get_event_mean_gaussian_parameters(kmer)
            prob = self.model.log_event_mean_gaussian_probability_match(50, kmer)
            self.assertAlmostEqual(prob, emissions_signal_logGaussPdf(50, mu, sigma))

    def test_log_event_sd_inv_gaussian_probability_match(self):
        def emissions_signal_logInvGaussPdf(eventNoise, modelNoiseMean, modelNoiseLambda):
            l_twoPi = 1.8378770664093453  # // log(2*pi)
            l_eventNoise = np.log(eventNoise)
            a = (eventNoise - modelNoiseMean) / modelNoiseMean
            l_modelNoseLambda = np.log(modelNoiseLambda)
            # // returns Log-space
            return (l_modelNoseLambda - l_twoPi - 3 * l_eventNoise - (modelNoiseLambda * a * a / eventNoise)) / 2

        for x in range(10):
            kmer = get_random_string(5, chars="ATGC")
            mu, lambda1 = self.model.get_event_sd_inv_gaussian_parameters(kmer)
            prob = self.model.log_event_sd_inv_gaussian_probability_match(2, kmer)
            self.assertAlmostEqual(prob, emissions_signal_logInvGaussPdf(2, mu, lambda1))

    def test_get_event_mean_gaussian_parameters(self):
        for x in range(10):
            kmer = get_random_string(5, chars="ATGC")
            mu, sigma = self.model.get_event_mean_gaussian_parameters(kmer)
        mu, sigma = self.model.get_event_mean_gaussian_parameters("TTTTT")
        self.assertEqual(self.model.event_model["means"][-1], mu)
        self.assertEqual(self.model.event_model["SDs"][-1], sigma)

    def test_get_event_sd_inv_gaussian_parameters(self):
        for x in range(10):
            kmer = get_random_string(5, chars="ATGC")
            mu, sigma = self.model.get_event_sd_inv_gaussian_parameters(kmer)
        mean, lambda1 = self.model.get_event_sd_inv_gaussian_parameters("TTTTT")
        self.assertEqual(self.model.event_model["noise_means"][-1], mean)
        self.assertEqual(self.model.event_model["noise_lambdas"][-1], lambda1)


if __name__ == '__main__':
    unittest.main()
