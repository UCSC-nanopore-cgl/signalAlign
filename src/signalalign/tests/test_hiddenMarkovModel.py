#!/usr/bin/env python3
"""
    Place unit tests for hiddenMarkovModel.py
"""
########################################################################
# File: test_hiddenMarkovModel.py
#  executable: test_hiddenMarkovModel.py
# Purpose: test hiddenMarkovModel
#
# Author: Andrew Bailey
# History: 5/31/2018 Created
########################################################################
import unittest
import os
import numpy as np
import threading
import time
import tempfile
from shutil import copyfile

from signalalign.hiddenMarkovModel import *
from signalalign.fast5 import Fast5
from py3helpers.utils import all_string_permutations, get_random_string


class HiddenMarkovTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(HiddenMarkovTests, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.model_file = os.path.join(cls.HOME, "models/testModelR9p4_5mer_acgt_RNA.model")
        cls.r9_model_file = os.path.join(cls.HOME, "models/testModelR9_acegt_complement.model")

        cls.model = HmmModel(model_file=cls.model_file)
        cls.expectation_file = os.path.join(cls.HOME,
                                            "tests/test_expectation_files/4f9a316c-8bb3-410a-8cfc-026061f7e8db.template.expectations.tsv")

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

    def test_HmmModel(self):
        hdp_model_file = os.path.join(self.HOME, "models/testModelR9p4_5mer_acgt_RNA.model")
        model = HmmModel(model_file=hdp_model_file)
        self.assertIsInstance(model, HmmModel)

        model = HmmModel(model_file=self.model_file)
        self.assertIsInstance(model, HmmModel)

    def test_add_expectations_file(self):
        model = HmmModel(model_file=self.r9_model_file)
        model.add_expectations_file(self.expectation_file)
        model = HmmModel(model_file=self.r9_model_file)
        model.add_expectations_file(self.expectation_file)
        self.assertRaises(AssertionError, self.model.add_expectations_file, self.expectation_file)

    def test_check_header_line(self):
        self.model.check_header_line(['3', '4', "ACGT", '5'], "path")
        self.assertRaises(AssertionError, self.model.check_header_line, ['1', '4', "ACGT", '5'], "ssomething")
        self.assertRaises(AssertionError, self.model.check_header_line, ['3', '2', "ACGT", '5'], "ssomething")

    def test_add_and_normalize_expectations(self):
        with tempfile.TemporaryDirectory() as tempdir:
            test_expecations_file = os.path.join(tempdir, "fake.expectations.tsv")
            copyfile(self.expectation_file, test_expecations_file)
            files = [test_expecations_file]
            model = HmmModel(model_file=self.r9_model_file)
            model.add_and_normalize_expectations(files, os.path.join(tempdir, "fake.hmm"))

    def test_normalize_transitions_expectations(self):
        hdp_model_file = os.path.join(self.HOME, "models/testModelR9_acegt_complement.model")
        model = HmmModel(model_file=hdp_model_file)
        model.add_expectations_file(self.expectation_file)
        model.add_expectations_file(self.expectation_file)
        model.add_expectations_file(self.expectation_file)
        model.normalize_transitions_expectations()
        for from_state in range(model.state_number):
            i = model.state_number * from_state
            self.assertAlmostEqual(sum(model.transitions_expectations[i:i + model.state_number]), 1)

    def test_write(self):
        with tempfile.TemporaryDirectory() as tempdir:
            test_hmm_file = os.path.join(tempdir, "fake.model.hmm")
            model = HmmModel(model_file=self.model_file)
            self.assertRaises(AssertionError, model.write, test_hmm_file)
            model.normalized = True
            model.write(test_hmm_file)

    def test_normalise(self):
        model = HmmModel(model_file=self.r9_model_file)
        model.add_expectations_file(self.expectation_file)
        model.add_expectations_file(self.expectation_file)
        model.add_expectations_file(self.expectation_file)
        model.normalize(update_transitions=True, update_emissions=False)
        model.normalize(update_emissions=True, update_transitions=True)
        self.assertTrue(model.normalized)

    def test_reset_assignments(self):
        model = HmmModel(model_file=self.r9_model_file)
        model.add_expectations_file(self.expectation_file)
        model.reset_assignments()
        self.assertSequenceEqual(model.event_assignments, [])
        self.assertSequenceEqual(model.kmer_assignments, [])


if __name__ == '__main__':
    unittest.main()
