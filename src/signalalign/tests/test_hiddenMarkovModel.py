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
import tempfile
from shutil import copyfile

from signalalign.hiddenMarkovModel import *
from py3helpers.utils import all_string_permutations, get_random_string, list_dir


class HiddenMarkovTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(HiddenMarkovTests, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.model_file = os.path.join(cls.HOME, "models/testModelR9p4_5mer_acgt_RNA.model")
        cls.r9_model_file = os.path.join(cls.HOME, "models/testModelR9_acegt_complement.model")

        cls.model = HmmModel(ont_model_file=cls.model_file)
        cls.expectation_file = os.path.join(cls.HOME,
                                            "tests/test_expectation_files/4f9a316c-8bb3-410a-8cfc-026061f7e8db.template.expectations.tsv")
        cls.nanopolish_model = os.path.join(cls.HOME, "models/r9.4_450bps.nucleotide.6mer.template.model")
        cls.cpg_nanopolish_model = os.path.join(cls.HOME, "models/r9.4_450bps.cpg.6mer.template.model")

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
        model = HmmModel(ont_model_file=hdp_model_file)
        self.assertIsInstance(model, HmmModel)

        model = HmmModel(ont_model_file=self.model_file)
        self.assertIsInstance(model, HmmModel)

    def test_add_expectations_file(self):
        model = HmmModel(ont_model_file=self.r9_model_file)
        model.add_expectations_file(self.expectation_file)
        model = HmmModel(ont_model_file=self.r9_model_file)
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
            model = HmmModel(ont_model_file=self.r9_model_file)
            model.add_and_normalize_expectations(files, os.path.join(tempdir, "fake.hmm"))

    def test_normalize_transitions_expectations(self):
        hdp_model_file = os.path.join(self.HOME, "models/testModelR9_acegt_complement.model")
        model = HmmModel(ont_model_file=hdp_model_file)
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
            model = HmmModel(ont_model_file=self.model_file)
            self.assertRaises(AssertionError, model.write, test_hmm_file)
            model.normalized = True
            model.write(test_hmm_file)

    def test_normalise(self):
        model = HmmModel(ont_model_file=self.r9_model_file)
        model.add_expectations_file(self.expectation_file)
        model.add_expectations_file(self.expectation_file)
        model.add_expectations_file(self.expectation_file)
        model.normalize(update_transitions=True, update_emissions=False)
        model.normalize(update_emissions=True, update_transitions=True)
        self.assertTrue(model.normalized)

    def test_reset_assignments(self):
        model = HmmModel(ont_model_file=self.r9_model_file)
        model.add_expectations_file(self.expectation_file)
        model.reset_assignments()
        self.assertSequenceEqual(model.event_assignments, [])
        self.assertSequenceEqual(model.kmer_assignments, [])

    def test_HDP_model_load(self):
        hdp_model = os.path.join(self.HOME, "models/template_RNA.singleLevelFixedCanonical.nhdp")
        hdp_handle = HmmModel(ont_model_file=self.model_file, hdp_model_file=hdp_model)
        kmer = "AACAT"
        kmer_id = 19
        self.assertEqual(kmer_id, hdp_handle.get_kmer_index(kmer))
        query_x = 83.674161662792542
        x = hdp_handle.linspace
        y = hdp_handle.all_posterior_pred[kmer_id]
        slope = hdp_handle.all_spline_slopes[kmer_id]
        length = hdp_handle.grid_length
        prob = hdp_handle.grid_spline_interp(query_x, x, y, slope, length)
        self.assertEqual(prob, 0.29228949476718646)
        query_x = 81.55860779063407
        prob = hdp_handle.grid_spline_interp(query_x, x, y, slope, length)
        self.assertEqual(prob, 0.12927539337648492)
        kmer = "CATTT"
        kmer_id = 319
        self.assertEqual(kmer_id, hdp_handle.get_kmer_index(kmer))
        y = hdp_handle.all_posterior_pred[kmer_id]
        slope = hdp_handle.all_spline_slopes[kmer_id]
        query_x = 80.605230545769458
        prob = hdp_handle.grid_spline_interp(query_x, x, y, slope, length)
        self.assertEqual(prob, 0.12328410496683605)

    def test_get_hdp_probability(self):
        hdp_model = os.path.join(self.HOME, "models/template_RNA.singleLevelFixedCanonical.nhdp")
        hdp_handle = HmmModel(ont_model_file=self.model_file, hdp_model_file=hdp_model)
        query_x = 83.674161662792542
        prob = hdp_handle.get_hdp_probability("AACAT", query_x)
        self.assertEqual(prob, 0.29228949476718646)
        query_x = 81.55860779063407
        prob = hdp_handle.get_hdp_probability("AACAT", query_x)
        self.assertEqual(prob, 0.12927539337648492)
        query_x = 80.605230545769458
        prob = hdp_handle.get_hdp_probability("CATTT", query_x)
        self.assertEqual(prob, 0.12328410496683605)

    def test_get_new_linspace_hdp_probability_distribution(self):
        hdp_model = os.path.join(self.HOME, "models/template_RNA.singleLevelFixedCanonical.nhdp")
        hdp_handle = HmmModel(ont_model_file=self.model_file, hdp_model_file=hdp_model)
        kmer = "AACAT"
        linspace = hdp_handle.linspace
        kmer_id = hdp_handle.get_kmer_index(kmer)
        y = hdp_handle.all_posterior_pred[kmer_id]
        new_y = hdp_handle.get_new_linspace_hdp_probability_distribution(kmer, linspace)
        self.assertSequenceEqual(new_y, y)

    def test_write_new_model(self):
        with tempfile.TemporaryDirectory() as tempdir:
            test_model_file = os.path.join(tempdir, "fake.hmm")

            hmm_handle = HmmModel(ont_model_file=self.model_file)
            hmm_handle.write_new_model(out_path=test_model_file, alphabet="ATGCF", replacement_base="A")
            hmm_handle2 = HmmModel(ont_model_file=test_model_file)
            self.assertEqual(hmm_handle.kmer_length, hmm_handle2.kmer_length)
            self.assertEqual(hmm_handle2.alphabet, "ACFGT")
            self.assertEqual(hmm_handle2.alphabet_size, 5)
            self.assertRaises(AssertionError, hmm_handle.write_new_model, test_model_file, "ATGCW", "A")

    def test_create_new_model(self):
        with tempfile.TemporaryDirectory() as tempdir:
            test_model_file = os.path.join(tempdir, "fake.hmm")
            new_model = create_new_model(self.model_file, test_model_file, (("A", "F"), ("A", "J")))
            self.assertEqual(new_model.kmer_length, 5)
            self.assertEqual(new_model.alphabet, "ACFGJT")
            self.assertEqual(new_model.alphabet_size, 6)
            mean1 = new_model.get_event_mean_gaussian_parameters("AAAAA")
            mean2 = new_model.get_event_mean_gaussian_parameters("AAAJA")
            mean3 = new_model.get_event_mean_gaussian_parameters("AAAFA")
            mean4 = new_model.get_event_mean_gaussian_parameters("AAJFA")
            mean5 = new_model.get_event_mean_gaussian_parameters("AAJJJ")
            mean6 = new_model.get_event_mean_gaussian_parameters("FFJJJ")
            self.assertEqual(mean1, mean2)
            self.assertEqual(mean2, mean3)
            self.assertEqual(mean3, mean4)
            self.assertEqual(mean4, mean5)
            self.assertEqual(mean5, mean6)
            new_model = create_new_model(self.model_file, test_model_file, [("A", "J")])

    def test_set_kmer_event_mean(self):
        hmm_handle = HmmModel(ont_model_file=self.model_file)
        hmm_handle.set_kmer_event_mean("AAAAA", 1000)
        mean, sd = hmm_handle.get_event_mean_gaussian_parameters("AAAAA")
        self.assertEqual(mean, 1000)

    def test_set_kmer_event_sd(self):
        hmm_handle = HmmModel(ont_model_file=self.model_file)
        hmm_handle.set_kmer_event_sd("AAAAA", 1000)
        mean, sd = hmm_handle.get_event_mean_gaussian_parameters("AAAAA")
        self.assertEqual(sd, 1000)

    def test_set_kmer_noise_means(self):
        hmm_handle = HmmModel(ont_model_file=self.model_file)
        hmm_handle.set_kmer_noise_means("AAAAA", 1000)
        mean, sd = hmm_handle.get_event_sd_inv_gaussian_parameters("AAAAA")
        self.assertEqual(mean, 1000)

    def test_set_kmer_noise_lambdas(self):
        hmm_handle = HmmModel(ont_model_file=self.model_file)
        hmm_handle.set_kmer_noise_lambdas("AAAAA", 1000)
        mean, sd = hmm_handle.get_event_sd_inv_gaussian_parameters("AAAAA")
        self.assertEqual(sd, 1000)

    def test_read_in_alignment_file(self):
        assignments_dir = os.path.join(self.HOME, "tests/test_alignments/ecoli1D_test_alignments_sm3")
        data = read_in_alignment_file(list_dir(assignments_dir)[0])
        self.assertEqual(len(data["contig"]), 16852)
        self.assertEqual(len(data["reference_index"]), 16852)
        self.assertEqual(len(data["reference_kmer"]), 16852)
        self.assertEqual(len(data["read_file"]), 16852)
        self.assertEqual(len(data["strand"]), 16852)
        self.assertEqual(len(data["event_index"]), 16852)
        self.assertEqual(len(data["event_mean"]), 16852)
        self.assertEqual(len(data["event_noise"]), 16852)
        self.assertEqual(len(data["event_duration"]), 16852)
        self.assertEqual(len(data["aligned_kmer"]), 16852)
        self.assertEqual(len(data["scaled_mean_current"]), 16852)
        self.assertEqual(len(data["scaled_noise"]), 16852)
        self.assertEqual(len(data["posterior_probability"]), 16852)
        self.assertEqual(len(data["descaled_event_mean"]), 16852)
        self.assertEqual(len(data["ont_model_mean"]), 16852)
        self.assertEqual(len(data["path_kmer"]), 16852)
        self.assertEqual(len(data), 16852)

    def test_load_nanopolish_model(self):
        # model = HmmModel(ont_model_file=self.model_file, nanopolish_model_file=nanopolish_model)
        model, alphabet, k = load_nanopolish_model(self.nanopolish_model)
        self.assertEqual(len(model["means"]), 4**6)
        self.assertEqual(alphabet, "ACGT")
        self.assertEqual(k, 6)

    def test_convert_nanopolish_model_to_signalalign(self):
        with tempfile.TemporaryDirectory() as tempdir:
            sa_file = os.path.join(tempdir, "testModelr9.4_450bps.nucleotide.6mer.template.model")
            convert_nanopolish_model_to_signalalign(self.nanopolish_model, self.model.transitions, sa_file)
            sa_model = HmmModel(sa_file)
            model_mean, model_sd = sa_model.get_event_mean_gaussian_parameters("AAAATG")
            self.assertEqual(model_mean, 75.943873)
            self.assertEqual(model_sd, 1.542528)

    def test_convert_and_edit_nanopolish_model_to_signalalign(self):
        with tempfile.TemporaryDirectory() as tempdir:
            sa_file = os.path.join(tempdir, "testModelR9.4_450bps.cpg.6mer.template.model")
            convert_and_edit_nanopolish_model_to_signalalign(self.cpg_nanopolish_model, self.model.transitions, sa_file)
            sa_model = HmmModel(sa_file)
            model_mean, model_sd = sa_model.get_event_mean_gaussian_parameters("AAAAEE")
            self.assertEqual(model_mean, 75.7063)
            self.assertEqual(model_sd, 2.70501)

    def test_parse_alignment_file(self):
        test_path = "tests/test_alignments/ecoli1D_test_alignments_RAW/5cc86bac-79fd-4897-8631-8f1c55954a45.sm.backward.tsv"
        full_path = os.path.join(self.HOME, test_path)
        data = parse_alignment_file(full_path)
        self.assertTrue(data["kmer"][0], "AGGTG")
        self.assertTrue(data["strand"][0], "t")
        self.assertTrue(data["level_mean"][0], 63.674513)
        self.assertTrue(data["prob"][0], 1)

if __name__ == '__main__':
    unittest.main()
