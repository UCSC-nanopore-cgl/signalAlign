#!/usr/bin/env python3
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
from signalalign.fast5 import Fast5
from py3helpers.utils import all_string_permutations, get_random_string


class HiddenMarkovTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(HiddenMarkovTests, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.model_file = os.path.join(cls.HOME, "models/testModelR9p4_5mer_acgt_RNA.model")
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


class BandedAlignmentTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(BandedAlignmentTests, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.rna_model_file = os.path.join(cls.HOME, "models/testModelR9p4_5mer_acgt_RNA.model")
        cls.dna_template_model_file = os.path.join(cls.HOME, "models/testModelR9p4_5mer_acegt_template.model")
        cls.rna_model = get_model(model_type="threeState", model_file=cls.rna_model_file)
        cls.dna_model = get_model(model_type="threeState", model_file=cls.dna_template_model_file)

        cls.dna_fast5_path = os.path.join(cls.HOME,
                                          "tests/minion_test_reads/1D/LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch112_read108_strand.fast5")
        cls.dna_handle = Fast5(cls.dna_fast5_path)

    def test_simple_banded_alignment1(self):
        # event_table = self.dna_handle.get_basecall_data()
        # fastq = self.dna_handle.get_fastq(analysis="Basecall_1D", section='template')
        # print(fastq)
        # nuc_sequence = fastq.split('\n')[1]
        # print(nuc_sequence)
        event_table = np.empty(3, dtype=[('start', float), ('length', float),
                                         ('mean', float), ('stdv', float),
                                         ('model_state', 'S5'), ('move', '<i4'),
                                         ('p_model_state', float)])

        event_table["mean"] = [86.9, 75.3, 80.2]
        nuc_sequence = "AAAAATA"

        events, sum_emission = simple_banded_event_align(event_table, self.dna_model, nuc_sequence)
        self.assertSequenceEqual(events['mean'].tolist(), [86.9, 75.3, 80.2])
        self.assertSequenceEqual([bytes.decode(x) for x in events['model_state']], ["AAAAA", "AAAAT", "AAATA"])
        self.assertSequenceEqual(events['move'].tolist(), [0, 1, 1])

    def test_simple_banded_alignment2(self):
        # event_table = self.dna_handle.get_basecall_data()
        # fastq = self.dna_handle.get_fastq(analysis="Basecall_1D", section='template')
        # print(fastq)
        # nuc_sequence = fastq.split('\n')[1]
        # print(nuc_sequence)
        event_table = np.empty(3, dtype=[('start', float), ('length', float),
                                         ('mean', float), ('stdv', float),
                                         ('model_state', 'S5'), ('move', '<i4'),
                                         ('p_model_state', float)])

        event_table["mean"] = [86.9, 75.3, 80.2]
        nuc_sequence = "AAAAATA"

        events, sum_emission = simple_banded_event_align(event_table, self.dna_model, nuc_sequence)
        self.assertSequenceEqual(events['mean'].tolist(), [86.9, 75.3, 80.2])
        self.assertSequenceEqual([bytes.decode(x) for x in events['model_state']], ["AAAAA", "AAAAT", "AAATA"])
        self.assertSequenceEqual(events['move'].tolist(), [0, 1, 1])

    def test_adaptive_banded_alignment1(self):
        # event_table = self.dna_handle.get_basecall_data()
        # fastq = self.dna_handle.get_fastq(analysis="Basecall_1D", section='template')
        # print(fastq)
        # nuc_sequence = fastq.split('\n')[1]
        # print(nuc_sequence)
        event_table = np.empty(3, dtype=[('start', float), ('length', float),
                                         ('mean', float), ('stdv', float),
                                         ('model_state', 'S5'), ('move', '<i4'),
                                         ('p_model_state', float)])

        event_table["mean"] = [86.9, 75.3, 80.2]
        nuc_sequence = "AAAAATA"

        events, sum_emission = adaptive_banded_simple_event_align(event_table, self.dna_model, nuc_sequence)
        self.assertSequenceEqual(events['mean'].tolist(), [86.9, 75.3, 80.2])
        self.assertSequenceEqual([bytes.decode(x) for x in events['model_state']], ["AAAAA", "AAAAT", "AAATA"])
        self.assertSequenceEqual(events['move'].tolist(), [0, 1, 1])


if __name__ == '__main__':
    unittest.main()
