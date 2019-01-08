#!/usr/bin/env python
"""Utility functions and classes for creating mixture models"""
########################################################################
# File: mixture_model.py
#
# Author: Andrew Bailey
# History: 1/7/19 Created
########################################################################


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import os
import tempfile
from sklearn.mixture import GaussianMixture
from sklearn.datasets.samples_generator import make_blobs
from signalalign.hiddenMarkovModel import HmmModel, create_new_model
from signalalign.signalAlignment import SignalAlignment
from signalalign.train.trainModels import make_master_assignment_table


def get_nanopore_gauss_mixture(event_means, n_models):
    """Create a mixture model from an array of event means and fit some number of models to the data
    :param event_means: array of event means
    :param n_models: number of gaussians to fit
    """
    model = GaussianMixture(n_models).fit(event_means)
    assert model.converged_, "Model has not converged"
    return model


def find_best_1d_gaussian_fit(x, max_n, aic=True):
    """
    :param x: input data
    :param max_n: max number of gaussians to try to fit to data
    :param aic: boolean option to use aic or bic. if false use bic as selection criterion
    :return:
    """
    N = np.arange(1, max_n)
    models = [None for i in range(len(N))]

    for i in range(len(N)):
        models[i] = GaussianMixture(N[i]).fit(x)

    # use AIC or BIC for model selection
    if aic:
        aic = [m.aic(x) for m in models]
        m_best = models[np.argmin(aic)]
    else:
        bic = [m.bic(x) for m in models]
        m_best = models[np.argmin(bic)]

    return m_best


# TODO check if this is right! I dont think using the weights like this works
def get_mus_and_sigmas_1d(gaussian_model):
    """Get the mean and stdv of each normal curve given a GaussianMixture model
    :param gaussian_model: an already converged GaussianMixture model
    :return: list of tuples with tup[0] = mu and tup[1] = sigma
    """
    assert gaussian_model.converged_, "Model has not converged"
    normals = []
    for i, mu in enumerate(gaussian_model.means_):
        assert len(gaussian_model.covariances_[i]) == 1, "This function only works for 1D gaussian mixture models"
        sigma = np.sqrt(gaussian_model.covariances_[i][0])
        # sigma = sigma / gaussian_model.weights_[i]
        normals.append((mu, sigma))

    return normals


def closest_to_canonical(mixture_normals, canonical_mu):
    """Find the normal distribution closet to canonical mu"""
    min_index = 0
    min_distance = 1000
    for i in range(len(mixture_normals)):
        mu = mixture_normals[i][0]
        distance = abs(mu - canonical_mu)
        if distance < min_distance:
            min_index = i
            min_distance = distance

    match = mixture_normals.pop(min_index)
    return match, mixture_normals


def fit_model_to_kmer_dist(all_assignments, kmer, n_normals=2):
    """Return a mixture model from the distribution of event means for a given kmer
    :param all_assignments: master table of assignments (must have fields "k-mer" and "descaled_event_mean"
    :param kmer: str that must be in the assignments table
    :param n_normals: number of normal gaussians to fit to distirbution
    """
    kmer = str.encode(kmer)
    samples = np.array(all_assignments[all_assignments["k-mer"] == kmer]["descaled_event_mean"]).reshape(-1, 1)
    model = get_nanopore_gauss_mixture(samples, n_normals)
    return model


def main():
    base_model = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/models/testModelR9_acgt_template.model"
    # assignment_dir = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/tests/test_assignment_files"
    # alphabet = "ATGCJL"
    built_model = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/models/buildAlignment.tsv"
    model_h = HmmModel(base_model)
    assignments = SignalAlignment.read_in_signal_align_tsv(built_model, "assignments")
    # assignments = make_master_assignment_table(assignment_dir)
    # X = assignments['descaled_event_mean'].reshape((-1, 1))
    main_kmer = "TAATT"
    new_kmer = "TAJTT"

    # print(assignments[assignments["k-mer"] == main_kmer]["descaled_event_mean"])



    fit_model_to_kmer_dist(assignments, main_kmer, n_normals=2)
    with tempfile.TemporaryDirectory() as tempdir:
        test_model_file = os.path.join(tempdir, "fake.hmm")
        new_model = create_new_model(base_model, test_model_file, (("A", "J")))
    new_model.plot_kmer_distribution(main_kmer)
    new_model.plot_kmer_distribution(new_kmer)


    ##################################################
    ##################################################
    ##################################################

    # X = np.array(assignments[assignments["k-mer"] == b'TAATT']["descaled_event_mean"]).reshape(-1, 1)
    #
    # min_range = 10
    # max_range = 140
    # original_mu = 100
    # # ------------------------------------------------------------
    # # X, y_true = make_blobs(n_samples=400, n_features=1, centers=4,
    # #                        cluster_std=4, random_state=0, center_box=(min_range, max_range))
    # # ------------------------------------------------------------
    # # Learn the best-fit GMM models
    # #  Here we'll use GMM in the standard way: the fit() method
    # #  uses an Expectation-Maximization approach to find the best
    # #  mixture of Gaussians for the data
    #
    # # model = find_best_1d_gaussian_fit(X, 3)
    # N = np.arange(1, 8)
    # models = [None for i in range(len(N))]
    #
    # for i in range(len(N)):
    #     models[i] = GaussianMixture(N[i]).fit(X)
    #
    # # use AIC or BIC for model selection
    # aic = [m.aic(X) for m in models]
    # bic = [m.bic(X) for m in models]
    #
    # fig = plt.figure(figsize=(10, 3))
    # fig.subplots_adjust(left=0.12, right=0.97,
    #                     bottom=0.21, top=0.9, wspace=0.5)
    #
    # # plot 1: data + best-fit mixture
    # ax = fig.add_subplot(131)
    # M_best = models[np.argmin(aic)]
    # print(M_best.means_)
    # print(M_best.covariances_)
    # print(M_best.weights_)
    # x = np.linspace(min_range, max_range, 1000).reshape(1000, 1)
    # responsibilities = M_best.predict_proba(x)
    # logprob = M_best.score_samples(x)
    #
    # pdf = np.exp(logprob)
    # pdf_individual = responsibilities * pdf[:, np.newaxis]
    #
    # ax.hist(X, 30, normed=True, histtype='stepfilled', alpha=0.4)
    # ax.plot(x, pdf, '-k', label="mixture pdf")
    # ax.plot(x, pdf_individual, '--k', label="individual pdf")
    # ax.set_title("Best-fit Mixture")
    # ax.set_xlabel('$x$')
    # ax.set_ylabel('$p(x)$')
    #
    # mixture_normals = get_mus_and_sigmas_1d(M_best)
    # closest, the_rest = closest_to_canonical(mixture_normals, original_mu)
    # # plot the closest to original
    # x = np.linspace(closest[0] - 3 * closest[1], closest[0] + 3 * closest[1], 100)
    # ax.plot(x, mlab.normpdf(x, closest[0], closest[1]), label="Closest to Original")
    #
    # for mu, sigma in the_rest:
    #     x = np.linspace(mu - 3 * sigma, mu + 3 * sigma, 100)
    #     ax.plot(x, mlab.normpdf(x, mu, sigma), label="New_models")
    # ax.legend(loc=2)
    #
    # # plot 2: AIC and BIC
    # ax = fig.add_subplot(132)
    # ax.plot(N, aic, '-k', label='AIC')
    # ax.plot(N, bic, '--k', label='BIC')
    # ax.set_xlabel('n. components')
    # ax.set_ylabel('information criterion')
    # ax.legend(loc=2)
    #
    # plt.show()


if __name__ == '__main__':
    main()
