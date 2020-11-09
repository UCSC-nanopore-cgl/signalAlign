#!/usr/bin/env python
"""Utility functions and classes for creating mixture models"""
########################################################################
# File: mixture_model.py
#
# Author: Andrew Bailey
# History: 1/7/19 Created
########################################################################


import os
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser
from sklearn.mixture import GaussianMixture
from timeit import default_timer as timer
from py3helpers.utils import load_json, create_dot_dict
from scipy.stats import norm
from sklearn.neighbors import KernelDensity

from signalalign.hiddenMarkovModel import HmmModel, parse_assignment_file, parse_alignment_file
from signalalign.utils.sequenceTools import get_motif_kmers, find_modification_index_and_character

import tempfile
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from sklearn.datasets.samples_generator import make_blobs


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--config', '-c', required=True, action='store',
                        dest='config', type=str, default=None,
                        help="Path to config file")

    args = parser.parse_args()
    return args


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
    return match, mixture_normals, min_distance


def fit_model_to_kmer_dist(all_assignments, kmer, n_normals=2):
    """Return a mixture model from the distribution of event means for a given kmer
    :param all_assignments: master table of assignments (must have fields "k-mer" and "descaled_event_mean"
    :param kmer: str that must be in the assignments table
    :param n_normals: number of normal gaussians to fit to distirbution
    """
    samples = all_assignments[all_assignments["kmer"] == kmer]["level_mean"].values.reshape(-1, 1)
    model = False
    if len(samples) == 0:
        print("No alignments found for kmer: {}".format(kmer))
    else:
        model = get_nanopore_gauss_mixture(samples, n_normals)
    return model


def generate_gaussian_mixture_model_for_motifs(model_h, assignments, all_kmer_pars, strand,
                                               output_dir, plot=False, name="", target_model=None, show=False,
                                               close=True):
    """Generate new hmm model using mixture model of assignment data for each required kmer given the set of motifs
    :param model_h: HmmModel
    :param strand: 't' for template or 'c' for complement
    :param plot: plot model data
    :param assignments: assignment DataFrame with "strand", "kmer" and "level_mean"
    :param all_kmer_pars: list of list of [canonical, modified] kmers
    :param output_dir: path to save figures, models and log file
    :param name: optional argument for naming the mixture model
    :param target_model: use for plotting expected distribution for modified kmer
    """
    assert strand in ('t', 'c'), "Strand must be either 'c' or 't'. strand = {}".format(strand)
    assignments = assignments[assignments["strand"] == strand]
    canonical_mixture_components_comparison = []
    if name is not "":
        name += "_"
    output_model_path = os.path.join(output_dir, "{}_{}mixture_model.hmm".format(strand, name))

    for kmer_pair in all_kmer_pars:
        old_kmer = kmer_pair[0]
        new_kmer = kmer_pair[1]

        # fit
        mixture_model = fit_model_to_kmer_dist(assignments, old_kmer, n_normals=2)
        if mixture_model:
            mixture_normals = get_mus_and_sigmas_1d(mixture_model)
            kmer_mean, kmer_sd = model_h.get_event_mean_gaussian_parameters(old_kmer)
            match, other, distance = closest_to_canonical(mixture_normals, kmer_mean)
            # set parameters
            model_h.set_kmer_event_mean_params(new_kmer, other[0][0][0], other[0][1][0])
            canonical_mixture_components_comparison.append(
                [old_kmer, kmer_mean, kmer_sd, match[0][0], match[1][0], other[0][0][0],
                 other[0][1][0], distance[0], strand])
            print(old_kmer, mixture_normals)
            if plot:
                # model_h.plot_kmer_distributions([old_kmer, new_kmer],
                #                                 alignment_file_data=assignments,
                #                                 savefig_dir=output_dir,
                #                                 name=strand)
                plot_output_dir = output_dir
                if show:
                    plot_output_dir = None
                plot_mixture_model_distribution(old_kmer, new_kmer, kmer_mean, kmer_sd, match[0][0], match[1][0],
                                                other[0][0][0], other[0][1][0], strand, mixture_model=mixture_model,
                                                kmer_assignments=assignments,
                                                save_fig_dir=plot_output_dir,
                                                target_model=target_model, close=close)
    model_h.normalize(False, False)
    model_h.write(output_model_path)

    data = pd.DataFrame(canonical_mixture_components_comparison, columns=["kmer",
                                                                          "canonical_model_mean",
                                                                          "canonical_model_sd",
                                                                          "canonical_mixture_mean",
                                                                          "canonical_mixture_sd",
                                                                          "modified_mixture_mean",
                                                                          "modified_mixture_sd",
                                                                          "distance",
                                                                          "strand"])
    data.sort_values("distance", inplace=True, ascending=False)
    log_file = os.path.join(output_dir, "{}_distances.tsv".format(strand))
    data.to_csv(log_file, sep="\t", index=False)
    return data


def get_motif_kmer_pairs(motif_pair, k, alphabet="ATGC"):
    """Given a motif pair, create a list of all kmers which contain modification """
    all_kmer_pars = []
    motif_kmers = get_motif_kmers(motif_pair, k, alphabet=alphabet)
    pos, old_char, new_char = find_modification_index_and_character(motif_pair[0], motif_pair[1])

    for new_kmer in motif_kmers:
        # get original kmer
        pos = new_kmer.find(new_char)
        old_kmer = new_kmer[0:pos] + old_char + new_kmer[pos + 1:]
        all_kmer_pars.append([old_kmer, new_kmer])
    return all_kmer_pars


def plot_mixture_model_distribution(canonical_kmer, modified_kmer, canonical_model_mean, canonical_model_sd,
                                    canonical_mixture_mean,
                                    canonical_mixture_sd, modified_mixture_mean, modified_mixture_sd,
                                    strand, mixture_model=None, target_model=None,
                                    kmer_assignments=None, save_fig_dir=None, close=False):
    """Plot normal distributions from mixture model and compare with original canonical model
    :param canonical_model_mean: canonical_model_mean
    :param canonical_model_sd: canonical_model_sd
    :param canonical_mixture_mean: canonical_mixture_mean
    :param canonical_mixture_sd: canonical_mixture_sd
    :param modified_mixture_mean: modified_mixture_mean
    :param modified_mixture_sd: modified_mixture_sd
    :param save_fig_dir: optional path to save figure
    :param strand: template or complement ('t' or 'c')
    :param canonical_kmer: kmer to plot
    :param modified_kmer: modified kmer
    :param target_model: model to compare the mixture to
    :param mixture_model: an already fit GaussianMixture model
    :param kmer_assignments: assignments with ("level_mean" and "kmer") named columns of DataFrame
    """
    fig = plt.figure(figsize=(12, 8))
    panel1 = plt.axes([0.1, 0.1, .6, .8])
    panel1.set_xlabel('pA')
    panel1.set_ylabel('Density')
    panel1.grid(color='black', linestyle='-', linewidth=1, alpha=0.5)
    panel1.set_title("Mixture Model Comparison: {}".format(canonical_kmer))

    # original canonical model
    x = np.linspace(canonical_model_mean - 4 * canonical_model_sd, canonical_model_mean + 4 * canonical_model_sd, 200)
    panel1.plot(x, norm.pdf(x, canonical_model_mean, canonical_model_sd), label="{} ONT model".format(canonical_kmer))

    # selected mixture model for canonical kmer
    x = np.linspace(canonical_mixture_mean - 4 * canonical_mixture_sd,
                    canonical_mixture_mean + 4 * canonical_mixture_sd, 200)
    panel1.plot(x, norm.pdf(x, canonical_mixture_mean, canonical_mixture_sd), label="{} mixture".format(canonical_kmer))

    # selected mixture model for modified kmer
    x = np.linspace(modified_mixture_mean - 4 * modified_mixture_sd,
                    modified_mixture_mean + 4 * modified_mixture_sd, 200)
    panel1.plot(x, norm.pdf(x, modified_mixture_mean, modified_mixture_sd), label="{} mixture".format(modified_kmer))
    x_min = min([canonical_mixture_mean - 4 * canonical_mixture_sd, modified_mixture_mean - 4 * modified_mixture_sd,
                 canonical_model_mean - 4 * canonical_model_sd])
    x_max = max([canonical_mixture_mean + 4 * canonical_mixture_sd, modified_mixture_mean + 4 * modified_mixture_sd,
                 canonical_model_mean + 4 * canonical_model_sd])
    panel1.set_xlim(x_min, x_max)

    if mixture_model is not None:
        x = np.linspace(x_min, x_max, 1000).reshape(1000, 1)
        responsibilities = mixture_model.predict_proba(x)
        logprob = mixture_model.score_samples(x)

        pdf = np.exp(logprob)
        pdf_individual = responsibilities * pdf[:, np.newaxis]

        panel1.plot(x, pdf, '-k', label="mixture pdf")
        panel1.plot(x, pdf_individual, '--k', label="individual pdf")

    if kmer_assignments is not None:
        kmer_assignments = kmer_assignments.loc[kmer_assignments['kmer'] == canonical_kmer]
        kmer_data = kmer_assignments["level_mean"]
        # get event means and linspace in correct format
        x = np.asarray(kmer_data).reshape(len(kmer_data), 1)
        x_plot = np.linspace(x_min, x_max, 200)[:, np.newaxis]
        # get estimate for data
        if len(kmer_data) > 0:

            kde = KernelDensity(kernel="gaussian", bandwidth=0.5).fit(x)
            # estimate across the linspace
            log_dens = kde.score_samples(x_plot)
            panel1.plot(x_plot[:, 0], np.exp(log_dens), '-', label="Gaussian KDE Estimate")
            panel1.plot(x[:, 0], -0.005 - 0.01 * np.random.random(x.shape[0]), '+k',
                        label="Event Means: {} points".format(len(kmer_data)))

    if target_model is not None:
        if target_model.has_hdp_model:
            # plot HDP predicted distribution
            kmer_id = target_model.get_kmer_index(modified_kmer)
            x = target_model.linspace
            hdp_y = target_model.all_posterior_pred[kmer_id]
            if len(hdp_y) == len(x):
                panel1.plot(x, hdp_y, '-', label="{} HDP Distribution".format(modified_kmer))
        else:
            normal_mean, normal_sd = target_model.get_event_mean_gaussian_parameters(modified_kmer)
            x = np.linspace(normal_mean - 4 * normal_sd, normal_mean + 4 * normal_sd, 200)
            panel1.plot(x, norm.pdf(x, normal_mean, normal_sd), label="{} ONT Distribution".format(modified_kmer))

    panel1.legend(loc='upper right', fancybox=True, shadow=True)
    # option to save figure or just show it
    if save_fig_dir:
        out_name = "{}_{}_{}_{}.png".format(canonical_kmer, modified_kmer, strand, "mixture_model")
        out_path = os.path.join(save_fig_dir, out_name)
        plt.savefig(out_path)
    else:
        plt.show()
    if close:
        plt.close(fig)


def main(config=None):
    """Plot event to reference labelled ONT nanopore reads"""
    start = timer()
    if config is None:
        args = parse_args()
        # load model files
        assert os.path.exists(args.config), "Config file does not exist: {}".format(args.config)
        config = load_json(args.config)

    args = create_dot_dict(config)
    # get assignments and load model
    try:
        assignments = parse_assignment_file(args.assignments)
    except ValueError:
        assignments = parse_alignment_file(args.assignments)

    model_h = HmmModel(args.model_path, rna=args.rna)
    target_model = None
    if args.target_hmm_model is not None:
        target_model = HmmModel(args.target_hmm_model, hdp_model_file=args.target_hdp_model, rna=args.rna)
    # generate kmers to match
    all_kmer_pairs = set()
    for motif in args.motifs:
        all_kmer_pairs |= set(tuple(row) for row in get_motif_kmer_pairs(motif_pair=motif, k=model_h.kmer_length))

    data = generate_gaussian_mixture_model_for_motifs(model_h, assignments, all_kmer_pairs, args.strand,
                                                      args.output_dir, plot=args.plot, name="ccwgg",
                                                      target_model=target_model, show=args.show)
    # data = pd.read_csv(os.path.join(args.output_dir, "t_distances.tsv"), delimiter="\t")
    # data = data.ix[0]
    # plot_mixture_model_distribution(data["kmer"], data["canonical_model_mean"], data["canonical_model_sd"],
    #                                 data["canonical_mixture_mean"],
    #                                 data["canonical_mixture_sd"], data["modified_mixture_mean"],
    #                                 data["modified_mixture_sd"],
    #                                 data["strand"], kmer_assignments=assignments, save_fig_dir=None)
    stop = timer()
    print("Running Time = {} seconds".format(stop - start), file=sys.stderr)

    ##################################################
    ##################################################
    ##################################################
    # t_base_model = "/Users/andrewbailey/data/ccwgg_new_em_trained_model/all_models/template_trained.hmm"
    # c_base_model = "/Users/andrewbailey/data/ccwgg_new_em_trained_model/all_models/complement_trained.hmm"
    # built_model = "/Users/andrewbailey/data/ccwgg_new_em_trained_model/buildAlignment_no_E2.tsv"
    # save_fig_dir = "/Users/andrewbailey/data/ccwgg_new_em_trained_model/mixture_models"
    # assignments = parse_assignment_file(built_model)
    # motifs = [["CCAGG", "CEAGG"], ["CCTGG", "CETGG"]]
    # generate_gaussian_mixture_model_for_motifs(t_base_model, assignments, motifs, "t",
    #                                            save_fig_dir, rna=rna, plot=plot, name="ccwgg")
    # generate_gaussian_mixture_model_for_motifs(c_base_model, assignments, motifs, "c",
    #                                            save_fig_dir, rna=rna, plot=plot, name="ccwgg")

    # X = assignments[assignments["kmer"] == main_kmer]["level_mean"].values.reshape(-1, 1)
    #
    # min_range = 60
    # max_range = 140
    # original_mu = new_model_h.get_event_mean_gaussian_parameters(main_kmer)[0]
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
    # ax = fig.add_subplot(121)
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
    # ax.set_title("Best-fit Mixture: {}".format(main_kmer))
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
    # ax = fig.add_subplot(122)
    # ax.plot(N, aic, '-k', label='AIC')
    # ax.plot(N, bic, '--k', label='BIC')
    # ax.set_xlabel('n. components')
    # ax.set_ylabel('information criterion')
    # ax.legend(loc=2)
    #
    # plt.show()


if __name__ == '__main__':
    main()
