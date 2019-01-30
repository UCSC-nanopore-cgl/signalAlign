#!/usr/bin/env python
"""Compare multiple hdp and ont trained models"""
########################################################################
# File: compare_trained_models.py
#  executable: compare_trained_models.py
#
# Author: Andrew Bailey
# History: 01/24/18 Created
########################################################################

import os
import numpy as np
from argparse import ArgumentParser
from itertools import zip_longest

from scipy.stats import norm
from sklearn.neighbors import KernelDensity

from py3helpers.utils import load_json, create_dot_dict
from signalalign.hiddenMarkovModel import HmmModel, parse_assignment_file, parse_alignment_file

import matplotlib as mpl

if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
mpl.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--config', '-c', action='store',
                        dest='config', required=True, type=str, default=None,
                        help="Path to json config file")

    args = parser.parse_args()
    return args


def plot_something(models, kmer_list, assignment_data=None, strands=None, savefig_dir=None):
    """Plot multiple kmer distribution onto a single plot with ONT and/or HDP distributions
    :param strands: list of 't' or 'c' for template or complement strands
    :param models: list of HmmModel models to plot
    :param kmer_list: list of kmers for plotting each model
    :param assignment_data: use alignment data if it has already been loaded in
    :param savefig_dir: path to plot save directory
    """
    if savefig_dir:
        assert os.path.exists(savefig_dir), "Save figure directory does not exist: {}".format(savefig_dir)
    assert len(kmer_list) == len(models), \
        "Must have same number of kmer lists: {} and models: {}".format(len(kmer_list), len(models))
    if assignment_data is not None:
        assert strands is not None, "Strand cannot be None if assignment data is passed in."
    # keep track of handles and text depending on which models are loaded
    handles1 = []
    legend_text1 = []
    handles2 = []
    legend_text2 = []
    plt.figure(figsize=(12, 9))
    panel1 = plt.axes([0.1, 0.37, .8, .6])
    panel1.set_xlabel('pA')
    panel1.set_ylabel('Density')
    panel1.grid(color='black', linestyle='-', linewidth=1, alpha=0.5)
    panel1.xaxis.set_major_locator(ticker.AutoLocator())
    panel1.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    min_x = 1000
    max_x = 0
    titles = []
    colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
    markers = ["+", '^', 'o', '1', 's']
    marker_index = 0
    color_index = 0

    if assignment_data is None:
        assignment_data = [None]
    for kmer, model, model_assignment_data, strand in zip_longest(kmer_list, models, assignment_data, strands):
        if kmer is not None:
            nuc_type = "RNA" if model.rna else "DNA"
            strand = "t" if strand is None else strand
            name = "_".join([model.name, nuc_type, strand, kmer])
            normal_mean, normal_sd = model.get_event_mean_gaussian_parameters(kmer)

            tmp_min_x = normal_mean - (5 * normal_sd)
            tmp_max_x = normal_mean + (5 * normal_sd)
            if min_x > tmp_min_x:
                min_x = tmp_min_x
            if max_x < tmp_max_x:
                max_x = tmp_max_x

            # plot ont normal distribution
            x = np.linspace(normal_mean - 4 * normal_sd, normal_mean + 4 * normal_sd, 200)
            ont_handle, = panel1.plot(x, norm.pdf(x, normal_mean, normal_sd), label=kmer, color=colors[color_index])
            color_index += 1
            if color_index > 6:
                color_index = 0
            # panel1.plot([normal_mean, normal_mean], [0, norm.pdf(normal_mean, normal_mean, normal_sd)], lw=2)
            ont_model_name = os.path.basename(model.ont_model_file)
            txt_handle1, = panel1.plot([], [], ' ')
            txt_handle2, = panel1.plot([], [], ' ')
            txt_handle3, = panel1.plot([], [], ' ')

            handles1.append(ont_handle)
            legend_text1.append("{} ONT Normal".format(name))

            handles2.extend([txt_handle1, txt_handle2, txt_handle3])
            print("{} ONT Model: {}".format(name, ont_model_name))
            print("{} ONT Event Mean: {}".format(name, normal_mean))
            print("{} ONT Event SD: {}".format(name, normal_sd))
            legend_text2.extend(["{} ONT Model: {}".format(name, ont_model_name),
                                 "{} ONT Event Mean: {}".format(name, normal_mean),
                                 "{} ONT Event SD: {}".format(name, normal_sd)])

            if model.has_hdp_model:
                # plot HDP predicted distribution
                kmer_id = model.get_kmer_index(kmer)
                x = model.linspace
                if min_x > min(x):
                    min_x = min(x)
                if max_x < max(x):
                    max_x = max(x)
                hdp_y = model.all_posterior_pred[kmer_id]
                if len(hdp_y) == len(x):
                    hdp_handle, = panel1.plot(x, hdp_y, '--', color=colors[color_index])
                    color_index += 1
                    if color_index > 6:
                        color_index = 0

                    handles1.append(hdp_handle)
                    legend_text1.append("{} HDP Distribution".format(name))

            if model_assignment_data is not None:
                kmer_assignments = model_assignment_data.loc[model_assignment_data['kmer'] == kmer]
                kmer_assignments = kmer_assignments.loc[kmer_assignments['strand'] == strand]
                kmer_data = kmer_assignments["level_mean"]
                kmer_prob = kmer_assignments["prob"]
                # get event means and linspace in correct format
                x = np.asarray(kmer_data).reshape(len(kmer_data), 1)
                alphas = np.asarray(kmer_prob).reshape(len(kmer_prob), 1)
                x_plot = model.linspace[:, np.newaxis]
                rgba_colors = np.zeros((len(kmer_data), 4))
                # for red the first column needs to be one
                if 0 < len(titles) < 4:
                    rgba_colors[:, len(titles)] = 1.0

                # the fourth column needs to be your alphas
                rgba_colors[:, 3] = alphas[:, 0]

                # get estimate for data
                if len(kmer_data) > 0:

                    kde = KernelDensity(kernel="gaussian", bandwidth=0.5).fit(x)
                    # estimate across the linspace
                    log_dens = kde.score_samples(x_plot)
                    kde_handle, = panel1.plot(x_plot[:, 0], np.exp(log_dens), '-')
                    raw_data_handle = panel1.scatter(x[:, 0], -0.005 - 0.01 * np.random.random(x.shape[0]),
                                                     marker=markers[marker_index],
                                                     c=rgba_colors)
                    marker_index += 1
                    if marker_index > 4:
                        marker_index = 0

                    # add to legend
                    handles1.extend([kde_handle, raw_data_handle])
                    legend_text1.extend(["Gaussian KDE Estimate: {}".format(name),
                                         "Event Means: {} points\nProb: mu: {}, sd:{}".format(len(kmer_data),
                                                                                              np.mean(alphas[:, 0]),
                                                                                              np.std(alphas[:, 0]))])

                else:
                    print("{} not found in alignment file".format(kmer))

            titles.append(name)
    # create legend
    first_legend = panel1.legend(handles1, legend_text1, bbox_to_anchor=(0.43, -.05))
    ax = plt.gca().add_artist(first_legend)

    panel1.legend(handles2, legend_text2, bbox_to_anchor=(1.1, -.07))

    panel1.set_xlim(min_x, max_x)
    panel1.set_title("Kmer distribution comparisons")

    # option to save figure or just show it
    if savefig_dir:
        base_name = "-".join(titles)
        name = "{}.png".format(base_name)
        out_path = os.path.join(savefig_dir, name)
        plt.savefig(out_path)
    else:
        plt.show()


def main(config=None):
    if config is None:
        args = parse_args()
        # load model files
        assert os.path.exists(args.config), "Config file does not exist: {}".format(args.config)
        config = load_json(args.config)

    args = create_dot_dict(config)
    # load model files
    models = []
    kmer_lists = []
    assignment_data = []
    strands = []
    max_plots = 0
    # create models and grab kmer lists
    for model in args.models:
        models.append(HmmModel(ont_model_file=model.ont_model, hdp_model_file=model.hdp_model, rna=model.rna,
                               name=model.name))
        model_kmer_list = model.kmers
        n_kmers_to_plot = len(model_kmer_list)
        kmer_lists.append(model_kmer_list)
        max_plots = n_kmers_to_plot if n_kmers_to_plot > max_plots else max_plots

        if model.builtAlignment_tsv is not None:
            assert os.path.exists(model.builtAlignment_tsv), \
                "builtAlignment_tsv does not exist: {}".format(model.builtAlignment_tsv)
            # read in both types of data
            try:
                assignment_data.append(parse_assignment_file(model.builtAlignment_tsv))
            except ValueError:
                assignment_data.append(parse_alignment_file(model.builtAlignment_tsv))
        else:
            assignment_data.append(None)
        strands.append(model.strand)

    # Start plotting
    for kmer_list in zip_longest(*kmer_lists):
        plot_something(models, kmer_list, assignment_data=assignment_data,
                       strands=strands, savefig_dir=args.save_fig_dir)


if __name__ == "__main__":
    main()
    raise SystemExit
