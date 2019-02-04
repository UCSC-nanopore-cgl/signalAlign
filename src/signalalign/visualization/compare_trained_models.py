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
import csv
import matplotlib as mpl
if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from argparse import ArgumentParser
from itertools import zip_longest
import itertools
from scipy.stats import norm
from sklearn.neighbors import KernelDensity
from scipy.stats import norm, invgauss, entropy
from scipy.spatial.distance import euclidean

from py3helpers.utils import load_json, create_dot_dict, save_json
from signalalign.hiddenMarkovModel import HmmModel, parse_assignment_file, parse_alignment_file, hellinger2




def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--config', '-c', action='store',
                        dest='config', required=True, type=str, default=None,
                        help="Path to json config file")

    args = parser.parse_args()
    return args


class MultipleModelHandler(object):

    def __init__(self, models, strands, assignment_data=None, savefig_dir=None):
        assert len(models) == len(strands), "Must have strand with each model. models = {} :: strands = {}".format(
            models, strands)
        if savefig_dir is not None:
            assert os.path.isdir(savefig_dir), "savefig_dir must be a directory. {}".format(savefig_dir)
        self.models = models
        self.assignment_data = assignment_data
        if self.assignment_data is None:
            self.assignment_data = [None]
        self.strands = strands
        self.savefig_dir = savefig_dir

    def plot_kmer_distribution(self, kmer_list_list):
        """Plot multiple kmer distribution onto a single plot with ONT and/or HDP distributions
        :param kmer_list_list: list of kmers for plotting each model
        """
        if self.savefig_dir:
            assert os.path.exists(self.savefig_dir), "Save figure directory does not exist: {}".format(self.savefig_dir)
        assert len(kmer_list_list) == len(self.models), \
            "Must have same number of kmer lists: {} and models: {}".format(len(kmer_list_list), len(self.models))
        # keep track of handles and text depending on which models are loaded
        handles1 = []
        legend_text1 = []
        handles2 = []
        legend_text2 = []
        plt.figure(figsize=(20, 9))
        panel1 = plt.axes([0.1, 0.3, .8, .6])
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

        for kmer, model, model_assignment_data, strand in zip_longest(kmer_list_list, self.models,
                                                                      self.assignment_data, self.strands):
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
                                                                                                  np.std(
                                                                                                      alphas[:, 0]))])

                    else:
                        print("{} not found in alignment file".format(kmer))

                titles.append(name)
        # create legend
        first_legend = panel1.legend(handles1, legend_text1, bbox_to_anchor=(0, -0.1), loc='upper left')
        ax = plt.gca().add_artist(first_legend)

        panel1.legend(handles2, legend_text2, bbox_to_anchor=(0.5, -.1), loc='upper left')

        panel1.set_xlim(min_x, max_x)
        panel1.set_title("Kmer distribution comparisons")

        # option to save figure or just show it
        if self.savefig_dir:
            base_name = "-".join(titles)
            name = "{}.png".format(base_name)
            out_path = os.path.join(self.savefig_dir, name)
            plt.savefig(out_path)
        else:
            plt.show()

    def plot_all_model_comparisons(self, write_log_file=True):
        """Plot every comparison between each model"""
        plt.figure(figsize=(10, 8))
        panel1 = plt.axes([0.1, 0.08, .85, .2])
        panel1.set_title("Kullback–Leibler Divergence between distributions", x=0.5, y=1.0)
        panel1.set_xlabel('KL Divergence Distance')
        panel1.set_ylabel('Count')
        panel1.grid(color='black', linestyle='-', linewidth=1, alpha=0.5)

        panel2 = plt.axes([0.1, 0.4, .85, .2])
        panel2.set_title("Hellinger Distance between distributions")
        panel2.set_xlabel('Hellinger Distance')
        panel2.set_ylabel('Count')
        panel2.grid(color='black', linestyle='-', linewidth=1, alpha=0.5)

        panel3 = plt.axes([0.1, 0.72, .85, .2])
        panel3.set_title("abs(Median Delta) between distributions")
        panel3.set_xlabel('abs(Median Delta)')
        panel3.set_ylabel('Count')
        panel3.grid(color='black', linestyle='-', linewidth=1, alpha=0.5)

        all_hellinger_distances = []
        all_kl_divergences = []
        all_median_deltas = []

        for model_pair in itertools.combinations(self.models, 2):
            hellinger_distances, kl_divergences, median_deltas = \
                self.compare_distributions_between_models(model_pair[0], model_pair[1])

            if write_log_file and self.savefig_dir:
                kmers = model_pair[0].sorted_kmer_tuple
                model_names = "{}_{}".format(model_pair[0].name, model_pair[1].name)
                hellinger_outpath = os.path.join(self.savefig_dir,
                                                 "{}_{}".format(model_names, "hellinger_distances.tsv"))
                kl_outpath = os.path.join(self.savefig_dir, "{}_{}".format(model_names, "kl_divergences.tsv"))
                delta_outpath = os.path.join(self.savefig_dir, "{}_{}".format(model_names, "median_deltas.tsv"))
                # write kmer_differences
                self.write_kmer_distribution_comparison_logfile(kmers, kl_divergences, hellinger_distances,
                                                                median_deltas, outfile=hellinger_outpath)

            kl_divergences = [x for x in kl_divergences if x is not None]
            all_hellinger_distances.append(hellinger_distances)
            all_kl_divergences.append(kl_divergences)
            all_median_deltas.append(median_deltas)

        max_hellinger = max([max(x) for x in all_hellinger_distances])
        max_kl = max([max(x) for x in all_kl_divergences])
        max_delta = max([max(x) for x in all_median_deltas])
        panel1_bins = np.linspace(0, max_hellinger, num=30)
        panel2_bins = np.linspace(0, max_kl, num=30)
        panel3_bins = np.linspace(0, max_delta, num=30)

        for i, model_pair in enumerate(itertools.combinations(self.models, 2)):
            panel1.hist(all_kl_divergences[i], bins=panel1_bins,
                        label="KL divergences: {} vs {}".format(model_pair[0].name,
                                                                model_pair[1].name))
            panel2.hist(all_hellinger_distances[i], bins=panel2_bins,
                        label="Hellinger distances: {} vs {}".format(model_pair[0].name,
                                                                     model_pair[1].name))
            panel3.hist(all_median_deltas[i], bins=panel3_bins,
                        label="Median Deltas: {} vs {}".format(model_pair[0].name,
                                                               model_pair[1].name))

            panel1.legend(loc='upper right', fancybox=True, shadow=True)
            panel2.legend(loc='upper right', fancybox=True, shadow=True)
            panel3.legend(loc='upper right', fancybox=True, shadow=True)

        if self.savefig_dir:
            plt.savefig(os.path.join(self.savefig_dir, "model_comparisons.png"))
        else:
            plt.show()

    @staticmethod
    def write_kmer_distribution_comparison_logfile(kmers, kl_divergences, hellinger_distances, median_deltas, outfile):
        """Write a sorted by divergence tsv of kmers"""
        assert len(kmers) == len(kl_divergences), \
            "Number of kmers and divergences must match. " \
            "n_kmers : {} != n_divergences: {}".format(len(kmers), len(kl_divergences))
        assert len(kmers) == len(hellinger_distances), \
            "Number of kmers and hellinger_distances must match. n_kmers : " \
            "{} != n_hellinger_distances: {}".format(len(kmers), len(hellinger_distances))
        assert len(kmers) == len(median_deltas), \
            "Number of kmers and median_deltas must match. " \
            "n_kmers : {} != n_median_deltas: {}".format(len(kmers), len(median_deltas))

        zipped_kmers = [(k, d1, d2, d3) for k, d1, d2, d3 in
                        zip(kmers, kl_divergences, hellinger_distances, median_deltas)
                        if d1 is not None]

        zipped_kmers.sort(key=lambda x: x[1], reverse=True)
        none_zipped_kmers = [(k, d1, d2, d3) for k, d1, d2, d3 in
                             zip(kmers, kl_divergences, hellinger_distances, median_deltas)
                             if d1 is None]

        with open(outfile, 'w') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t')
            writer.writerows(zipped_kmers)
            writer.writerows(none_zipped_kmers)

        return outfile

    @staticmethod
    def read_kmer_distribution_comparison_logfile(infile):
        """Read in kmer distribution comparison tsv logfile"""
        data = []
        with open(infile, newline='\n') as csvfile:
            spamreader = csv.reader(csvfile, delimiter='\t')
            for row in spamreader:
                # catch None's in tsv
                d1 = None if row[1] == '' else float(row[1])
                d2 = None if row[2] == '' else float(row[2])
                d3 = None if row[3] == '' else float(row[3])
                data.append([row[0], d1, d2, d3])
        return data

    def compare_distributions_between_models(self, model1, model2, hdp=True):
        """Calculate hellinger divergence and kl divergence between the hdp or hmm model between two models."""
        hellinger_distances = []
        kl_divergences = []
        median_deltas = []
        get_new_linspace = False
        if model1.has_hdp_model and model2.has_hdp_model and hdp:
            if np.array_equal(model1.linspace, model2.linspace):
                linspace = model1.linspace
            else:
                get_new_linspace = True
                linspace_min = max([model1.linspace[0], model2.linspace[0]])
                linspace_max = min([model1.linspace[-1], model2.linspace[-1]])
                linspace = np.linspace(linspace_min, linspace_max, 3000)
        elif model1.has_hdp_model:
            linspace = model1.linspace
        else:
            linspace = model2.linspace

        for kmer in model1.sorted_kmer_tuple:
            # if statements used if the HDP model does not have information on the kmer distribution
            if hdp and model1.has_hdp_model:
                m1_dist = self.get_hdp_kmer_posterior_prediction(model1, kmer, linspace, get_new_linspace)
                if m1_dist is None:
                    m1_dist = self.get_ont_kmer_posterior_prediction(model1, kmer, linspace)
            else:
                m1_dist = self.get_ont_kmer_posterior_prediction(model1, kmer, linspace)

            if hdp and model2.has_hdp_model:
                m2_dist = self.get_hdp_kmer_posterior_prediction(model2, kmer, linspace, get_new_linspace)
                if m2_dist is None:
                    m2_dist = self.get_ont_kmer_posterior_prediction(model2, kmer, linspace)
            else:
                m2_dist = self.get_ont_kmer_posterior_prediction(model2, kmer, linspace)

            kl_divergences.append(self.get_kl_divergence(m1_dist, m2_dist))
            hellinger_distances.append(self.get_hellinger_distance(m1_dist, m2_dist))
            median_deltas.append(self.get_median_delta(m1_dist, m2_dist, linspace))

        return hellinger_distances, kl_divergences, median_deltas

    @staticmethod
    def get_hdp_kmer_posterior_prediction(model, kmer, linspace, get_new_linspace=False):
        """For a given model, grab the posterior prediction distribution"""
        if model.has_hdp_model:
            if get_new_linspace:
                posterior_pred = model.get_new_linspace_hdp_probability_distribution(kmer, linspace)
            else:
                kmer_id = model.get_kmer_index(kmer)
                posterior_pred = model.all_posterior_pred[kmer_id]
                # print("[Kullback–Leibler divergence] No HDP data for {}".format(kmer))
            if posterior_pred is None:
                return None
            elif len(posterior_pred) == 0:
                return None
            return posterior_pred
        else:
            return None

    @staticmethod
    def get_ont_kmer_posterior_prediction(model, kmer, linspace):
        """For a given model, grab the posterior prediction distribution"""
        # print("[Kullback–Leibler divergence] No HDP data for {}".format(kmer))
        normal_mean, normal_sd = model.get_event_mean_gaussian_parameters(kmer)
        posterior_pred = norm.pdf(linspace, normal_mean, normal_sd)

        return posterior_pred

    @staticmethod
    def get_kl_divergence(dist1, dist2):
        """Get Kullback–Leibler divergence between the HDP and ONT models for a specific kmer"""
        if min(dist1) == 0:
            dist1[dist1 == 0] = 0.000001
        #     np.nextafter(0, 1)
        if min(dist2) == 0:
            dist2[dist2 == 0] = 0.000001

        kl_divergence = entropy(pk=dist1, qk=dist2, base=2)
        if kl_divergence == np.inf:
            # print("[Kullback–Leibler divergence] Zero probability for {}".format(kmer))
            return None
        return kl_divergence

    @staticmethod
    def get_hellinger_distance(dist1, dist2):
        """Get Hellinger distance between the HDP and ONT models for a specific kmer"""
        h_distance = hellinger2(p=dist1, q=dist2)
        return h_distance

    @staticmethod
    def get_median_delta(dist1, dist2, linspace):
        """Calculate the difference between the max value of HDP and ONT kmer distributions"""
        dist1 = list(dist1)
        dist2 = list(dist2)
        delta = linspace[dist1.index(max(dist1))] - linspace[dist2.index(max(dist2))]
        return abs(delta)


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

    mmh = MultipleModelHandler(models, strands=strands, assignment_data=assignment_data, savefig_dir=args.save_fig_dir)
    if args.summary_distance:
        mmh.plot_all_model_comparisons()
    # Start plotting
    for kmer_list in zip_longest(*kmer_lists):
        mmh.plot_kmer_distribution(kmer_list)

    if args.save_fig_dir:
        save_json(args, os.path.join(args.save_fig_dir, "compare_trained_models_config.json"))


if __name__ == "__main__":
    main()
    raise SystemExit
