#!/usr/bin/env python
"""Plot kmer distributions given ONT model file with options for the HDP model file and the buildAssignments.tsv file"""
########################################################################
# File: plot_kmer_distributions.py
#  executable: plot_kmer_distributions.py
#
# Author: Andrew Bailey
# History: 08/10/18 Created
########################################################################

import os
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from signalalign.hiddenMarkovModel import HmmModel


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--rna_ont_model', action='store',
                        dest='rna_ont_model', required=True, type=str, default=None,
                        help="Path to RNA ONT model file")

    parser.add_argument('--rna_hdp_model', action='store', default=None,
                        dest='rna_hdp_model', required=True, type=str,
                        help="Path to RNA HDP model file")

    parser.add_argument('--dna_ont_model', action='store',
                        dest='dna_ont_model', required=True, type=str, default=None,
                        help="Path to DNA ONT model file")

    parser.add_argument('--dna_hdp_model', action='store', default=None,
                        dest='dna_hdp_model', required=True, type=str,
                        help="Path to DNA HDP model file")

    parser.add_argument('--output_dir', action='store',
                        dest='output_dir', required=False, type=str,
                        help="If set, will write out to output directory otherwise the graph is just shown")

    args = parser.parse_args()
    return args


def plot_rna_dna_ont_hdp_comparison(rna_ont_model, rna_hdp_model, dna_ont_model, dna_hdp_model, savefig_dir=None):
    """Plot both the RNA and DNA distribution of Hellinger distances and Kullback–Leibler Divergences between
    the default ONT model and the corresponding HDP model file
    """
    hmm_handle = HmmModel(rna_ont_model, rna_hdp_model)
    # means = hmm_handle.event_model['means']
    dna_hmm_handle = HmmModel(dna_ont_model, dna_hdp_model)

    hellinger_distances, kl_divergences, median_deltas = hmm_handle.compare_distributions()
    dna_hellinger_distances, dna_kl_divergences, dna_median_deltas = dna_hmm_handle.compare_distributions()
    # print(dna_hellinger_distances)
    # print(dna_kl_divergences)
    # hellinger distances plot comparisons
    plt.figure(figsize=(8, 8))
    panel1 = plt.axes([0.05, 0.08, .9, .2])
    panel1.set_title("DNA and RNA Hellinger Distance between\nHDP distributions and ONT normal distributions", x=0.5, y=1.0)
    panel1.set_xlabel('Hellinger Distance')
    panel1.set_ylabel('Count')
    panel1.grid(color='black', linestyle='-', linewidth=1, alpha=0.5)
    # panel1.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
    panel1_bins = np.linspace(0, max(hellinger_distances+dna_hellinger_distances), num=30)
    panel1.hist(hellinger_distances, bins=panel1_bins, label="RNA Hellinger Distances")
    panel1.hist(dna_hellinger_distances, bins=panel1_bins, alpha=0.7, label="DNA Hellinger Distances")
    panel1.legend(loc='upper right', fancybox=True, shadow=True)

    # kl divergence plot comparisons
    panel2 = plt.axes([0.05, 0.4, .9, .2])
    panel2.set_title("DNA and RNA Kullback–Leibler Divergence between\nHDP distributions and ONT normal distributions")

    panel2.set_xlabel('Kullback–Leibler Divergence')
    panel2.set_ylabel('Count')
    panel2.grid(color='black', linestyle='-', linewidth=1, alpha=0.5)
    # panel2.xaxis.set_major_locator(ticker.MultipleLocator(1))
    panel2_bins = np.linspace(0, max(dna_kl_divergences+kl_divergences), num=30)

    panel2.hist(kl_divergences, bins=panel2_bins, label="RNA KL Divergence")
    panel2.hist(dna_kl_divergences, bins=panel2_bins, alpha=0.7, label="DNA KL Divergence")
    panel2.legend(loc='upper right', fancybox=True, shadow=True)

    panel3 = plt.axes([0.05, 0.72, .9, .2])
    panel3.set_title("DNA and RNA abs(Median Delta) between\nHDP distributions and ONT normal distributions")

    panel3.set_xlabel('abs(Median Delta)')
    panel3.set_ylabel('Count')
    panel3.grid(color='black', linestyle='-', linewidth=1, alpha=0.5)
    # panel2.xaxis.set_major_locator(ticker.MultipleLocator(1))
    panel3_bins = np.linspace(0, max(median_deltas+dna_median_deltas), num=30)

    panel3.hist(median_deltas, bins=panel3_bins, label="RNA Median Delta")
    panel3.hist(dna_median_deltas, bins=panel3_bins, alpha=0.7, label="DNA Median Delta")
    panel3.legend(loc='upper right', fancybox=True, shadow=True)

    if savefig_dir:
        plt.savefig(os.path.join(savefig_dir, "dna_rna_kl_hellinger_delta.png"))
    else:
        plt.show()


def main():
    args = parse_args()

    plot_rna_dna_ont_hdp_comparison(args.rna_ont_model, args.rna_hdp_model,
                                    args.dna_ont_model, args.dna_hdp_model,
                                    savefig_dir=args.output_dir)


if __name__ == "__main__":
    main()
    raise SystemExit
