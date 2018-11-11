#!/usr/bin/env python
"""Plot kmer distributions given ONT model file with options for the HDP model file and the buildAssignments.tsv file"""
########################################################################
# File: plot_kmer_distributions.py
#  executable: plot_kmer_distributions.py
#
# Author: Andrew Bailey
# History: 08/10/18 Created
########################################################################

from signalalign.hiddenMarkovModel import HmmModel, parse_alignment_file
from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--ont_model', '-m', action='store',
                        dest='ont_model', required=True, type=str, default=None,
                        help="Path to ONT model file")

    parser.add_argument('--hdp_model', action='store', default=None,
                        dest='hdp_model', required=False, type=str,
                        help="HDP model file")

    parser.add_argument('--build_alignment_file', action='store', default=None,
                        dest='build_alignment_file', required=False, type=str,
                        help="TSV alignment file which is where all the data for the HDP is stored")

    parser.add_argument('--output_dir', action='store',
                        dest='output_dir', required=False, type=str,
                        help="If set, will write out to output directory otherwise the graph is just shown")

    parser.add_argument('--rna', action='store_true', default=False,
                        dest='rna', required=False,
                        help="Boolean option for writing RNA onto graphs")

    parser.add_argument('--kmer', action='store', default=None,
                        dest='kmer', required=False,
                        help="Kmer you want to plot if you only want to plot a single kmer")

    parser.add_argument('--all_kmers', action='store_true', default=False,
                        dest='all_kmers', required=False,
                        help="Iterate all possible kmers and plot the distributions.")

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    # load model files
    model_handle = HmmModel(args.ont_model, args.hdp_model, rna=args.rna)

    # output a single plot or one for each kmer
    if args.kmer:
        model_handle.plot_kmer_distribution(args.kmer,
                                            alignment_file=args.build_alignment_file,
                                            savefig_dir=args.output_dir)
    else:
        assert args.all_kmers, "Must pick a single kmer to plot using --kmer or pass the flag --all_kmers"
        alignment_data = None
        if args.build_alignment_file:
            alignment_data = parse_alignment_file(args.build_alignment_file)
        for kmer in model_handle.sorted_kmer_tuple:
            model_handle.plot_kmer_distribution(kmer,
                                                alignment_file_data=alignment_data,
                                                savefig_dir=args.output_dir)


if __name__ == "__main__":
    main()
    raise SystemExit
