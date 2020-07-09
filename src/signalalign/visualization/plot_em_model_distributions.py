#!/usr/bin/env python
"""Plot distributions for positions covering modified bases"""
########################################################################
# File: plot_em_model_distributions.py
#  executable: plot_em_model_distributions.py
#
# Author: Andrew Bailey
# History: 03/16/20 Created
########################################################################

from argparse import ArgumentParser
from py3helpers.utils import load_json, create_dot_dict, time_it
from py3helpers.seq_tools import ReferenceHandler
from signalalign.utils.sequenceTools import CustomAmbiguityPositions, reverse_complement
from signalalign.hiddenMarkovModel import HmmModel
from signalalign.visualization.compare_trained_models import MultipleModelHandler
import os
import shutil


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--config', '-c', action='store',
                        dest='config', required=True, type=str, default=None,
                        help="Path to config file")
    args = parser.parse_args()
    return args


def get_covered_kmers(reference_handler: ReferenceHandler, chromosome_name: str, strand: str, pos: list,
                      variant_chars: list,
                      rna: bool = True, kmer_length: int = 5) -> list:
    """Get a list of lists of kmers covering a certain position. The first list will always be canonical,
    any subsequent kmer lists are the same kmers with the specified position replaced with non
    reference characters
    :param reference_handler: ReferenceHandler object
    :param chromosome_name: reference chromosome name
    :param strand: strand
    :param pos: 0 index position of reference
    :param variant_chars: all modified or canonical nucleotides
    :param rna: boolean option (default True)
    :param kmer_length: kmer length (default 5)
    """
    min_pos = min(pos)
    max_pos = max(pos)

    sequence = reference_handler.get_sequence(chromosome_name, max(0, min_pos - (kmer_length - 1)),
                                              max_pos + kmer_length)
    if strand == "-":
        sequence = reverse_complement(sequence)
    sequence = list(sequence)
    for pos, chars in zip(pos, variant_chars):
        assert sequence[kmer_length - 1 + (pos - min_pos)] in chars, \
            "Reference base is not in variant characters: pos{}:{} not in {}".format(pos,
                                                                                     sequence[kmer_length - 1 + (
                                                                                                 pos - min_pos)],
                                                                                     chars)
        sequence[kmer_length - 1 + (pos - min_pos)] = chars
    if rna:
        sequence = sequence[::-1]
    kmers = get_kmer_permutations(sequence)
    all_kmers = []
    for x in kmers:
        all_kmers.append(get_kmers_from_seq(x, kmer_length))
    return [set(x) for x in zip(*all_kmers)]


def get_kmer_permutations(bases: list):
    """Given a list of strings, generate a list of permutations of each string for each position

    ex: input: ["AT", "GC"] output-> ["AG", "AC", "TG", "GC"]
    :param bases: list of bases for each position
    """
    switch = [1]
    multiply = 1
    for x in bases[::-1]:
        switch.insert(0, multiply)
        multiply *= len(x)
    num_kmers = multiply
    i = 0
    kmers = []
    kmer = None
    while i < num_kmers:
        kmer = ""
        num = i
        for x, chars in enumerate(bases):
            # print(num, switch[x], num // switch[x], end=" ")
            # print(chars[num // switch[x]])
            kmer += chars[num // switch[x]]
            num -= (num // switch[x]) * switch[x]
        kmers.append(kmer)
        i += 1
    assert len(set(kmers)) == num_kmers, "Made wrong number of kmers: {} != {}".format(len(set(kmer)), num_kmers)
    return kmers


def get_kmers_from_seq(sequence, kmer_length=5):
    """Get list of kmers given a sequence and kmer length
    :param sequence: string
    :param kmer_length: kmer length (default 5)
    """
    assert kmer_length <= len(sequence), \
        "kmer length cannot be larger than sequence length: {} > {}".format(kmer_length, len(sequence))
    kmers = []
    x = 0
    while x <= len(sequence) - kmer_length:
        kmers.append(sequence[x:kmer_length + x])
        x += 1
    return kmers


def get_covered_bases(reference, positions, kmer_length=5, rna=True):
    """Get all covered kmers by position.
    :param reference: path to reference sequence
    :param positions: path to positions file
    :param kmer_length: kmer length
    :param rna: boolean option for rna
    :return list : [[contig, list of positions, list of variants, list of sets of kmers]]
    """
    ref_handler = ReferenceHandler(reference)
    positions_data = CustomAmbiguityPositions.parseAmbiguityFile(positions)
    contig = None
    pos = None
    strand = None
    # base = None
    variant_bases = None

    all_covered_bases = []
    all_pos = []
    all_variant_bases = []
    for i, row in positions_data.iterrows():
        if i > 0:
            all_pos.append(pos)
            all_variant_bases.append(variant_bases)
            if pos + kmer_length <= row[1] or strand != row[2] or contig != row[0]:
                all_covered_bases.append([contig, all_pos, all_variant_bases,
                                          get_covered_kmers(ref_handler, contig, strand, all_pos, all_variant_bases,
                                                            rna, kmer_length)])
                all_pos = []
                all_variant_bases = []
        contig = row[0]
        pos = row[1]
        strand = row[2]
        # base = row[3]
        variant_bases = row[4]
    all_pos.append(pos)
    all_variant_bases.append(variant_bases)
    all_covered_bases.append([contig, all_pos, all_variant_bases,
                              get_covered_kmers(ref_handler, contig, strand, all_pos, all_variant_bases,
                                                rna, kmer_length)])

    return all_covered_bases


def main():
    args = parse_args()
    # load model files
    assert os.path.exists(args.config), "Config file does not exist: {}".format(args.config)
    config = create_dot_dict(load_json(args.config))
    print("Output dir: ", config.save_fig_dir)
    print("Reference: ", config.reference)
    print("Assignments: ", config.assignments)
    print("Positions: ", config.positions)
    print("Models: ", config.models)
    print("RNA: ", config.rna)
    print("Save: ", config.save)
    print("Scatter: ", config.scatter)
    assert os.path.isdir(config.save_fig_dir), "save_fig_dir does not exist: {}".format(config.save_fig_dir)
    assert len(config.assignments) == len(config.models), "Number of models has to equal number of assignment files"
    kmer_data = []
    for path in config.assignments:
        if path:
            kmer_data.append(MultipleModelHandler.read_in_alignment_data(path))
        else:
            kmer_data.append(None)

    if config.rna:
        kmer_length = 5
    else:
        kmer_length = 6

    all_covered_bases = get_covered_bases(config.reference, config.positions, kmer_length, config.rna)
    models = [HmmModel(x, rna=config.rna, name="model" + str(i)) for i, x in enumerate(config.models)]
    mmh = MultipleModelHandler(models, ["t"]*len(models), assignment_data=kmer_data)
    # #
    for contig, pos, vars1, kmers in all_covered_bases:
        dir_name = os.path.join(config.save_fig_dir, "_".join([contig + '_' + str(x) + y for x, y in zip(pos, vars1)]))
        if os.path.exists(dir_name):
            shutil.rmtree(dir_name)
        os.mkdir(dir_name)
        print(contig, pos, vars1, kmers)
        for kmer_group in kmers:
            output_file = None
            if config.save:
                output_file = os.path.join(dir_name, "_".join(kmer_group)+".gif")
            mmh.animate_kmer_distribution(sorted(list(kmer_group)), output_file=output_file, scatter=config.scatter)


if __name__ == '__main__':
    _, time = time_it(main)
    print("Running Time: {} seconds".format(time))
