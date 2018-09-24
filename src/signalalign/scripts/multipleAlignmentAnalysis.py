#!/usr/bin/env python3
"""Run signal-to-reference alignments
"""
from __future__ import print_function

import numpy as np
import glob
import shutil
import subprocess
import os
import sys
from argparse import ArgumentParser
import scipy
import math
from signalalign.utils.sequenceTools import reverse_complement

# signal align indices
CONTIG_IDX=0
REF_POS_IDX=1
REF_KMER_IDX=2
PATH_KMER_IDX=9
POSTERIOR_PROB_IDX=12
KMER_IDX=15

# data structure keys
REF_CHAR_KEY='r'
ALIGNED_FILE_COUNT_KEY='afc'
ALIGNED_KMER_COUNT_KEY='akc'
POSTERIORS_KEY='p'
CHARACTERS_KEY='c'


def parse_args(args=None):
    parser = ArgumentParser(description=__doc__)

    parser.add_argument('--alignment_glob', '-i', action='store',
                        dest='alignment_glob', required=True, type=str, default=None,
                        help="glob matching signalAlign alignments")
    parser.add_argument('--output_file', '-o', action='store',
                        dest='output_file', required=False, type=str, default=None,
                        help="output filename; defaults to stdout")

    args = parser.parse_args(args)

    return args


def log(msg):
    print(msg, file=sys.stderr)


def is_forward_alignment(filename):
    with open(filename) as file:
        for line in file:
            line = line.split("\t")
            ref_kmer = line[REF_KMER_IDX].upper()
            path_kmer = line[PATH_KMER_IDX].upper()

            if ref_kmer != path_kmer:
                # must be aligned backwards
                return False
            elif reverse_complement(ref_kmer) != path_kmer:
                # they are the same, and we know they are not aligned backwards and palindromic
                return True
            else:
                # we can't be sure, we must keep looking
                continue
    log("Cannot infer direction of alignment: {}".format(filename))
    return None


def get_empty_pos_entry(ref_char=None):
    return {
        REF_CHAR_KEY: ref_char,
        ALIGNED_FILE_COUNT_KEY: 0,
        ALIGNED_KMER_COUNT_KEY: 0,
        POSTERIORS_KEY: list(),
        CHARACTERS_KEY: dict(),
    }


# def shannon_entropy(probabilities):
#     total_prob = sum(probabilities)
#     probabilities = list(map(lambda x: x / total_prob, probabilities))
#     entropy = sum(map(lambda x: x * math.log2(x), probabilities))
#     if entropy == 0.0: return entropy
#     entropy /= -2.0
#     return entropy


# def shannon_entropy(probabilities):
#     entropy = sum(map(lambda x: x * math.log2(x), probabilities))
#     if entropy == 0.0: return entropy
#     entropy /= -2.0
#     entropy /= len(probabilities)
#     return entropy


def shannon_entropy(probabilities):
    if len(probabilities) == 0: return 0.0
    entropy = np.mean(list(map(lambda x: 0.0 if x == 0 or x == 1 else x * math.log2(x) + (1.0-x) * math.log2(1.0-x), probabilities)))
    if entropy == 0.0: return entropy
    entropy /= -2.0
    return entropy


def qualify_entropies(entropies):
    max_entropy = max(entropies)
    buckets = {i:0 for i in range(16+1)}
    for e in entropies:
        b = int(16.0 * e / max_entropy)
        buckets[b] += 1
    max_bucket = max(buckets.values())
    for i in range(16):
        log("%2.6f: %s" % (i / max_entropy, "#" * int(32.0 * buckets[i] / max_bucket)))


def main(args):

    # parse args
    args = parse_args(args)
    files = glob.glob(args.alignment_glob)
    file_count = len(files)
    if file_count == 0:
        log("No files matching {}".format(args.alignment_glob))
        return 1

    """
    { 
        chrom : {
            ref_pos : {
                ref_char : ''
                aligned_file_count : 0
                aligned_kmer_count : 0
                posteriors : []
                characters : {A:[], C:[], .. }
            }
        }
    }
    """
    alignment_data = dict()
    all_characters = set()
    for fcnt, file in enumerate(files):
        forward = is_forward_alignment(file)
        log("%5d (%2d%%): %s %s" % (fcnt, int(100.0 * fcnt/file_count), 'f' if forward else 'b', file))

        min_pos = sys.maxsize
        max_pos = 0
        chrom = None
        with open(file) as input:
            for line in input:
                line = line.split("\t")
                chrom = line[CONTIG_IDX]
                ref_pos = int(line[REF_POS_IDX])
                kmer = line[KMER_IDX].strip().upper()
                path_kmer = line[PATH_KMER_IDX].strip().upper()
                prob = float(line[POSTERIOR_PROB_IDX])

                # store data
                if chrom not in alignment_data: alignment_data[chrom] = dict()
                for i, c in enumerate(kmer):
                    # get characters and index
                    if forward:
                        char_pos = ref_pos + i
                        char = c
                        ref_char = path_kmer[i]
                    else:
                        char_pos = ref_pos - i + len(kmer) - 1
                        char = reverse_complement(c)
                        ref_char = reverse_complement(path_kmer[i])
                    all_characters.add(char)

                    # positional stuff
                    min_pos = min(min_pos, char_pos)
                    max_pos = max(max_pos, char_pos)

                    # ensure data structure
                    if char_pos not in alignment_data[chrom]:
                        alignment_data[chrom][char_pos] = get_empty_pos_entry(ref_char)
                    if alignment_data[chrom][char_pos][REF_CHAR_KEY] is None:
                        alignment_data[chrom][char_pos][REF_CHAR_KEY] = ref_char
                    if char not in alignment_data[chrom][char_pos][CHARACTERS_KEY]:
                        alignment_data[chrom][char_pos][CHARACTERS_KEY][char] = list()
                    # sanity check
                    if alignment_data[chrom][char_pos][REF_CHAR_KEY] != ref_char:
                        log("PROGRAMMER ERROR regarding ref position, likely due to fwd/bkwd orientation")
                        return 1

                    # store data
                    alignment_data[chrom][char_pos][POSTERIORS_KEY].append(prob)
                    alignment_data[chrom][char_pos][ALIGNED_KMER_COUNT_KEY] += 1
                    alignment_data[chrom][char_pos][CHARACTERS_KEY][char].append(prob)

            if chrom is None:
                log("\tNo alignments found!")
            else:
                max_pos += 1
                for i in range(min_pos, max_pos):
                    if i not in alignment_data[chrom]:
                        alignment_data[chrom][i] = get_empty_pos_entry()
                    alignment_data[chrom][i][ALIGNED_FILE_COUNT_KEY] += 1

    # data
    output = None
    try:
        all_characters = list(all_characters)
        all_characters.sort()
        entropies = list()

        # header
        output = open(args.output_file, 'w') if args.output_file is not None else sys.stdout
        output.write("#contig\tposition\tref_char\tshannon_entropy\tkmer_aln_count\tread_aln_count")
        for char in all_characters:
            output.write("\tp{}".format(char))
        output.write("\n")

        # save all the data
        chroms = list(alignment_data.keys())
        chroms.sort()
        for chrom in chroms:
            positions = list(alignment_data[chrom].keys())
            positions.sort()
            for pos in positions:
                ref_char = alignment_data[chrom][pos][REF_CHAR_KEY]
                if ref_char is None: ref_char = ""
                entropy = shannon_entropy(alignment_data[chrom][pos][POSTERIORS_KEY])
                kmer_aln_count = alignment_data[chrom][pos][ALIGNED_KMER_COUNT_KEY]
                read_aln_count = alignment_data[chrom][pos][ALIGNED_FILE_COUNT_KEY]

                output.write("{contig}\t{position}\t{ref_char}\t{shannon_entropy}\t{kmer_aln_count}\t{read_aln_count}"
                             .format(contig=chrom, position=pos, ref_char=ref_char, shannon_entropy=entropy,
                                     kmer_aln_count=kmer_aln_count, read_aln_count=read_aln_count))

                chars_at_pos = alignment_data[chrom][pos][CHARACTERS_KEY]
                total_char_prob = sum(map(sum, chars_at_pos.values())) \
                    if alignment_data[chrom][pos][REF_CHAR_KEY] is not None else 1.0
                for char in all_characters:
                    output.write("\t{}".format(0.0 if char not in chars_at_pos else
                                               sum(chars_at_pos[char]) / total_char_prob))

                output.write("\n")
                entropies.append(entropy)

        qualify_entropies(entropies)

    finally:
        if args.output_file is not None and output is not None: output.close()




if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
