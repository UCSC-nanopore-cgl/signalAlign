"""Functions that are useful for exploring data
"""
import os
import random
import pandas as pd
import numpy as np
from collections import Counter
from signalalign.utils.parsers import read_fasta
from signalalign.utils.sequenceTools import reverse_complement
from signalalign.utils import kmer_iterator


def get_all_sequence_kmers(seq, k=5):
    kmers = Counter()
    for kmer in kmer_iterator(seq, k):
        kmers[kmer] += 1
        kmers[reverse_complement(kmer, reverse=True, complement=True)] += 1
    return kmers


def parse_alignment_file(alignment_file):
    data = pd.read_table(alignment_file, usecols=(1, 4, 5, 9, 12, 13),
                         dtype={'ref_pos': np.int64,
                                'strand': np.str,
                                'event_index': np.int64,
                                'kmer': np.str,
                                'posterior_prob': np.float64,
                                'event_mean': np.float64},
                         header=None,
                         names=['ref_pos', 'strand', 'event_index', 'kmer', 'posterior_prob', 'event_mean'])
    return data


def get_first_seq(fasta_file):
    assert os.path.isfile(fasta_file)
    i = 0
    seq = ""
    for t, c, s in read_fasta(fasta_file):
        seq += s
        i += 1
        if i > 0:
            break
    return seq


def find_occurences(seq, ch):
    return [i for i, letter in enumerate(seq) if letter == ch]


def parse_motif_file(motif_file):
    starts = []
    for line in open(motif_file, 'r'):
        line = line.split()
        start = int(line[0]) + 5
        starts.append(start)
    return starts


def check_starts(starts, seq, target):
    for start in starts:
        assert seq[start] == target
    return True


def find_ccwgg_motifs(seq):
    motifs = ["CCAGG", "CCTGG"]
    motif_length = len(motifs[0])
    for i, _ in enumerate(seq):
        if seq[i:i + motif_length] in motifs:
            yield i + 1  # + 1 because we want to label the second C in CCWGG


def make_CCWGG_positions_file(fasta, out_file):
    def check_starts(starts, seq, target):
        for start in starts:
            assert seq[start] == target
        return True

    if os.path.exists(out_file):
        print("Outfile you're trying to make already exists")
        return
    seq = get_first_seq(fasta)
    starts = [x for x in find_ccwgg_motifs(seq)]
    fH = open(out_file, 'w')
    if check_starts(starts, seq, "C"):
        fH.write("X\t")
        for start in starts:
            fH.write("{}\t".format(start))
    fH.write("\n")
    # + 2 because the starts are for the second C in CCWGG and on the complement we want to mark the the C in the motif
    rc_starts = [x + 2 for x in starts]
    if check_starts(rc_starts, seq, "G"):
        fH.write("X\t")
        for start in rc_starts:
            fH.write("{}\t".format(start))
    fH.write("\n")
    fH.close()
    return out_file


# this function just marks all of the base pairs in starts
def make_substitute_file(starts, seq, out_file):
    if os.path.exists(out_file):
        print("Outfile you're trying to make already exists")
        return
    fH = open(out_file, 'w')

    fH.write("X\t")
    for start in starts:
        fH.write("{}\t".format(start))
    fH.write("\n")

    fH.write("X\t")
    for start in starts:
        fH.write("{}\t".format(start))
    fH.write("\n")

    fH.close()

# makes a motif file that is passed to SignalAlign with the -q flag
def make_motif_file(starts, seq, outfile):
    if os.path.exists(outfile):
        print("Outfile you're trying to make already exists")
        return
    fH = open(outfile, 'w')
    for start in starts:
        s = start - 5
        e = start + 6
        q = seq[s:e]
        fH.write("{start}\t{end}\t{motif}\n".format(start=s, end=e, motif=q))


def random_positions(seq, n, start, end):
    """
    seq: sequence, needed for length
    n: number of random ints to get
    """
    N = []
    l = len(seq)
    for i in range(n):
        N.append(random.randint(start, end))
    print("Made {} random ints".format(len(N)))
    return list(set(N))


def group_sites_in_window2(sites, window=6):
    def collect_group(start):
        i = start
        g = [sites[start]]
        while sites[i + 1] - sites[i] < window:
            g.append(sites[i + 1])
            i += 1
            if len(sites) <= i + 1:
                break
        return g, i + 1

    sites.sort()
    groups = []
    i = 0
    while i + 1 < len(sites):
        g, i = collect_group(i)
        groups.append(g)
    return groups


def parse_substitution_file(substitution_file):
    fH = open(substitution_file, 'r')
    line = fH.readline().split()
    forward_sub = line[0]
    forward_pos = list(map(np.int64, line[1:]))
    line = fH.readline().split()
    backward_sub = line[0]
    backward_pos = list(map(np.int64, line[1:]))
    return (forward_sub, forward_pos), (backward_sub, backward_pos)
