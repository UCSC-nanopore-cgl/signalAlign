#!/usr/bin/env python

import os
import sys
import glob
import pandas as pd
import numpy as np
import string
import yaml
from argparse import ArgumentParser
from subprocess import Popen
from itertools import product
from random import shuffle
from shutil import copyfile
from py3helpers.utils import create_dot_dict, merge_dicts, merge_lists, all_string_permutations
from signalalign.utils.fileHandlers import FolderHandler
from signalalign.utils.parsers import read_fasta
from signalalign.train.trainModels import trainHMM
from signalalign.signalAlignment import multithread_signal_alignment_samples, create_signalAlignment_args, \
    SignalAlignSample
from signalalign.utils.sequenceTools import find_gatc_motifs, count_all_sequence_kmers, get_motif_kmers, \
    get_sequence_kmers
from signalalign.utils.commonFunctions import get_first_seq, make_motif_file, \
    make_CCWGG_positions_file, \
    find_ccwgg_motifs

PATH_TO_SIGNALALIGN = os.path.abspath("../../signalAlign/")
PATH_TO_BINS = PATH_TO_SIGNALALIGN + "/bin/"


def make_alignment_line(strand, kmer, prob, event):
    """Convert strand, kmer, probability and event to a correctly formatted alignment file line
    :param strand: 't' or 'c' representing template or complement strand of read
    :param kmer: nucleotide kmer
    :param prob: probability of kmer coming from certain event
    :param event: mean of the corresponding event
    :return: correctly formatted alignment line
    """
    assert strand in ['c', 't'], "Strand must be either 'c' or 't'. strand: {}".format(strand)
    assert isinstance(prob, float), "probability must be a float: prob {}".format(prob)
    assert isinstance(kmer, str), "kmer must be a string: kmer {}".format(kmer)
    assert isinstance(event, float), "event must be a float: event {}".format(event)
    entry_line = "blank\t0\tblank\tblank\t{strand}\t0\t0.0\t0.0\t0.0\t{kmer}\t0.0\t0.0\t{prob}\t{event}\t0.0\n"
    return entry_line.format(strand=strand, kmer=kmer, prob=prob, event=event)


def train_test_split_fofn(pcr_reads_dir, genomic_reads_dir, working_directory, split=0.5):
    def split_files(files):
        shuffle(files)
        split_point = int(split * len(files))
        return files[:split_point], files[split_point:]

    def write_fofn(outfile, files):
        with open(outfile, "w") as fH:
            for f in files:
                fH.write("{filepath}\n".format(filepath=f))
        return outfile

    pcr_read_files = [pcr_reads_dir + x for x in os.listdir(pcr_reads_dir) if x.endswith(".fast5")]
    genomic_read_files = [genomic_reads_dir + x for x in os.listdir(genomic_reads_dir) if x.endswith(".fast5")]
    assert len(pcr_read_files) > 0 and len(genomic_read_files) > 0
    train_pcr_files, test_pcr_files = split_files(pcr_read_files)
    train_gen_files, test_gen_files = split_files(genomic_read_files)

    pcr_train_fofn = write_fofn(working_directory + "/pcr_train.fofn", train_pcr_files)
    pcr_test_fofn = write_fofn(working_directory + "/pcr_test.fofn", test_pcr_files)

    gen_train_fofn = write_fofn(working_directory + "/gen_train.fofn", train_gen_files)
    gen_test_fofn = write_fofn(working_directory + "/gen_test.fofn", test_gen_files)

    return (pcr_train_fofn, pcr_test_fofn), (gen_train_fofn, gen_test_fofn)


def kmer_length_from_model(model_file):
    with open(model_file, "r") as model:
        line = model.readline().split()
        kmer_length = int(line[-1])
        model.close()
    assert kmer_length == 5 or kmer_length == 6
    return kmer_length


def get_methyl_char(degenerate):
    if degenerate == "adenosine":
        return "I"
    else:
        return "E"


def make_positions_file(fasta, degenerate, outfile):
    if degenerate == "adenosine":
        return make_gatc_position_file(fasta, outfile)
    else:
        return make_CCWGG_positions_file(fasta, outfile)


def gatc_kmers(sequence_kmers, kmerlength):
    assert kmerlength == 5 or kmerlength == 6, "only works with kmer lengths 5 and 6"
    # NNNNGATCNNN
    methyl_core = "GITC"
    normal_core = "GATC"
    nucleotides = "ACGT"

    fourmers = [''.join(x) for x in product(nucleotides, repeat=4)]
    threemers = [''.join(x) for x in product(nucleotides, repeat=3)]
    twomers = [''.join(x) for x in product(nucleotides, repeat=2)]

    labeled_kmers = []

    # add NNNNGA*
    if kmerlength == 6:
        for fourmer in fourmers:
            labeled_kmer = (fourmer + methyl_core)[:kmerlength]
            normal_kmer = (fourmer + normal_core)[:kmerlength]
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
    # add NNNGA*T and NNNGA*
    for threemer in threemers:
        labeled_kmer = (threemer + methyl_core)[:kmerlength]
        normal_kmer = (threemer + normal_core)[:kmerlength]
        if normal_kmer in sequence_kmers:
            labeled_kmers.append(labeled_kmer)
        # A*TCNNN
        if kmerlength == 6:
            labeled_kmer = (methyl_core + threemer)[1:]
            normal_kmer = (normal_core + threemer)[1:]
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
    # add NNGA*TC and NNGA*T
    for twomer in twomers:
        labeled_kmer = (twomer + methyl_core)[:kmerlength]
        normal_kmer = (twomer + normal_core)[:kmerlength]
        if normal_kmer in sequence_kmers:
            labeled_kmers.append(labeled_kmer)
        # A*TCNN
        if kmerlength == 5:
            labeled_kmer = (methyl_core + twomer)[1:]
            normal_kmer = (normal_core + twomer)[1:]
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
        # NGA*TCN
        if kmerlength == 6:
            labeled_kmer = (twomer[0] + methyl_core + twomer[1])
            normal_kmer = (twomer[0] + normal_core + twomer[1])
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
            # TODO forgot GA*TCNN for 6mers!
            labeled_kmer = (methyl_core + twomer)
            normal_kmer = (normal_core + twomer)
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
    if kmerlength == 5:
        for onemer in "ACTG":
            labeled_kmer = onemer + methyl_core
            normal_kmer = onemer + normal_core
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
            labeled_kmer = methyl_core + onemer
            normal_kmer = normal_core + onemer
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
    return set(labeled_kmers)


def ctag_kmers(sequence_kmers, kmerlength):
    assert kmerlength == 5 or kmerlength == 6, "only works with kmer lengths 5 and 6"
    # NNNCTAGNNNN
    methyl_core = "CTIG"
    normal_core = "CTAG"
    nucleotides = "ACGT"

    fourmers = [''.join(x) for x in product(nucleotides, repeat=4)]
    threemers = [''.join(x) for x in product(nucleotides, repeat=3)]
    twomers = [''.join(x) for x in product(nucleotides, repeat=2)]

    labeled_kmers = []

    # add A*GNNNN
    if kmerlength == 6:
        for fourmer in fourmers:
            labeled_kmer = (methyl_core + fourmer)[2:]
            normal_kmer = (normal_core + fourmer)[2:]
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
    # add NNNCTA*
    for threemer in threemers:
        if kmerlength == 6:
            labeled_kmer = (threemer + methyl_core)[:kmerlength]
            normal_kmer = (threemer + normal_core)[:kmerlength]
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
            labeled_kmer = (methyl_core + threemer)[1:]
            normal_kmer = (normal_core + threemer)[1:]
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
        # A*GNNN
        if kmerlength == 5:
            labeled_kmer = (methyl_core + threemer)[2:]
            normal_kmer = (normal_core + threemer)[2:]
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)

    # add NNCTA*G and NNCTA*
    for twomer in twomers:
        labeled_kmer = (twomer + methyl_core)[:kmerlength]
        normal_kmer = (twomer + normal_core)[:kmerlength]
        if normal_kmer in sequence_kmers:
            labeled_kmers.append(labeled_kmer)
        # CTA*GNN
        if kmerlength == 6:
            labeled_kmer = (methyl_core + twomer)[:kmerlength]
            normal_kmer = (normal_core + twomer)[:kmerlength]
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
        # TA*GNN
        if kmerlength == 5:
            labeled_kmer = (methyl_core + twomer)[1:]
            normal_kmer = (normal_core + twomer)[1:]
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)

    if kmerlength == 5:
        for onemer in nucleotides:
            labeled_kmer = onemer + methyl_core
            normal_kmer = onemer + normal_core
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
            labeled_kmer = methyl_core + onemer
            normal_kmer = normal_core + onemer
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
    return set(labeled_kmers)


def ggwcc_kmers(sequence_kmers, kmer_length):
    def check_and_add(methyl_kmer):
        normal_kmer = string.translate(methyl_kmer, demethylate)
        if normal_kmer in sequence_kmers:
            labeled_kmers.append(methyl_kmer)

    labeled_kmers = []

    methyl_core1 = "GGAEC"
    methyl_core2 = "GGTEC"
    demethylate = string.maketrans("E", "C")

    nucleotides = "ACGT"
    fourmers = [''.join(x) for x in product(nucleotides, repeat=4)]
    threemers = [''.join(x) for x in product(nucleotides, repeat=3)]
    twomers = [''.join(x) for x in product(nucleotides, repeat=2)]

    # NNGGWC*CNNN

    # C*CNNNN
    for fourmer in fourmers:
        labeled_kmer1 = (methyl_core1 + fourmer)[3:]
        labeled_kmer2 = (methyl_core2 + fourmer)[3:]
        check_and_add(labeled_kmer1)
        check_and_add(labeled_kmer2)

    # WC*CNNN and C*CNNN
    for threemer in threemers:
        labeled_kmer1 = (methyl_core1 + threemer)[2:] if kmer_length == 6 else (methyl_core1 + threemer)[3:]
        labeled_kmer2 = (methyl_core2 + threemer)[2:] if kmer_length == 6 else (methyl_core2 + threemer)[3:]
        check_and_add(labeled_kmer1)
        check_and_add(labeled_kmer2)

    # GWC*CNN and WC*CNN
    for twomer in twomers:
        labeled_kmer1 = (methyl_core1 + twomer)[1:] if kmer_length == 6 else (methyl_core1 + twomer)[2:]
        labeled_kmer2 = (methyl_core2 + twomer)[1:] if kmer_length == 6 else (methyl_core2 + twomer)[2:]
        check_and_add(labeled_kmer1)
        check_and_add(labeled_kmer2)
        # NNGGWC*
        if kmer_length == 6:
            labeled_kmer1 = (twomer + methyl_core1)[:kmer_length]
            labeled_kmer2 = (twomer + methyl_core2)[:kmer_length]
            check_and_add(labeled_kmer1)
            check_and_add(labeled_kmer2)

    for onemer in nucleotides:
        # NGGWC* and NGGWC*C
        labeled_kmer1 = (onemer + methyl_core1)[:kmer_length]
        labeled_kmer2 = (onemer + methyl_core2)[:kmer_length]
        check_and_add(labeled_kmer1)
        check_and_add(labeled_kmer2)
        # GGWC*CN GWC*CN
        labeled_kmer1 = methyl_core1 + onemer if kmer_length == 6 else (methyl_core1 + onemer)[1:]
        labeled_kmer2 = methyl_core2 + onemer if kmer_length == 6 else (methyl_core2 + onemer)[1:]
        check_and_add(labeled_kmer1)
        check_and_add(labeled_kmer2)

    if kmer_length == 5:
        check_and_add(methyl_core1)
        check_and_add(methyl_core2)

    return set(labeled_kmers)


def motif_kmers(core, kmer_length=5, multiplier=5):
    motifs = []
    repeat = kmer_length - len(core)
    if repeat == 0:
        return [core] * multiplier
    else:
        for k in product("ACGT", repeat=repeat):
            fix = ''.join(k)
            motifs.append(fix + core)
            motifs.append(core + fix)
        return motifs * multiplier


def nuc_position(seq_str, char):
    """Finds all positions of specific character
        withing sequence"""
    motif_position = [m.start() for m in re.finditer(char, seq_str)]
    return motif_position


def make_gatc_position_file(fasta, outfile):
    outfile = os.path.abspath(outfile)
    fH = open(outfile, 'w')
    fH.write("X\t")

    seq = get_first_seq(fasta)
    for i in find_gatc_motifs(seq):
        assert seq[i] == "A"
        fH.write("{}\t".format(i))
    fH.write("\n")
    fH.write("X\t")
    for i in find_gatc_motifs(seq):
        t_pos = i + 1
        assert seq[t_pos] == "T"
        fH.write("{}\t".format(t_pos))  # + 1 because this is for the reverse complement
    fH.write("\n")
    fH.close()
    return outfile


def make_gatc_or_ccwgg_motif_file(fasta, degenerate, outfile):
    if degenerate == "adenosine":
        motif_finder = find_gatc_motifs
    else:
        motif_finder = find_ccwgg_motifs
    outfile = os.path.abspath(outfile)
    seq = get_first_seq(fasta)
    positions = [x for x in motif_finder(seq)]
    make_motif_file(positions, seq, outfile)
    return outfile


def run_guide_alignment(fasta, pcr_fofn, genomic_fofn, jobs, positions_file, motif_file, t_model, c_model,
                        outpath, n, degenerate, em_iteration="", t_hdp=None, c_hdp=None):
    def get_labels():
        if degenerate == "adenosine":
            return ["A", "I"]
        else:
            return ["C", "E"]

    working_path = os.path.abspath(outpath)
    commands = []
    c = PATH_TO_BINS + "runSignalAlign -fofn={fofn} -r={fasta} -T={tModel} -C={cModel} -f=assignments " \
                       "-o={outpath} -p={positions} -q={targetFile} -X={sub} -n={n} -j={jobs} "

    if t_hdp is not None and c_hdp is not None:
        c += "-tH={tHdp} -cH={cHdp} ".format(tHdp=os.path.abspath(t_hdp), cHdp=os.path.abspath(c_hdp))

    read_sets = [pcr_fofn, genomic_fofn]
    labels = get_labels()
    working_directories = [working_path + "/pcr_{}".format(em_iteration),
                           working_path + "/genomic_{}".format(em_iteration)]

    assert os.path.exists(t_model), "Didn't find template model, looked {}".format(t_model)
    assert os.path.exists(c_model), "Didn't find complement model, looked {}".format(c_model)

    for fofn, label, working_directory in zip(read_sets, labels, working_directories):
        # assemble the command
        command = c.format(fofn=fofn, fasta=fasta, tModel=t_model, cModel=c_model, outpath=working_directory,
                           positions=positions_file, targetFile=motif_file, sub=label, n=n, jobs=int(jobs / 2))
        commands.append(command)

    os.chdir(PATH_TO_BINS)
    procs = [Popen(x.split(), stdout=sys.stdout, stderr=sys.stderr) for x in commands]
    status = [p.wait() for p in procs]

    os.chdir(working_path)
    working_directories = [d + "tempFiles_alignment/*.assignments" for d in working_directories]

    return working_directories


def run_variant_calling_experiment(fasta, pcr_fofn, genomic_fofn, jobs, positions_file, motif_file, t_model, c_model,
                                   outpath, n, degenerate, t_hdp, c_hdp):
    working_path = os.path.abspath(outpath)
    commands = []
    c = PATH_TO_BINS + "runSignalAlign -fofn={fofn} -r={fasta} -T={tModel} -C={cModel} -f=variantCaller " \
                       "-o={outpath} -p={positions} -q={targetFile} -n={n} -j={jobs} -x={degenerate} " \
                       "-tH={tHdp} -cH={cHdp} -smt=threeStateHdp"
    read_sets = [pcr_fofn, genomic_fofn]
    working_directories = [working_path + "/{}_pcr_".format(degenerate),
                           working_path + "/{}_genomic_".format(degenerate)]

    assert os.path.exists(t_model), "Didn't find template model, looked {}".format(t_model)
    assert os.path.exists(c_model), "Didn't find complement model, looked {}".format(c_model)
    assert os.path.exists(t_hdp), "Didn't find template model, looked {}".format(t_hdp)
    assert os.path.exists(c_hdp), "Didn't find complement model, looked {}".format(c_hdp)

    for fofn, working_directory in zip(read_sets, working_directories):
        # assemble the command
        command = c.format(fofn=fofn, fasta=fasta, tModel=t_model, cModel=c_model, outpath=working_directory,
                           positions=positions_file, targetFile=motif_file, n=n, jobs=int(jobs / 2), tHdp=t_hdp,
                           cHdp=c_hdp, degenerate=degenerate)
        commands.append(command)

    os.chdir(PATH_TO_BINS)
    procs = [Popen(x.split(), stdout=sys.stdout, stderr=sys.stderr) for x in commands]
    status = [p.wait() for p in procs]
    os.chdir(working_path)
    return


def parse_assignment_file(file_path):
    """Parse the .assignments.tsv output file from signalAlign:

    :param file_path: path to assignments file
    :return: panda DataFrame with column names "kmer", "strand", "level_mean", "prob"
    """
    data = pd.read_table(file_path,
                         usecols=(0, 1, 2, 3),
                         names=["kmer", "strand", "level_mean", "prob"],
                         dtype={"kmer": np.str, "strand": np.str, "level_mean": np.float64, "prob": np.float64},
                         header=None
                         )
    return data


def make_master_assignment_table(list_of_assignment_paths):
    """Create a master assignment table from a list of assignment paths

    :param list_of_assignment_paths: list of all paths to assignment.tsv files to concat
    :return: pandas DataFrame of all assignments
    """
    assignment_dfs = []
    for f in list_of_assignment_paths:
        assignment_dfs.append(parse_assignment_file(f))
    return pd.concat(assignment_dfs)


class CreateHdpTrainingData(object):
    """Process the assignment files created from SignalAlign for the HDP distribution estimation"""

    def __init__(self, samples, out_file_path, k, template=True, complement=False, verbose=True):
        """
        Control how each kmer/event assignment is processed given a set of samples and the parameters associated with
        each sample

        :param samples:
        :param out_file:
        :param template: generate kmers for template read strand: default: True
        :param complement: generate kmers for complement read strand: default: True
        :param min_probability: the minimum probability to use for assigning kmers
        :param verbose: option to print update statements
        """
        self.strands = []
        if template:
            self.strands.append('t')
        if complement:
            self.strands.append('c')
        assert self.strands != [], 'template or complement need to be set to True. ' \
                                   'complement: {}, template: {}'.format(complement, template)

        for sample in samples:
            assert isinstance(sample, SignalAlignSample)

        self.samples = samples
        self.out_file_path = out_file_path
        self.template = template
        self.complement = complement
        self.verbose = verbose
        self.k = k
        self.master_assignment_table = \
            make_master_assignment_table(merge_lists([sample.analysis_files for sample in self.samples]))

    def generate_hdp_training_lines_wrapper(self, kmer_list, max_assignments=100, min_probability=0.8):
        """Convert assignments to alignment line format for HDP training.

        Filter assignments on a minimum probability, read strand, and a max number of kmer assignments

        :param kmer_list: list of kmers to write to alignment file
        :param max_assignments: max number of assignments to process for each kmer
        :param min_probability: the minimum probability to use for assigning kmers
        """
        # loop through for each strand in the assignments
        return self._generate_hdp_training_lines(self.master_assignment_table, kmer_list,
                                                 max_assignments=max_assignments,
                                                 strands=self.strands, min_probability=min_probability,
                                                 verbose=self.verbose)

    def write_hdp_training_file(self):
        """Write a hdp training file to a specified location"""
        with open(self.out_file_path, 'w') as out_file:
            for sample in self.samples:
                kmers = self.get_sample_kmers(sample)
                for line in self.generate_hdp_training_lines_wrapper(kmers,
                                                                     max_assignments=sample.number_of_kmer_assignments,
                                                                     min_probability=sample.probability_threshold):
                    out_file.write(line)
        return self.out_file_path

    def get_sample_kmers(self, sample):
        """Get all kmers from a sample, either from the reference sequence, all possible kmers from an alphabet or
            all kmers that cover a modified nucelotide

        :param sample: AbstractSamples object
        :return: set of desired kmers
        """
        kmers = set()
        # if motifs is present, process for all motifs with modified base
        if sample.motifs:
            for motif in sample.motifs:
                kmers |= get_motif_kmers(motif, self.k, alphabet="ATGC")
        # if we want to limit kmers which were seen in reference sequence
        elif sample.kmers_from_reference:
            for _, _, sequence in read_fasta(sample.bwa_reference):
                kmers |= get_sequence_kmers(sequence, k=self.k, rev_comp=True)
        else:
            kmers |= {x for x in all_string_permutations("ATGC", length=self.k)}

        return kmers

    @staticmethod
    def _generate_hdp_training_lines(assignments, kmer_list, max_assignments=10,
                                     strands=['t', 'c'], min_probability=0.8, verbose=False):
        """Convert assignments to alignment line format for HDP training.

        Filter assignments on a minimum probability, read strand, and a max number of kmer assignments

        :param assignments: pandas array of assignments to search
        :param template: generate kmers for template read strand: default: True
        :param complement: generate kmers for complement read strand: default: True
        :param verbose: option to print update statements
        :param kmer_list: list of kmers to write to alignment file
        :param max_assignments: max number of assignments to process for each kmer
        :param min_probability: the minimum probability to use for assigning kmers
        """
        # loop through for each strand in the assignments
        assert isinstance(strands, list) and len(strands) > 0, "strands must be a list and not be empty. strands: {}" \
                                                               "".format(strands)
        for strand in strands:
            by_strand = assignments.loc[(assignments['strand'] == strand)
                                        & (assignments['prob'] >= min_probability)]

            for k in kmer_list:
                kmer_assignments = by_strand.loc[by_strand['kmer'] == k]
                if kmer_assignments.empty and verbose:
                    print("missing kmer {}, continuing".format(k), file=sys.stderr)
                    continue
                kmer_assignments = kmer_assignments.sort_values(['prob'], ascending=0)
                n = 0
                for _, r in kmer_assignments.iterrows():
                    yield make_alignment_line(strand=r['strand'], kmer=r['kmer'], event=r['level_mean'], prob=r['prob'])
                    n += 1
                    if n >= max_assignments:
                        break
                if n < max_assignments and verbose:
                    print("WARNING didn't find {max} requested assignments for {kmer} only found {found}"
                          "".format(max=max_assignments, kmer=k, found=n))


def make_bulk_build_alignment(assignments, degenerate, n_canonical_assignments, n_methyl_assignments, threshold,
                              outfile, strands=["t", "c"]):
    def write_bulk_assignments(assignment_df, max_assignments):
        for strand in strands:
            by_strand = assignment_df.ix[(assignment_df['strand'] == strand) & (assignment_df['prob'] >= threshold)]
            n = 0
            for i, r in by_strand.iterrows():
                fH.write(
                    entry_line.format(strand=r['strand'], kmer=r['kmer'], event=r['level_mean'], prob=r['prob']))
                n += 1
                if n >= max_assignments:
                    break
            if n < max_assignments:
                print("WARNING: only found {found} assignments when {requested} were requested"
                      "".format(found=n, requested=max_assignments))

    fH = open(outfile, "w")
    methyl_char = get_methyl_char(degenerate)
    entry_line = "blank\t0\tblank\tblank\t{strand}\t0\t0.0\t0.0\t0.0\t{kmer}\t0.0\t0.0\t{prob}\t{event}\t0.0\n"
    labeled = assignments.loc[assignments["kmer"].str.contains(methyl_char)]
    canonical = assignments.loc[~(assignments["kmer"].str.contains(methyl_char))]
    write_bulk_assignments(labeled, n_methyl_assignments)
    write_bulk_assignments(canonical, n_canonical_assignments)
    fH.close()
    return outfile


def build_hdp(build_alignment_path, template_model, complement_model, outpath, hdp_type, path_to_bin, samples=15000,
              em_iteration=""):
    """Python wrapper to the hdp_pipeline script to train the HDP models for the kmer distributions"""
    assert os.path.exists(build_alignment_path), "Build alignment path does not exist. build_alignment_path: " \
                                                 "{}".format(build_alignment_path)
    assert os.path.exists(template_model), "template_model does not exist. complement_model: {}".format(template_model)
    assert os.path.exists(complement_model), "complement_model does not exist. complement_model: {}".format(
        complement_model)
    assert os.path.exists(outpath), "outpath does not exist. outpath: {}".format(outpath)

    hdp_pipeline_dir = os.path.join(outpath, "hdpPipeline{}/".format(em_iteration))
    os.makedirs(hdp_pipeline_dir)
    os.chdir(path_to_bin)
    hdp_pipeline = os.path.join(path_to_bin, "hdp_pipeline")
    c = "{hdp_pipeline} --build_alignment={build} -tM={tModel} -cM={cModel} -Ba=1 -Bb=1 -Ma=1 -Mb=1 -La=1 -Lb=1 " \
        "-s={samples} --verbose --grid_start=50 --grid_end=140 --grid_length=1800 --verbose -o={out} " \
        "--hdp_type={hdp_type}".format(hdp_pipeline=hdp_pipeline, build=build_alignment_path, tModel=template_model,
                                       cModel=complement_model, samples=samples,
                                       out=hdp_pipeline_dir, hdp_type=hdp_type)
    os.system(c)
    os.chdir(outpath)

    template_hdp_model = os.path.join(hdp_pipeline_dir, "template.{}.nhdp".format(hdp_type))
    complement_hdp_model = os.path.join(hdp_pipeline_dir, "complement.{}.nhdp".format(hdp_type))
    assert os.path.exists(complement_hdp_model), "Complement hdp model does not exist".format(complement_hdp_model)
    assert os.path.exists(template_hdp_model), "Template hdp model does not exist. {}".format(template_hdp_model)

    return [template_hdp_model, complement_hdp_model]


def count_lines_txt_file(txt_file):
    """Count number of lines in a file
    :param txt_file: path to text file of any type
    :return: number of lines in file
    """
    count = 0
    with open(txt_file, 'r') as txt_fh:
        for _ in txt_fh:
            count += 1
    return count


def kmer_length_from_2_models(template_model_file, complement_model_file):
    def get_kmer_length(model_file):
        with open(model_file, "r") as fH:
            line = fH.readline().split()
            assert len(line) == 4, "HdpPipeline ERROR: wrong header in model file {}".format(model_file)
            kmer_length = int(line[3])
            fH.close()
        return kmer_length
    template_kmer_length = get_kmer_length(template_model_file)
    complement_kmer_length = get_kmer_length(complement_model_file)
    assert template_kmer_length == complement_kmer_length
    return template_kmer_length


def get_hdp_type(requested_type):
    """Get the integer representing the model type for the buildHdpUtils.c program

    :param requested_type: string associated with each type of hdp model
    :return: integer associated with that
    """
    hdp_types = {
        "singleLevelFixed": 0,
        "singleLevelPrior": 1,
        "multisetFixed": 2,
        "multisetPrior": 3,
        "compFixed": 4,
        "compPrior": 5,
        "middleNtsFixed": 6,
        "middleNtsPrior": 7,
        "groupMultisetFixed": 8,
        "groupMultisetPrior": 9,
        "singleLevelPrior2": 10,
        "multisetPrior2": 11,
        "multisetPriorEcoli": 12,
        "singleLevelPriorEcoli": 13
    }
    assert (requested_type in list(hdp_types.keys())), "Requested HDP type is invalid, got {}".format(requested_type)
    return hdp_types[requested_type]


def get_initial_hdp_args(args, hdp_type):
    """For each type of HDP we need to set certain parameters.

    :param args: arguments from create_hdp_training_args
    :param hdp_type: integer representation of hdp model
    :return: commands for buildHdpUtil.c with correct hdp parameters
    """
    # if we're making a HDP with fixed concentration parameters
    if hdp_type in [0, 2, 4, 6, 8]:
        assert None not in [args.base_gamma, args.leaf_gamma], \
            "ERROR: need to specify concentration parameters for type {}".format(hdp_type)
        if hdp_type == 0: # singleLevelFixed
            return "-B {base} -L {leaf} ".format(base=args.base_gamma, leaf=args.leaf_gamma)
        else:  # multisetFixed/ compFixed / middleNtsFixed / groupMultisetFixed
            assert args.middle_gamma is not None, "ERROR: need to specify middle concentration param"
            return "-B {base} -M {middle} -L {leaf} ".format(base=args.base_gamma, middle=args.middle_gamma,
                                                             leaf=args.leaf_gamma)
    else: # 'Prior' type models
        assert None not in [args.base_alpha, args.base_beta, args.leaf_alpha, args.leaf_beta], \
            "ERROR: missing Gamma prior hyper parameters"
        if hdp_type == 1 or hdp_type == 10:  # singleLevelPrior, singleLevelPrior2
            return "-g {Ba} -r {Bb} -i {La} -u {Lb} ".format(Ba=args.base_alpha, Bb=args.base_beta,
                                                             La=args.leaf_alpha, Lb=args.leaf_beta)
        else:  # multisetPrior, compPrior, middleNtsPrior, groupMultisetPrior, multisetPrior2,
            #    singleLevelPriorEcoli, multisetPriorEcoli
            assert None not in [args.middle_alpha, args.middle_beta], "ERROR: need middle hyper parameters"
            return "-g {Ba} -r {Bb} -j {Ma} -y {Mb} -i {La} -u {Lb} ".format(Ba=args.base_alpha, Bb=args.base_beta,
                                                                             Ma=args.middle_alpha, Mb=args.middle_beta,
                                                                             La=args.leaf_alpha, Lb=args.leaf_beta)


# globals
HDP_TYPES = [
    ("singleLevelFixed", 0),
    ("singleLevelPrior", 1),
    ("multisetFixed", 2),
    ("multisetPrior", 3),
    ("compFixed", 4),
    ("compPrior", 5),
    ("middleNtsFixed", 6),
    ("middleNtsPrior", 7),
    ("groupMultisetFixed", 8),
    ("groupMultisetPrior", 9),
]

HDP_TYPES_2 = [
    ("singleLevelPrior2", 10),
    ("multisetPrior2", 11),
]

HDP_TYPES_ECOLI = [
    ("multisetPriorEcoli", 12),
    ("singleLevelPriorEcoli", 13),
]


def create_hdp_training_args(C_alignments=None, mC_alignments=None, hmC_alignments=None,
                             number_of_assignments=10000, build_alignment=None, threshold=0.0,
                             hdp_type="ecoli", template_model=None, complement_model=None,
                             base_gamma=1.0, middle_gamma=1.0, leaf_gamma=1.0, base_alpha=1.0,
                             base_beta=1.0, middle_alpha=1.0, middle_beta=1.0, leaf_alpha=1.0,
                             leaf_beta=1.0, samples=10000, thinning=100, verbose=True, grid_start=30.0,
                             grid_end=90.0, grid_length=1200, outpath=None, path_to_bin='', num_alignments=None,
                             twoD=True):
    """
    :param outpath: output file path
    :param C_alignments:
    :param mC_alignments:
    :param hmC_alignments:
    :param number_of_assignments: total number of assignments to collect FOR EACH GROUP
    :param build_alignment: path to alignment file
    :param num_alignments: number of alignments in alignment file
    :param threshold:
    :param verbose:
    :param path_to_bin
    :param twoD:
    :param hdp_type: Build Hdp, specify type, options: "Prior, Fixed, twoWay. twoWay is a Prior-type model (recommended)"
    # initial HDP
    :param template_model: Input template lookup table
    :param complement_model: Input complement lookup table
    # fixed concentration models
    :param base_gamma:
    :param middle_gamma:
    :param leaf_gamma:
    # gamma prior models
    :param base_alpha:
    :param base_beta:
    :param middle_alpha:
    :param middle_beta:
    :param leaf_alpha:
    :param leaf_beta:
    # gibbs sampling
    :param samples: number of gibbs samples
    :param thinning: how many thinning draws?
    # sample grid
    :param grid_start:
    :param grid_end:
    :param grid_length:
    :return: dictionary of hdp training options
    """

    args = dict(C_alignments=C_alignments, mC_alignments=mC_alignments, hmC_alignments=hmC_alignments,
                number_of_assignments=number_of_assignments, build_alignment=build_alignment, threshold=threshold,
                hdp_type=hdp_type, template_model=template_model, complement_model=complement_model,
                base_gamma=base_gamma, middle_gamma=middle_gamma, leaf_gamma=leaf_gamma, base_alpha=base_alpha,
                base_beta=base_beta, middle_alpha=middle_alpha, middle_beta=middle_beta, leaf_alpha=leaf_alpha,
                leaf_beta=leaf_beta, samples=samples, thinning=thinning, verbose=verbose, grid_start=grid_start,
                grid_end=grid_end, grid_length=grid_length, outpath=outpath, path_to_bin=path_to_bin,
                num_alignments=num_alignments, twoD=twoD)
    return create_dot_dict(args)


def train_hdp(args):
    """Format arguments for the buildHdpUtil.c program which trains a HDP to model kmer distributions.

    arguments are outlined in 'create_hdp_training_args'
    """
    # Pipeline Script
    working_directory = os.path.abspath(args.outpath)  # this is the directory we will use for everything
    assert os.path.isdir(working_directory), "ERROR: the working directory you specified doesn't exist." \
                                             "args.out: {}".format(args.out)
    # move alignments to working directory
    build_alignment_location = os.path.join(working_directory, "buildAlignment.tsv")
    assert os.path.isfile(args.build_alignment), "ERROR: Didn't find input BuildAlignment"
    copyfile(args.build_alignment, build_alignment_location)
    # count alignments
    if args.num_alignments:
        approx_total_build_assignments = args.num_alignments
    else:
        approx_total_build_assignments = count_lines_txt_file(build_alignment_location)
    # initial HDP
    assert (os.path.isfile(build_alignment_location)), "ERROR: Didn't find build alignment"
    buildHdpUtil = os.path.join(args.path_to_bin, "./buildHdpUtil")
    assert (os.path.exists(buildHdpUtil)), "ERROR: Didn't find buildHdpUtil. {}".format(buildHdpUtil)
    print("[hdp_pipeline] NOTICE: Making initial HDP of type {}\n".format(args.hdp_type), file=sys.stderr)

    kmer_length = kmer_length_from_2_models(args.template_model, args.complement_model)
    template_lookup_table = " -T" + args.template_model
    complement_lookup_table = " -C" + args.complement_model
    verbose_flag = "--verbose " if args.verbose is True else ""



    if args.hdp_type == "cytosine2":
        hdp_types = HDP_TYPES_2
    elif args.hdp_type == "ecoli":
        hdp_types = HDP_TYPES_ECOLI
    elif args.hdp_type == "Prior":
        hdp_types = HDP_TYPES[1::2]
    else:
        hdp_types = [(args.hdp_type, get_hdp_type(args.hdp_type))]

    if args.twoD:
        oneD = None
    else:
        oneD = '--oneD'
        assert set(hdp_types) <= set(HDP_TYPES_2)


    build_commands = []
    for hdp_type, i, in hdp_types:
        template_hdp_location = os.path.join(working_directory, "template." + hdp_type + ".nhdp")
        complement_hdp_location = os.path.join(working_directory, "complement." + hdp_type + ".nhdp")

        build_initial_hdp_command = "{buildHdpUtil} {verbose}-p {hdpType} -v {tHdpLoc} -w {cHdpLoc} -l {buildAln} " \
                                    "-a {kmerLength} -n {samples} -I {burnIn} -t {thin} -s {start} -e {end} " \
                                    "-k {len}{tL}{cL} {oneD} " \
                                    "".format(buildHdpUtil=buildHdpUtil, hdpType=i, tHdpLoc=template_hdp_location,
                                              cHdpLoc=complement_hdp_location, buildAln=build_alignment_location,
                                              samples=args.samples, burnIn=32 * approx_total_build_assignments,
                                              thin=args.thinning, start=args.grid_start, end=args.grid_end,
                                              len=args.grid_length, verbose=verbose_flag, tL=template_lookup_table,
                                              cL=complement_lookup_table, kmerLength=kmer_length, oneD=oneD)
        build_initial_hdp_command += get_initial_hdp_args(args=args, hdp_type=i)
        build_commands.append(build_initial_hdp_command)
        print("[hdp_pipeline] Command: {}\n".format(build_initial_hdp_command))

    # initial_hdp_build_out = open(os.path.join(working_directory + "build_initial_hdp.out"), 'w')
    # initial_hdp_build_err = open(os.path.join(working_directory + "build_initial_hdp.err"), 'w')
    procs = [Popen(x.split(), stdout=sys.stdout, stderr=sys.stderr) for x in build_commands]

    # procs = [Popen(x.split(), stdout=initial_hdp_build_out, stderr=initial_hdp_build_err) for x in build_commands]
    status = [p.wait() for p in procs]

    # initial_hdp_build_out.close()
    # initial_hdp_build_err.close()

    print("[pipeline] DONE.", file=sys.stderr)
    template_hdp_models = [os.path.join(working_directory, "template.{}.nhdp".format(hdp_type)) for hdp_type, i, in hdp_types]
    for template_hdp_model in template_hdp_models:
        assert os.path.exists(template_hdp_model), "Template hdp model does not exist. {}".format(template_hdp_model)
    complement_hdp_models = None

    if args.twoD:
        complement_hdp_models = [os.path.join(working_directory, "complement.{}.nhdp".format(hdp_type)) for hdp_type, i, in hdp_types]
        for complement_hdp_model in complement_hdp_models:
            assert os.path.exists(complement_hdp_model), "Complement hdp model does not exist".format(complement_hdp_model)

    return [template_hdp_models, complement_hdp_models]


def HDP_EM(ref_fasta, pcr_fofn, gen_fofn, degenerate, jobs, positions_file, motif_file, n_assignment_alns,
           n_canonical_assns, n_methyl_assns, iterations, batch_size, working_path, start_hdps, threshold,
           hdp_type, start_temp_hmm, start_comp_hmm, n_iterations, gibbs_samples, bulk):
    template_hdp = start_hdps[0]
    complement_hdp = start_hdps[1]
    template_hmm = start_temp_hmm
    complement_hmm = start_comp_hmm
    for i in range(n_iterations):
        # first train the model transitions
        hdp_models = train_model_transitions(fasta=ref_fasta,
                                             pcr_fofn=pcr_fofn,
                                             genomic_fofn=gen_fofn,
                                             degenerate=degenerate,
                                             jobs=jobs,
                                             positions_file=positions_file,
                                             iterations=iterations,
                                             batch_size=batch_size,
                                             outpath=working_path,
                                             em_iteration="{}_".format(i),
                                             stateMachine="threeStateHdp",
                                             t_hdp=template_hdp,
                                             c_hdp=complement_hdp,
                                             t_model=template_hmm,
                                             c_model=complement_hmm,
                                             hdp_type=hdp_type)
        # next get assignments
        assignment_dirs = run_guide_alignment(fasta=ref_fasta,
                                              pcr_fofn=pcr_fofn,
                                              genomic_fofn=gen_fofn,
                                              jobs=jobs,
                                              positions_file=positions_file,
                                              motif_file=motif_file,
                                              n=n_assignment_alns,
                                              em_iteration="{}_".format(i),
                                              degenerate=degenerate,
                                              t_model=hdp_models[0],
                                              c_model=hdp_models[1],
                                              t_hdp=hdp_models[2],
                                              c_hdp=hdp_models[3],
                                              outpath=working_path)

        assert kmer_length_from_model(hdp_models[0]) == kmer_length_from_model(
            hdp_models[1]), "Models had different kmer lengths"

        # assemble them into a big table
        master = make_master_assignment_table(assignment_directories=assignment_dirs)
        # make the build alignment of assignments
        if bulk is True:
            build_alignment = make_bulk_build_alignment(assignments=master,
                                                        degenerate=degenerate,
                                                        n_canonical_assignments=n_canonical_assns,
                                                        n_methyl_assignments=n_methyl_assns,
                                                        threshold=threshold,
                                                        outfile=working_path + "/buildAlignment_{}.tsv".format(i))
        else:
            build_alignment = make_build_alignment(assignments=master,
                                                   degenerate=degenerate,
                                                   kmer_length=kmer_length_from_model(hdp_models[0]),
                                                   ref_fasta=ref_fasta,
                                                   n_canonical_assignments=n_canonical_assns,
                                                   n_methyl_assignments=n_methyl_assns,
                                                   outfile=working_path + "/buildAlignment_{}.tsv".format(i),
                                                   threshold=threshold)
        # make new hdps
        new_hdps = build_hdp(build_alignment_path=build_alignment,
                             template_model=hdp_models[0],
                             complement_model=hdp_models[1],
                             hdp_type=hdp_type,
                             outpath=working_path,
                             samples=gibbs_samples,
                             em_iteration="_{}".format(i))
        template_hdp, template_hmm = new_hdps[0], hdp_models[0]
        complement_hdp, complement_hmm = new_hdps[1], hdp_models[1]

    return template_hmm, complement_hmm, template_hdp, complement_hdp


def check_config(config):
    """Make sure training configuration file is correctly filled out

    :param config: path to BuildModels yaml configuration file
    """
    params = create_dot_dict({x.replace('-', '_'): y for x, y in yaml.load(open(config).read()).items()})
    # assert os.path.isfile(args.reference)
    assert os.path.exists(params.output_dir)
    assert os.path.isfile(params.in_T_Hmm)
    assert os.path.isfile(params.in_C_Hmm)

    return params


def main(args):
    # def parse_args():
    #     parser = ArgumentParser(description=__doc__)
    #     parser.add_argument("-r", action="store", dest="reference", required=True)
    #     parser.add_argument("-pcr", action="store", dest="pcr_reads", required=True)
    #     parser.add_argument("-gen", action="store", dest="genomic_reads", required=True)
    #     parser.add_argument('--in_template_hmm', '-T', action='store', dest='in_T_Hmm', required=True, type=str)
    #     parser.add_argument('--in_complement_hmm', '-C', action='store', dest='in_C_Hmm', required=True, type=str)
    #     parser.add_argument("-o", action="store", dest="outpath", required=True)
    #     parser.add_argument("-x", action="store", dest="degenerate", required=True)
    #     parser.add_argument('--positions', action='append', dest='positions_file', required=False, default=None)
    #     parser.add_argument('--motif', action='append', dest='motif_file', required=False, default=None)
    #     parser.add_argument('--bulk', action='store_true', dest='bulk', required=False, default=False)
    #     parser.add_argument('--hdp_type', action='store', dest='hdp_type', required=False, default="multiset")
    #     parser.add_argument("-j", action="store", dest="jobs", required=False, default=4, type=int)
    #     parser.add_argument("-i", action="store", dest="iterations", required=False, type=int, default=20)
    #     parser.add_argument("-a", action="store", dest="batch", required=False, type=int, default=15000)
    #     parser.add_argument("-s", action="store", dest="assignments", required=False, type=int, default=30)
    #     parser.add_argument("-c", action="store", dest="methyl_assignments", required=False, type=int, default=200)
    #     parser.add_argument("-n", action="store", dest="n_aligns", required=False, type=int, default=1000)
    #     parser.add_argument("-e", action="store", dest="n_test_alns", required=False, type=int, default=1000)
    #     parser.add_argument("-t", action="store", dest="assignment_threshold", required=False, type=float, default=0.8)
    #     parser.add_argument("-g", action="store", dest="samples", required=False, type=int, default=15000)
    #     parser.add_argument("--hdp_em", action="store", dest="HDP_EM", required=False, type=int, default=None)
    #     parser.add_argument("--split", action="store", dest="split", required=False, type=float, default=0.5)
    #     args = parser.parse_args()
    #     return args

    # command_line = " ".join(sys.argv[:])
    # print("Command Line: {cmdLine}\n".format(cmdLine=command_line), file=sys.stderr)

    def parse_args():
        parser = ArgumentParser()
        subparsers = parser.add_subparsers(dest="command")

        # parsers for running the full pipeline
        run_parser = subparsers.add_parser("run", help="runs full workflow ")
        run_parser.add_argument('--config', type=str,
                                help='Path to the (filled in) config file, generated with "generate".')
        return parser.parse_args()

    args = parse_args()

    if not os.path.exists(args.config):
        print("{config} not found run generate-config".format(config=args.config))
        exit(1)
    # Parse config
    args = check_config(args.config)
    # create working directory
    working_folder = FolderHandler()
    working_path = working_folder.open_folder(os.path.join(args.output_dir, "tempFiles_trainModels"))
    # process reference sequences and samples
    samples = [SignalAlignSample(working_folder=working_folder, **s) for s in args.samples]

    train_transitions_config = {
        "output_dir": working_path,
        "bwa_reference": args.bwa_reference,
        "in_T_Hmm": args.in_T_Hmm,
        "in_C_Hmm": args.in_C_Hmm,
        "templateHdp": args.templateHdp,
        "complementHdp": args.complementHdp,
        "twoD": args.twoD,
        "alignment_file": args.alignment_file,
        "stateMachineType": args.stateMachineType,
        "training_bases": args.train_transitions_options.training_bases,
        "job_count": args.job_count,
        "iterations": args.train_transitions_options.iterations,
        "diagonal_expansion": args.train_transitions_options.diagonal_expansion,
        "constraint_trim": args.train_transitions_options.constraint_trim,
        "test": args.test,
        "debug": args.debug,
        "path_to_bin": args.path_to_bin}

    # train transitions if we don't know whats up with a new chemistry
    # this really should be training BOTH transitions and emissions using Normal and inv-gaussian distributions
    # TODO create MLE predictors for the emission distributions
    models = trainHMM(samples, working_folder, train_transitions_config)
    # alignment args are the parameters to the HMM/HDP model, and don't change

    template_hmm = models[0]
    complement_hmm = models[1]
    template_hdp = models[2]
    complement_hdp = models[3]

    # Create initial arguments for building assignments with a trained model
    alignment_args = create_signalAlignment_args(alignment_file=args.alignment_file,
                                                 bwa_reference=args.bwa_reference,
                                                 destination=working_path,
                                                 stateMachineType=args.stateMachineType,
                                                 in_templateHmm=template_hmm,
                                                 in_complementHmm=complement_hmm,
                                                 in_templateHdp=template_hdp,
                                                 in_complementHdp=complement_hdp,
                                                 twoD_chemistry=args.twoD,
                                                 output_format="assignments",
                                                 path_to_bin=args.path_to_bin)

    samples = multithread_signal_alignment_samples(samples, alignment_args, args.job_count)
    # concatenate the assignments into table
    template = True
    complement = False
    if args.twoD:
        complement = True

    kmer_length = kmer_length_from_model(models[0])

    build_alignment_path = CreateHdpTrainingData(samples, os.path.join(working_path, "buildAlignment.tsv"), kmer_length,
                                                 template=template,
                                                 complement=complement,
                                                 verbose=True).write_hdp_training_file()
    print(build_alignment_path)
    # master = make_master_assignment_table(merge_lists(sample_files.keys()))
    # if args.bulk is True:
    #     build_alignment = make_bulk_build_alignment(assignments=master,
    #                                                 degenerate=args.degenerate,
    #                                                 n_canonical_assignments=args.assignments,
    #                                                 n_methyl_assignments=args.methyl_assignments,
    #                                                 threshold=args.assignment_threshold,
    #                                                 outfile=working_path + "/buildAlignment.tsv")
    # else:
    #     build_alignment = make_build_alignment(assignments=master,
    #                                            degenerate=args.degenerate,
    #                                            kmer_length=kmer_length_from_model(models[0]),
    #                                            ref_fasta=reference_location,
    #                                            n_canonical_assignments=args.assignments,
    #                                            n_methyl_assignments=args.methyl_assignments,
    #                                            outfile=working_path + "/buildAlignment.tsv",
    #                                            threshold=args.assignment_threshold)

    # build hdp
    hdps = build_hdp(build_alignment_path=build_alignment_path,
                     template_model=models[0],
                     complement_model=models[1],
                     hdp_type=args.train_HDP_expectations_options.hdp_type,
                     outpath=working_path,
                     samples=10,
                     path_to_bin=args.path_to_bin)

    # if args.HDP_EM is not None:
    #     hdp_models = HDP_EM(ref_fasta=reference_location,
    #                         pcr_fofn=pcr_fofns[0],
    #                         gen_fofn=gen_fofns[0],
    #                         degenerate=args.degenerate,
    #                         jobs=args.jobs,
    #                         positions_file=positions_file,
    #                         motif_file=motif_file,
    #                         n_assignment_alns=args.n_aligns,
    #                         n_canonical_assns=args.assignments,
    #                         n_methyl_assns=args.methyl_assignments,
    #                         iterations=args.iterations,
    #                         batch_size=args.batch,
    #                         working_path=working_path,
    #                         start_hdps=hdps,
    #                         threshold=args.assignment_threshold,
    #                         start_temp_hmm=models[0],
    #                         start_comp_hmm=models[1],
    #                         n_iterations=args.HDP_EM,
    #                         gibbs_samples=args.samples,
    #                         bulk=args.bulk,
    #                         hdp_type=HDP_type)
    # else:
    #     # train HMM/HDP
    #     hdp_models = train_model_transitions(fasta=reference_location,
    #                                          pcr_fofn=pcr_fofns[0],
    #                                          genomic_fofn=gen_fofns[0],
    #                                          degenerate=args.degenerate,
    #                                          jobs=args.jobs,
    #                                          positions_file=positions_file,
    #                                          iterations=args.iterations,
    #                                          batch_size=args.batch,
    #                                          outpath=working_path,
    #                                          stateMachine="threeStateHdp",
    #                                          t_hdp=hdps[0],
    #                                          c_hdp=hdps[1],
    #                                          hdp_type=HDP_type,
    #                                          t_model=os.path.abspath(args.in_T_Hmm),
    #                                          c_model=os.path.abspath(args.in_C_Hmm))
    # # run methylation variant calling experiment
    # run_variant_calling_experiment(fasta=reference_location,
    #                                pcr_fofn=pcr_fofns[1],
    #                                genomic_fofn=gen_fofns[1],
    #                                jobs=args.jobs,
    #                                positions_file=test_positions,
    #                                motif_file=test_motifs,
    #                                t_model=hdp_models[0],
    #                                c_model=hdp_models[1],
    #                                outpath=working_path,
    #                                n=args.n_test_alns,
    #                                degenerate=args.degenerate,
    #                                t_hdp=hdp_models[2],
    #                                c_hdp=hdp_models[3])
    # run the control experiment
    # run_variant_calling_experiment(fasta=os.path.abspath(args.reference),
    #                               pcr_reads=os.path.abspath(args.pcr_reads) + "/",
    #                               genomic_reads=os.path.abspath(args.genomic_reads) + "/",
    #                               jobs=args.jobs,
    #                               positions_file=positions_file,
    #                               motif_file=motif_file,
    #                               t_model=hdp_models[0],
    #                               c_model=hdp_models[1],
    #                               outpath=working_path,
    #                               n=args.n_test_alns,
    #                               degenerate="variant",
    #                               t_hdp=hdp_models[2],
    #                               c_hdp=hdp_models[3])


if __name__ == "__main__":
    sys.exit(main(sys.argv))
