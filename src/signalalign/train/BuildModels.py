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
from signalalign.utils.commonFunctions import get_first_seq, make_motif_file, get_all_sequence_kmers, make_CCWGG_positions_file, \
    find_ccwgg_motifs
from py3helpers.utils import create_dot_dict, merge_dicts
from signalalign.utils.fileHandlers import FolderHandler
from signalalign.train.trainModels import trainHMMTransitions
from signalalign.signalAlignment import multithread_signal_alignment

PATH_TO_SIGNALALIGN = os.path.abspath("../../signalAlign/")
PATH_TO_BINS = PATH_TO_SIGNALALIGN + "/bin/"
ENTRY_LINE = "blank\t0\tblank\tblank\t{strand}\t0\t0.0\t0.0\t0.0\t{kmer}\t0.0\t0.0\t{prob}\t{event}\t0.0\n"


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


def ccwgg_kmers(sequence_kmers, kmer_length):
    def check_and_add(methyl_kmer):
        normal_kmer = string.translate(methyl_kmer, demethylate)
        if normal_kmer in sequence_kmers:
            labeled_kmers.append(methyl_kmer)

    labeled_kmers = []

    methyl_core1 = "CEAGG"
    methyl_core2 = "CETGG"
    demethylate = string.maketrans("E", "C")

    nucleotides = "ACGT"
    fourmers = [''.join(x) for x in product(nucleotides, repeat=4)]
    threemers = [''.join(x) for x in product(nucleotides, repeat=3)]
    twomers = [''.join(x) for x in product(nucleotides, repeat=2)]
    # NNNNCC*WGGNN

    # NNNNCC*
    if kmer_length == 6:
        for fourmer in fourmers:
            labeled_kmer1 = (fourmer + methyl_core1)[:kmer_length]
            labeled_kmer2 = (fourmer + methyl_core2)[:kmer_length]
            check_and_add(labeled_kmer1)
            check_and_add(labeled_kmer2)

    # NNNCC*W and NNNCC*
    for threemer in threemers:
        labeled_kmer1 = (threemer + methyl_core1)[:kmer_length]
        labeled_kmer2 = (threemer + methyl_core2)[:kmer_length]
        check_and_add(labeled_kmer1)
        check_and_add(labeled_kmer2)

    # NNCC*WG and NNCC*W
    for twomer in twomers:
        labeled_kmer1 = (twomer + methyl_core1)[:kmer_length]
        labeled_kmer2 = (twomer + methyl_core2)[:kmer_length]
        check_and_add(labeled_kmer1)
        check_and_add(labeled_kmer2)
        # C*WGGNN
        if kmer_length == 6:
            labeled_kmer1 = (methyl_core1 + twomer)[1:]
            labeled_kmer2 = (methyl_core2 + twomer)[1:]
            check_and_add(labeled_kmer1)
            check_and_add(labeled_kmer2)

    for onemer in nucleotides:
        # CC*WGGN and C*WGGN
        labeled_kmer1 = methyl_core1 + onemer
        labeled_kmer2 = methyl_core2 + onemer
        if kmer_length == 6:
            check_and_add(labeled_kmer1)
            check_and_add(labeled_kmer2)
        if kmer_length == 5:
            check_and_add(labeled_kmer1[1:])
            check_and_add(labeled_kmer2[1:])
        labeled_kmer1 = (onemer + methyl_core1)[:kmer_length]
        labeled_kmer2 = (onemer + methyl_core2)[:kmer_length]
        check_and_add(labeled_kmer1)
        check_and_add(labeled_kmer2)

    if kmer_length == 5:
        check_and_add(methyl_core1)
        check_and_add(methyl_core2)

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


def find_gatc_motifs(sequence):
    """Find 'GATC' motif in a nucleotide sequence.  """
    motif = "GATC"
    motif_length = len(motif)
    for i, _ in enumerate(sequence):
        if sequence[i:i+motif_length] == motif:
            yield i + 1

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
                           positions=positions_file, targetFile=motif_file, sub=label, n=n, jobs=int(jobs/2))
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
                           positions=positions_file, targetFile=motif_file, n=n, jobs=int(jobs/2), tHdp=t_hdp,
                           cHdp=c_hdp, degenerate=degenerate)
        commands.append(command)

    os.chdir(PATH_TO_BINS)
    procs = [Popen(x.split(), stdout=sys.stdout, stderr=sys.stderr) for x in commands]
    status = [p.wait() for p in procs]
    os.chdir(working_path)
    return


def make_master_assignment_table(assignment_directories):
    def parse_assignment_file(file):
        data = pd.read_table(file,
                             usecols=(0, 1, 2, 3),
                             names=["kmer", "strand", "level_mean", "prob"],
                             dtype={"kmer": np.str, "strand": np.str, "level_mean": np.float64, "prob": np.float64},
                             header=None
                             )
        return data

    assignments = []
    for d in assignment_directories:
        assignments += [x for x in glob.glob(d) if os.stat(x).st_size != 0]
    assignment_dfs = []
    for f in assignments:
        assignment_dfs.append(parse_assignment_file(f))
    return pd.concat(assignment_dfs)


def train_model_transitions(fasta, pcr_fofn, genomic_fofn, degenerate, jobs, positions_file, iterations, batch_size,
                            outpath, t_model, c_model, hdp_type, stateMachine="threeState", t_hdp=None, c_hdp=None,
                            em_iteration=""):
    working_path = os.path.abspath(outpath) + "/"
    model_directory = working_path + "{sm}_{it}".format(sm=stateMachine, it=em_iteration)
    assert os.path.exists(t_model), "Didn't find template model, looked {}".format(t_model)
    assert os.path.exists(c_model), "Didn't find complement model, looked {}".format(c_model)

    methyl_char = get_methyl_char(degenerate)

    os.chdir(PATH_TO_BINS)
    c = "trainModels -fofn={pcr} -fofn={genomic} -X=C -X={methylChar} -r={fasta} -i={iter} -a={batch} --transitions " \
        "-smt={smt} -T={tModel} -C={cModel} -j={jobs} -p={positions} -o={out} " \
        "".format(pcr=pcr_fofn, genomic=genomic_fofn, fasta=fasta, iter=iterations, batch=batch_size,
                  smt=stateMachine, tModel=t_model, cModel=c_model, jobs=jobs, methylChar=methyl_char,
                  positions=positions_file, out=model_directory)
    if t_hdp is not None and c_hdp is not None:
        c += "-tH={tHdp} -cH={cHdp} ".format(tHdp=os.path.abspath(t_hdp), cHdp=os.path.abspath(c_hdp))
    c = PATH_TO_BINS + c
    os.system(c)
    models = [model_directory + "tempFiles_expectations/template_trained.hmm",
              model_directory + "tempFiles_expectations/complement_trained.hmm"]
    if t_hdp is not None and c_hdp is not None:
        models.append(model_directory + "tempFiles_expectations/template.{}.nhdp".format(hdp_type))
        models.append(model_directory + "tempFiles_expectations/complement.{}.nhdp".format(hdp_type))
    os.chdir(working_path)
    return models


def write_kmers(assignments, threshold, max_assignments, kmer_list, entry_line, fH, strands=["t", "c"]):
    for strand in strands:
        by_strand = assignments.ix[(assignments['strand'] == strand) & (assignments['prob'] >= threshold)]
        for k in kmer_list:
            kmer_assignments = by_strand.ix[by_strand['kmer'] == k]
            if kmer_assignments.empty:
                print("missing kmer {}, continuing".format(k))
                continue
            kmer_assignments = kmer_assignments.sort_values(['prob'], ascending=0)
            n = 0
            for _, r in kmer_assignments.iterrows():
                fH.write(
                    entry_line.format(strand=r['strand'], kmer=r['kmer'], event=r['level_mean'], prob=r['prob']))
                n += 1
                if n >= max_assignments:
                    break
            if n < max_assignments:
                print("WARNING didn't find {max} requested assignments for {kmer} only found {found}"
                      "".format(max=max_assignments, kmer=k, found=n))


def make_build_alignment(assignments, degenerate, kmer_length, ref_fasta, n_canonical_assignments,
                         n_methyl_assignments, outfile, threshold):
    seq = get_first_seq(ref_fasta)
    sequence_kmers = list(get_all_sequence_kmers(seq, kmer_length).keys())
    if degenerate == "adenosine":
        methyl_kmers = list(gatc_kmers(sequence_kmers=sequence_kmers, kmerlength=kmer_length))
        methyl_kmers += list(ctag_kmers(sequence_kmers=sequence_kmers, kmerlength=kmer_length))
    else:
        methyl_kmers = list(ccwgg_kmers(sequence_kmers=sequence_kmers, kmer_length=kmer_length))
        methyl_kmers += list(ggwcc_kmers(sequence_kmers=sequence_kmers, kmer_length=kmer_length))
    fH = open(outfile, "w")
    write_kmers(assignments, threshold, n_canonical_assignments, sequence_kmers, ENTRY_LINE, fH)
    write_kmers(assignments, threshold, n_methyl_assignments, methyl_kmers, ENTRY_LINE, fH)
    fH.close()
    return outfile


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


def build_hdp(build_alignment_path, template_model, complement_model, outpath, hdp_type, samples=15000,
              em_iteration=""):
    working_path = os.path.abspath(outpath) + "/"
    build_alignment = os.path.abspath(build_alignment_path)
    t_model = os.path.abspath(template_model)
    c_model = os.path.abspath(complement_model)
    outpath = os.path.abspath(outpath) + "/"
    hdp_pipeline_dir = outpath + "hdpPipeline{}/".format(em_iteration)
    os.makedirs(hdp_pipeline_dir)
    os.chdir(PATH_TO_BINS)
    c = "hdp_pipeline --build_alignment={build} -tM={tModel} -cM={cModel} -Ba=1 -Bb=1 -Ma=1 -Mb=1 -La=1 -Lb=1 " \
        "-s={samples} --verbose --grid_start=50 --grid_end=140 --grid_length=1800 --verbose -o={out} " \
        "--hdp_type=ecoli".format(build=build_alignment, tModel=t_model, cModel=c_model, samples=samples,
                                  out=hdp_pipeline_dir)
    c = PATH_TO_BINS + c
    os.system(c)
    os.chdir(working_path)

    return [hdp_pipeline_dir + "template.{}.nhdp".format(hdp_type),
            hdp_pipeline_dir + "complement.{}.nhdp".format(hdp_type)]


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
        run_parser.add_argument('--config',type=str, help='Path to the (filled in) config file, generated with "generate".')
        return parser.parse_args()

    args = parse_args()

    if not os.path.exists(args.config):
        print("{config} not found run generate-config".format(config=args.config))
        exit(1)
    # Parse config
    args = check_config(args.config)
    # create working directory
    temp_folder = FolderHandler()
    working_path = temp_folder.open_folder(os.path.join(args.output_dir, "tempFiles_trainModels"))

    train_transitions_config = {
        "output_dir": working_path,
        "samples": args.samples,
        "bwa_reference": args.bwa_reference,
        "in_T_Hmm": args.in_T_Hmm,
        "in_C_Hmm": args.in_C_Hmm,
        "templateHdp": args.templateHdp,
        "complementHdp": args.complementHdp,
        "twoD": args.twoD,
        "alignment_file": args.alignment_file,
        "stateMachineType": args.stateMachineType,
        "training_bases": args.train_transitions_options.training_bases,
        "job_count": args.train_transitions_options.job_count,
        "iterations": args.train_transitions_options.iterations,
        "diagonal_expansion": args.train_transitions_options.diagonal_expansion,
        "constraint_trim": args.train_transitions_options.constraint_trim,
        "test": args.test,
        "debug": args.debug}

    # train transitions if we don't know whats up with a new chemistry
    # this really should be training BOTH transitions and emissions using Normal and inv-gaussian distributions
    # TODO create MLE predictors for the emission distributions
    models = trainHMMTransitions(train_transitions_config)
    # alignment args are the parameters to the HMM/HDP model, and don't change

    # make the positions and motif file
    # if args.positions_file is not None and args.motif_file is not None:
    #     assert len(args.positions_file) == 2 and len(args.motif_file) == 2, "need to give training and testing " \
    #                                                                         "positions/motif files"
    #     for i in range(2):
    #         assert os.path.exists(args.positions_file[i]), "Didn't find positions file, looked " \
    #                                                        "{}".format(args.positions_file)
    #         assert os.path.exists(args.motif_file[i]), "Didn't find motif file, looked {}".format(args.motif_file)
    #     positions_file = args.positions_file[0]
    #     motif_file = args.motif_file[0]
    #     test_positions = args.positions_file[1]
    #     test_motifs = args.motif_file[1]
    # else:
    #     # make the positions file
    #     positions_file = make_positions_file(fasta=args.reference,
    #                                          degenerate=args.degenerate,
    #                                          outfile=working_path + "/{}_positions.positions".format(args.degenerate))
    #
    #     # make the motif file
    #     motif_file = make_gatc_or_ccwgg_motif_file(fasta=args.reference,
    #                                                degenerate=args.degenerate,
    #                                                outfile=working_path + "/{}_target.target".format(args.degenerate))
    #     test_positions = positions_file
    #     test_motifs = motif_file

    # make the fofns for training and testing
    # pcr_fofns, gen_fofns = train_test_split_fofn(pcr_reads_dir=args.pcr_reads,
    #                                              genomic_reads_dir=args.genomic_reads,
    #                                              working_directory=working_path,
    #                                              split=args.split)
    # if args.hdp_type == "multiset":
    #     HDP_type = "multisetPriorEcoli"
    # else:
    #     HDP_type = "singleLevelPriorEcoli"

    template_hmm = models[0]
    complement_hmm = models[1]
    template_hdp = models[2]
    complement_hdp = models[3]


    # Create initial arguments for building assignments with a trained model
    alignment_args = {
        "backward_reference": sample.bw_fasta_path,
        "forward_reference": sample.fw_fasta_path,
        "alignment_file": alignment_file,
        "bwa_reference": bwa_reference,
        "destination": working_folder_path,
        "stateMachineType": state_machine_type,
        "in_templateHmm": template_hmm,
        "in_complementHmm": complement_hmm,
        "in_templateHdp": template_hdp,
        "in_complementHdp": complement_hdp,
        "threshold": 0.01,
        "diagonal_expansion": diagonal_expansion,
        "constraint_trim": constraint_trim,
        "target_regions": None,
        "degenerate": None,
        "twoD_chemistry": twoD,
        'track_memory_usage': False,
        'get_expectations': False,
        'output_format': "assignments"
    }

    multithread_signal_alignment(alignment_args, list_of_fast5s, workers)
    alignments = [x for x in glob.glob(os.path.join(working_folder_path.path, "*.tsv")) if os.stat(x).st_size != 0]

    working_directories = [d + "tempFiles_alignment/*.assignments" for d in working_directories]
    assignment_dirs = run_guide_alignment(fasta=reference_location,
                                          pcr_fofn=pcr_fofns[0],
                                          genomic_fofn=gen_fofns[0],
                                          jobs=args.jobs,
                                          positions_file=positions_file,
                                          motif_file=motif_file,
                                          n=args.n_aligns,
                                          degenerate=args.degenerate,
                                          t_model=models[0],
                                          c_model=models[1],
                                          outpath=working_path)
    assert kmer_length_from_model(models[0]) == kmer_length_from_model(models[1]), "Models had different kmer lengths"
    # concatenate the assignments into table
    master = make_master_assignment_table(assignment_dirs)
    if args.bulk is True:
        build_alignment = make_bulk_build_alignment(assignments=master,
                                                    degenerate=args.degenerate,
                                                    n_canonical_assignments=args.assignments,
                                                    n_methyl_assignments=args.methyl_assignments,
                                                    threshold=args.assignment_threshold,
                                                    outfile=working_path + "/buildAlignment.tsv")
    else:
        build_alignment = make_build_alignment(assignments=master,
                                               degenerate=args.degenerate,
                                               kmer_length=kmer_length_from_model(models[0]),
                                               ref_fasta=reference_location,
                                               n_canonical_assignments=args.assignments,
                                               n_methyl_assignments=args.methyl_assignments,
                                               outfile=working_path + "/buildAlignment.tsv",
                                               threshold=args.assignment_threshold)
    # build hdp
    hdps = build_hdp(build_alignment_path=build_alignment,
                     template_model=models[0],
                     complement_model=models[1],
                     hdp_type=HDP_type,
                     outpath=working_path,
                     samples=args.samples)

    if args.HDP_EM is not None:
        hdp_models = HDP_EM(ref_fasta=reference_location,
                            pcr_fofn=pcr_fofns[0],
                            gen_fofn=gen_fofns[0],
                            degenerate=args.degenerate,
                            jobs=args.jobs,
                            positions_file=positions_file,
                            motif_file=motif_file,
                            n_assignment_alns=args.n_aligns,
                            n_canonical_assns=args.assignments,
                            n_methyl_assns=args.methyl_assignments,
                            iterations=args.iterations,
                            batch_size=args.batch,
                            working_path=working_path,
                            start_hdps=hdps,
                            threshold=args.assignment_threshold,
                            start_temp_hmm=models[0],
                            start_comp_hmm=models[1],
                            n_iterations=args.HDP_EM,
                            gibbs_samples=args.samples,
                            bulk=args.bulk,
                            hdp_type=HDP_type)
    else:
        # train HMM/HDP
        hdp_models = train_model_transitions(fasta=reference_location,
                                             pcr_fofn=pcr_fofns[0],
                                             genomic_fofn=gen_fofns[0],
                                             degenerate=args.degenerate,
                                             jobs=args.jobs,
                                             positions_file=positions_file,
                                             iterations=args.iterations,
                                             batch_size=args.batch,
                                             outpath=working_path,
                                             stateMachine="threeStateHdp",
                                             t_hdp=hdps[0],
                                             c_hdp=hdps[1],
                                             hdp_type=HDP_type,
                                             t_model=os.path.abspath(args.in_T_Hmm),
                                             c_model=os.path.abspath(args.in_C_Hmm))
    # run methylation variant calling experiment
    run_variant_calling_experiment(fasta=reference_location,
                                   pcr_fofn=pcr_fofns[1],
                                   genomic_fofn=gen_fofns[1],
                                   jobs=args.jobs,
                                   positions_file=test_positions,
                                   motif_file=test_motifs,
                                   t_model=hdp_models[0],
                                   c_model=hdp_models[1],
                                   outpath=working_path,
                                   n=args.n_test_alns,
                                   degenerate=args.degenerate,
                                   t_hdp=hdp_models[2],
                                   c_hdp=hdp_models[3])
    # run the control experiment
    #run_variant_calling_experiment(fasta=os.path.abspath(args.reference),
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
