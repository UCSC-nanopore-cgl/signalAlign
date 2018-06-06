#!/usr/bin/env python
"""This script takes a directory of assignments, which are tab-seperated files containing k-mer to ionic
current pairs with a posterior probability it consolidates them, and reformats the whole thing for 
HdpBuildUtil and runs it"""

import sys
import os
import re
import subprocess
from itertools import product
from argparse import ArgumentParser
from signalalign.train.BuildModels import \
    make_master_assignment_table,\
    make_bulk_build_alignment,\
    kmer_length_from_model,\
    write_kmers,\
    ENTRY_LINE

PATH_TO_SIGNALALIGN = "../../signalAlign/bin"
PATH_TO_5MER_MODEL  = "../../signalAlign/models/testModelR9p4_5mer_acegt_template.model"
NUCLEOTIDES         = "ACGT"
NUCLEOTIDES_METHYLC = "ACGTE"


def getKmers(kmer_length):
    return ["".join(k) for k in product(NUCLEOTIDES, repeat=kmer_length)]


def count_lines_in_build_alignment(build_alignment_path):
    count = 0
    with open(build_alignment_path, "r") as fH:
        for line in fH:
            count += 1
    return count


def getCGMotifKmers(kmer_length=5):
    def find_CGs(k):
        return [m.start() for m in re.finditer("CG", k)]

    def find_GCs(k):
        return [m.start() + 1 for m in re.finditer("GC", k)]

    def methylate_kmers(kmers, positions, methyl="E"):
        for kmer, pos in zip(kmers, positions):
            l_kmer = list(kmer)
            for p in pos:
                l_kmer[p] = methyl
            yield ''.join(l_kmer)

    # handle CG containing kmers
    kmers     = getKmers(kmer_length)
    cg_kmers  = [k for k in kmers if "CG" in k]
    CGs       = [find_CGs(k) for k in cg_kmers]
    mcg_kmers = [k for k in methylate_kmers(cg_kmers, CGs)]
    # handle the complement (GC) kmers
    gc_kmers  = [k for k in kmers if "GC" in k]
    GCs       = [find_GCs(k) for k in gc_kmers]
    mgc_kmers = [k for k in methylate_kmers(gc_kmers, GCs)]
    return mcg_kmers + mgc_kmers


def make_cg_build_alignment(assignments, n_assignments, outfile, threshold=0.8, kmer_length=5):
    kmers = ["".join(k) for k in product(NUCLEOTIDES_METHYLC, repeat=kmer_length)]
    fH = open(outfile, "w")
    strands = ["t"]
    write_kmers(assignments, threshold, n_assignments, kmers, ENTRY_LINE, fH, strands=strands)
    fH.close()


def main(args):
    def parse_args():
        parser = ArgumentParser(description=__doc__)
        parser.add_argument("-pcr", action="store", dest="pcr", required=True)
        parser.add_argument("-gen", action="store", dest="gen", required=True)
        parser.add_argument("-o", action="store", dest="out_dir", required=True)

        parser.add_argument("--bulk", action="store_true", dest="bulk", required=False, default=False)
        parser.add_argument("-s", action="store", dest="assignments", required=False, type=int, default=30)
        parser.add_argument("-c", action="store", dest="methyl_assignments", required=False, type=int, default=200)
        parser.add_argument("-d", action="store", dest="assignment_threshold", required=False, type=float, default=0.8)
        parser.add_argument("-g", action="store", dest="samples", required=False, type=int, default=15000)
        parser.add_argument('-t', action='store', type=int, default=100, dest='thinning')

        return parser.parse_args()

    args = parse_args()
    print("Commandline: %s " % " ".join(sys.argv[:]), file=sys.stderr)

    workdir = os.path.abspath(args.out_dir)

    # routines to make the `build alignment` that is input to the HDP builder
    assignments = make_master_assignment_table([args.pcr, args.gen])
    out_aln_filename    = "build_alignment_CG_{mC}_{C}_{t}.tsv".format(mC=args.methyl_assignments,
                                                                       C=args.assignments,
                                                                       t=args.assignment_threshold)
    out_build_alignment = os.path.join(workdir, out_aln_filename)
    if args.bulk:
        make_bulk_build_alignment(assignments, "cytosine2",
                                  n_canonical_assignments=args.assignments,
                                  n_methyl_assignments=args.methyl_assignments,
                                  threshold=args.assignment_threshold,
                                  outfile=out_build_alignment,
                                  strands=["t"])
        assert(os.path.exists(out_build_alignment))
    else:
        make_cg_build_alignment(assignments=assignments,
                                n_assignments=args.assignments,
                                outfile=out_build_alignment,
                                threshold=args.assignment_threshold)

    # routines to run the HDP builder
    binary        = os.path.join(os.path.abspath(PATH_TO_SIGNALALIGN), "buildHdpUtil")
    kmer_model    = os.path.abspath(PATH_TO_5MER_MODEL)
    hdp_type      = 10  # singleLevelPrior
    grid_start    = 50
    grid_end      = 140
    grid_length   = 1800
    burn_in       = 32 * count_lines_in_build_alignment(out_build_alignment)
    out_hdp       = os.path.join(workdir, "template.singleLevelPrior.nhdp")
    build_command = " ".join([binary,
                              "--verbose",
                              "--oneD",
                              "-p %s" % hdp_type,
                              "-v %s" % out_hdp,
                              "-l %s" % out_build_alignment,
                              "-a %s" % kmer_length_from_model(kmer_model),
                              "-n %s" % args.samples,
                              "-I %s" % burn_in,
                              "-t %s" % args.thinning,
                              "-s %s" % grid_start,
                              "-e %s" % grid_end,
                              "-k %s" % grid_length,
                              "-T %s" % kmer_model,
                              "-g 1",
                              "-r 1",
                              "-j 1",
                              "-y 1",
                              "-i 1",
                              "-u 1"])
    subprocess.check_call(build_command.split(), stdout=sys.stdout, stderr=sys.stderr)

    sys.exit(0)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
