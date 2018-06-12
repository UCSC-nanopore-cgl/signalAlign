#!/usr/bin/env python3
"""Train HMM transitions using the BuildModels.py and
"""
########################################################################
# File: TrainSignalAlign.py
#  executable: TrainSignalAlign.py
#
# Author: Andrew Bailey
# History: 5/21/18 Created
########################################################################

import os
import sys
import glob
import pandas as pd
import numpy as np
import string
from argparse import ArgumentParser
from subprocess import Popen
from itertools import product
from random import shuffle
from signalalign.utils.commonFunctions import get_first_seq, make_motif_file, get_all_sequence_kmers, make_CCWGG_positions_file, \
    find_ccwgg_motifs
from signalalign.train.BuildModels import *


def main():
    def parse_args():
        parser = ArgumentParser(description=__doc__)
        parser.add_argument("-r", action="store", dest="reference", required=True)
        parser.add_argument("-pcr", action="store", dest="pcr_reads", required=True)
        parser.add_argument("-gen", action="store", dest="genomic_reads", required=True)
        parser.add_argument('--in_template_hmm', '-T', action='store', dest='in_T_Hmm', required=True, type=str)
        parser.add_argument('--in_complement_hmm', '-C', action='store', dest='in_C_Hmm', required=True, type=str)
        parser.add_argument("-o", action="store", dest="outpath", required=True)
        parser.add_argument("-x", action="store", dest="degenerate", required=True)
        parser.add_argument('--positions', action='append', dest='positions_file', required=False, default=None)
        parser.add_argument('--motif', action='append', dest='motif_file', required=False, default=None)
        parser.add_argument('--bulk', action='store_true', dest='bulk', required=False, default=False)
        parser.add_argument('--hdp_type', action='store', dest='hdp_type', required=False, default="multiset")
        parser.add_argument("-j", action="store", dest="jobs", required=False, default=4, type=int)
        parser.add_argument("-i", action="store", dest="iterations", required=False, type=int, default=20)
        parser.add_argument("-a", action="store", dest="batch", required=False, type=int, default=15000)
        parser.add_argument("-s", action="store", dest="assignments", required=False, type=int, default=30)
        parser.add_argument("-c", action="store", dest="methyl_assignments", required=False, type=int, default=200)
        parser.add_argument("-n", action="store", dest="n_aligns", required=False, type=int, default=1000)
        parser.add_argument("-e", action="store", dest="n_test_alns", required=False, type=int, default=1000)
        parser.add_argument("-t", action="store", dest="assignment_threshold", required=False, type=float, default=0.8)
        parser.add_argument("-g", action="store", dest="samples", required=False, type=int, default=15000)
        parser.add_argument("--hdp_em", action="store", dest="HDP_EM", required=False, type=int, default=None)
        parser.add_argument("--split", action="store", dest="split", required=False, type=float, default=0.5)
        args = parser.parse_args()
        return args

    command_line = " ".join(sys.argv[:])
    print("Command Line: {cmdLine}\n".format(cmdLine=command_line), file=sys.stderr)

    args = parse_args()

    # decide which type of HDP to carry through the pipeline
    if args.hdp_type == "multiset":
        HDP_type = "multisetPriorEcoli"
    else:
        HDP_type = "singleLevelPriorEcoli"

    working_path = os.path.abspath(args.outpath)

    reference_location = os.path.abspath(args.reference)

    # make the positions and motif file
    if args.positions_file is not None and args.motif_file is not None:
        assert len(args.positions_file) == 2 and len(args.motif_file) == 2, "need to give training and testing " \
                                                                            "positions/motif files"
        for i in range(2):
            assert os.path.exists(args.positions_file[i]), "Didn't find positions file, looked " \
                                                           "{}".format(args.positions_file)
            assert os.path.exists(args.motif_file[i]), "Didn't find motif file, looked {}".format(args.motif_file)
        positions_file = args.positions_file[0]
        motif_file = args.motif_file[0]
        test_positions = args.positions_file[1]
        test_motifs = args.motif_file[1]
    else:
        # make the positions file
        positions_file = make_positions_file(fasta=args.reference,
                                             degenerate=args.degenerate,
                                             outfile=working_path + "/{}_positions.positions".format(args.degenerate))

        # make the motif file
        motif_file = make_gatc_or_ccwgg_motif_file(fasta=args.reference,
                                                   degenerate=args.degenerate,
                                                   outfile=working_path + "/{}_target.target".format(args.degenerate))
        test_positions = positions_file
        test_motifs = motif_file

    # make the fofns for training and testing
    pcr_fofns, gen_fofns = train_test_split_fofn(pcr_reads_dir=args.pcr_reads,
                                                 genomic_reads_dir=args.genomic_reads,
                                                 working_directory=working_path,
                                                 split=args.split)

    # train the transitions
    models = train_model_transitions(fasta=reference_location,
                                     pcr_fofn=pcr_fofns[0],
                                     genomic_fofn=gen_fofns[0],
                                     degenerate=args.degenerate,
                                     jobs=args.jobs,
                                     positions_file=positions_file,
                                     iterations=args.iterations,
                                     batch_size=args.batch,
                                     outpath=working_path,
                                     hdp_type=HDP_type,
                                     t_model=os.path.abspath(args.in_T_Hmm),
                                     c_model=os.path.abspath(args.in_C_Hmm))
    # do the initial alignments
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
    sys.exit(main())
