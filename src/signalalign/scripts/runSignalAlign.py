#!/usr/bin/env python
"""Main driver script for running an ionic current-to-sequence alignment on a single machine.
"""

from __future__ import print_function
import sys
import os

# import pysam

from argparse import ArgumentParser
from random import shuffle
from multiprocessing import Process, current_process, Manager

from signalalign.signalAlignment import multithread_signal_alignment, create_signalAlignment_args
from signalalign.utils.sequenceTools import processReferenceFasta
from signalalign.utils.fileHandlers import FolderHandler
from signalalign.utils.bwaWrapper import buildBwaIndex
from signalalign.motif import getDegenerateEnum
from py3helpers.utils import time_it
from signalalign import parseFofn


def signalAlignSourceDir():
    current_file = os.path.abspath(__file__)
    # this is not a great solution, but it enables us to run this in an ide
    scripts_path = "src/signalalign/scripts/runSignalAlign.py"
    bin_path = "bin/runSignalAlign"
    if current_file.endswith(scripts_path):
        current_file = current_file.replace(scripts_path, '')
    elif current_file.endswith(bin_path):
        current_file = current_file.replace(bin_path, '')
    else:
        print("[runSignalAlign] unexpected script location: {}".format(current_file))

    # returns signalAlign root directory or current directory
    return current_file if os.path.isdir(current_file) else "/".join(current_file.split("/")[:-1])


def resolvePath(p):
    if p is None:
        return None
    elif p.startswith("/"):
        return p
    else:
        return os.path.abspath(p)


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--file_directory', '-d', action='store',
                        dest='files_dir', required=True, type=str, default=None,
                        help="directory with MinION fast5 reads to align")
    parser.add_argument('--forward_ref', action='store',
                        dest='forward_ref', required=False, type=str,
                        help="forward reference sequence for SignalAlignment align to, in FASTA")
    parser.add_argument('--backward_ref', action='store',
                        dest='backward_ref', required=False, type=str,
                        help="backward reference sequence for SignalAlignment align to, in FASTA")

    parser.add_argument('--alignment_file', action='store',
                        dest='alignment_file', required=False, type=str,
                        help="SAM or BAM file of alignments if FAST5s do not have fasta info or events")

    parser.add_argument('--bwa_reference', '-r', action='store',
                        dest='bwa_reference', required=True, type=str,
                        help="Reference sequence required for generating guide alignment")
    parser.add_argument('--output_location', '-o', action='store', dest='out',
                        required=True, type=str, default=None,
                        help="directory to put the alignments")
    # optional arguments
    parser.add_argument("--2d", action='store_true', dest="twoD", default=False, help="flag, specify if using 2D reads")
    parser.add_argument('--in_template_hmm', '-T', action='store', dest='in_T_Hmm',
                        required=False, type=str, default=None,
                        help="input HMM for template events, if you don't want the default")
    parser.add_argument('--in_complement_hmm', '-C', action='store', dest='in_C_Hmm',
                        required=False, type=str, default=None,
                        help="input HMM for complement events, if you don't want the default")
    parser.add_argument('--template_hdp', '-tH', action='store', dest='templateHDP', default=None,
                        help="template serialized HDP file")
    parser.add_argument('--complement_hdp', '-cH', action='store', dest='complementHDP', default=None,
                        help="complement serialized HDP file")
    parser.add_argument('--degenerate', '-x', action='store', dest='degenerate', default="variant",
                        help="Specify degenerate nucleotide options: "
                             "variant -> {ACGT}, cytosine2 -> {CE} cytosine3 -> {CEO} adenosine -> {AI}")
    parser.add_argument('--stateMachineType', '-smt', action='store', dest='stateMachineType', type=str,
                        default="threeState", help="decide which model to use, threeState by default")
    parser.add_argument('--file_of_files', '-fofn', action='store', dest='fofn', required=False, type=str, default=None,
                        help="text file containing absolute paths to files to use")
    parser.add_argument('--threshold', '-t', action='store', dest='threshold', type=float, required=False,
                        default=0.01, help="posterior match probability threshold, Default: 0.01")
    parser.add_argument('--diagonalExpansion', '-e', action='store', dest='diag_expansion', type=int,
                        required=False, default=None,
                        help="number of diagonals to expand around each anchor default: 50")
    parser.add_argument('--constraintTrim', '-m', action='store', dest='constraint_trim', type=int,
                        required=False, default=None, help='amount to remove from an anchor constraint')
    parser.add_argument('--target_regions', '-q', action='store', dest='target_regions', type=str,
                        required=False, default=None, help="tab separated table with regions to align to")
    parser.add_argument("--motifs", action="store", dest="motifs", default=None, help="Motif find and replace "
                                                                                      "must be in specific list within "
                                                                                      "a list format. eg: "
                                                                                      "[['CCAGG', 'CEAGG'], "
                                                                                      "['CCTGG', 'CETGG']]")
    parser.add_argument('--ambiguity_positions', '-p', action='store', required=False, default=None,
                        dest='ambiguity_positions', help="Ambiguity positions")
    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False,
                        default=4, type=int, help="number of jobs to run in parallel")
    parser.add_argument('--nb_files', '-n', action='store', dest='nb_files', required=False,
                        default=500, type=int, help="maximum number of reads to align")
    parser.add_argument('--output_format', '-f', action='store', default="full", dest='outFmt',
                        help="output format: full, variantCaller, or assignments. Default: full")
    parser.add_argument('--debug', action='store_true', dest="DEBUG", default=False,
                        help="Disable multiprocessing to gather error messages")
    parser.add_argument('--embed', action='store_true', dest="embed", default=False,
                        help="Embed full output into fast5 file")
    parser.add_argument('--event_table', action='store', dest="event_table", default=False,
                        help="Specify event table")
    parser.add_argument('--force_kmer_event_alignment', action='store_true', dest="perform_kmer_event_alignment",
                        default=None, help="If passed, force SignalAlign to infer kmer to event alignment. "
                                           "Must include alignment_file")
    parser.add_argument('--allow_unsupported_nanopore_read_versions', action='store_false',
                        dest="enforce_supported_versions", default=True,
                        help="Will attempt to complete execution with unsupported nanopore read versions")

    args = parser.parse_args()
    return args


def concat_variant_call_files(path):
    concat_command = "cat {path}/*.tsv > {path}/probs.tsv".format(path=path)
    os.system(concat_command)
    return


def main(args):
    # parse args
    args = parse_args()

    command_line = " ".join(sys.argv[:])
    print(os.getcwd())

    print("Command Line: {cmdLine}\n".format(cmdLine=command_line), file=sys.stderr)
    # get absolute paths to inputs
    args.files_dir = resolvePath(args.files_dir)
    args.forward_reference = resolvePath(args.forward_ref)
    args.backward_reference = resolvePath(args.backward_ref)
    args.out = resolvePath(args.out)
    args.bwa_reference = resolvePath(args.bwa_reference)
    args.in_T_Hmm = resolvePath(args.in_T_Hmm)
    args.in_C_Hmm = resolvePath(args.in_C_Hmm)
    args.templateHDP = resolvePath(args.templateHDP)
    args.complementHDP = resolvePath(args.complementHDP)
    args.fofn = resolvePath(args.fofn)
    args.target_regions = resolvePath(args.target_regions)
    args.ambiguity_positions = resolvePath(args.ambiguity_positions)
    args.alignment_file = resolvePath(args.alignment_file)
    start_message = """
#   Starting Signal Align
#   Aligning files from: {fileDir}
#   Aligning to reference: {reference}
#   Aligning maximum of {nbFiles} files
#   Using model: {model}
#   Using banding: True
#   Aligning to regions in: {regions}
#   Non-default template HMM: {inThmm}
#   Non-default complement HMM: {inChmm}
#   Template HDP: {tHdp}
#   Complement HDP: {cHdp}
    """.format(fileDir=args.files_dir, reference=args.bwa_reference, nbFiles=args.nb_files,
               inThmm=args.in_T_Hmm, inChmm=args.in_C_Hmm, model=args.stateMachineType, regions=args.target_regions,
               tHdp=args.templateHDP, cHdp=args.complementHDP)

    print(start_message, file=sys.stdout)

    if args.files_dir is None and args.fofn is None:
        print("Need to provide directory with .fast5 files of fofn", file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(args.bwa_reference):
        print("Did not find valid reference file, looked for it {here}".format(here=args.bwa_reference),
              file=sys.stderr)
        sys.exit(1)

    # make directory to put temporary files
    if not os.path.isdir(args.out):
        print("Creating output directory: {}".format(args.out), file=sys.stdout)
        os.mkdir(args.out)
    temp_folder = FolderHandler()
    temp_dir_path = temp_folder.open_folder(os.path.join(os.path.abspath(args.out), "tempFiles_alignment"))
    temp_dir_path = resolvePath(temp_dir_path)
    print(args.out)
    print(temp_dir_path)

    # generate reference sequence if not specified
    if not args.forward_reference or not args.backward_reference:
        args.forward_reference, args.backward_reference = processReferenceFasta(fasta=args.bwa_reference,
                                                                                motifs=args.motifs,
                                                                                work_folder=temp_folder,
                                                                                positions_file=args.ambiguity_positions)

    # list of read files
    if args.fofn is not None:
        fast5s = [x for x in parseFofn(args.fofn) if x.endswith(".fast5")]
    else:
        fast5s = ["/".join([args.files_dir, x]) for x in os.listdir(args.files_dir) if x.endswith(".fast5")]

    nb_files = args.nb_files
    if nb_files < len(fast5s):
        shuffle(fast5s)
        fast5s = fast5s[:nb_files]

    # get alignment_args
    alignment_args = {
        "destination": temp_dir_path,
        "stateMachineType": args.stateMachineType,
        "bwa_reference": args.bwa_reference,
        "in_templateHmm": args.in_T_Hmm,
        "in_complementHmm": args.in_C_Hmm,
        "in_templateHdp": args.templateHDP,
        "in_complementHdp": args.complementHDP,
        "output_format": args.outFmt,
        "threshold": args.threshold,
        "diagonal_expansion": args.diag_expansion,
        "constraint_trim": args.constraint_trim,
        "degenerate": getDegenerateEnum(args.degenerate),
        "twoD_chemistry": args.twoD,
        "target_regions": args.target_regions,
        "embed": args.embed,
        "event_table": args.event_table,
        "backward_reference": args.backward_reference,
        "forward_reference": args.forward_reference,
        "alignment_file": args.alignment_file,
        "check_for_temp_file_existance": True,
        "track_memory_usage": False,
        "get_expectations": False,
        "perform_kmer_event_alignment": args.perform_kmer_event_alignment,
        "enforce_supported_versions": args.enforce_supported_versions,
        "path_to_bin": "./bin"
    }

    print("[runSignalAlign]:NOTICE: Got {} files to align".format(len(fast5s)), file=sys.stdout)
    # setup workers for multiprocessing
    os.chdir(signalAlignSourceDir())
    multithread_signal_alignment(alignment_args, fast5s, args.nb_jobs, debug=args.DEBUG)

    print("\n#  signalAlign - finished alignments\n", file=sys.stderr)
    print("\n#  signalAlign - finished alignments\n", file=sys.stdout)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
