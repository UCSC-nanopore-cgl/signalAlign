#!/usr/bin/env python
"""Filter reads based on minimum average base quality and alignment"""
########################################################################
# File: filter_reads.py
#  executable: filter_reads.py
#
# Author: Andrew Bailey
# History: 10/15/18 Created
########################################################################

import os
import pysam
import shutil
import numpy as np
from collections import defaultdict
from contextlib import closing
from argparse import ArgumentParser
from signalalign.fast5 import Fast5
from signalalign.nanoporeRead import NanoporeRead
from signalalign.utils.sequenceTools import get_full_nucleotide_read_from_alignment
from py3helpers.utils import list_dir, get_all_sub_directories


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--alignment_file', action='store',
                        dest='alignment_file', required=True, type=str, default=None,
                        help="Sam/Bam file with all alignment data")

    parser.add_argument('--fast5_dir', action='store', default=None,
                        dest='fast5_dir', required=True, type=str,
                        help="Directory of all fast5 files")

    parser.add_argument('--readdb', action='store', default=None,
                        dest='readdb', required=True, type=str,
                        help="Path to readdb file")

    parser.add_argument('--pass_output_dir', action='store', default=None,
                        dest='pass_output_dir', required=True, type=str,
                        help="Location where all pass reads will be moved")

    parser.add_argument('--quality_threshold', action='store', default=7,
                        dest='quality_threshold', required=False, type=float,
                        help="Minimum average base quality threshold. Default = 7")

    parser.add_argument('--recursive', action='store_true', default=False,
                        dest='recursive', required=False,
                        help="Search directory recursively to find fast5 files")

    parser.add_argument('--copy_dir_structure', action='store_true', default=False,
                        dest='copy_dir_structure', required=False,
                        help="Step into directory and copy directory structure for output files")

    parser.add_argument('--trim', action='store', default=False,
                        dest='trim', required=False, type=int,
                        help="Only move as many files which contain a total of some number of bases set by trim. ")

    args = parser.parse_args()
    return args


def parse_readdb(readdb, directories):
    """Parse readdb file

    :param readdb: path to readdb file
    :param directories: path to directories of where the reads are
    """
    assert readdb.endswith("readdb"), "readdb file must end with .readdb: {}".format(readdb)
    with open(readdb, 'r') as fh:
        for line in fh:
            split_line = line.split()
            for dir_path in directories:
                full_path = os.path.join(dir_path, split_line[1])
                if os.path.exists(full_path):
                    yield split_line[0], full_path


def parse_seq_summary(seq_summary, directories):
    """Parse seq_summary file

    :param seq_summary: path to seq_summary file
    :param directories: path to directories of where the reads are
    """
    assert seq_summary.endswith("tsv"), "seq_summary file must end with .tsv: {}".format(seq_summary)
    with open(seq_summary, 'r') as fh:
        for line in fh:
            split_line = line.split()
            for dir_path in directories:
                full_path = os.path.join(dir_path, split_line[0])
                if os.path.exists(full_path):
                    yield split_line[1], full_path


def parse_read_name_map_file(read_map, directories, recursive=False):
    """Parse either a seq summary file or a readdb file"""
    if read_map.endswith("readdb"):
        name_index = 0
        path_index = 1
    else:
        name_index = 1
        path_index = 0

    with open(read_map, 'r') as fh:
        for line in fh:
            split_line = line.split()
            for dir_path in directories:
                if recursive:
                    directories2 = get_all_sub_directories(dir_path)
                    for dir_path2 in directories2:
                        full_path = os.path.join(dir_path2, split_line[path_index])
                        if os.path.exists(full_path):
                            yield split_line[name_index], full_path
                else:
                    full_path = os.path.join(dir_path, split_line[path_index])
                    if os.path.exists(full_path):
                        yield split_line[name_index], full_path


def filter_reads(alignment_file, readdb, read_dirs, quality_threshold=7, recursive=False, trim=False):
    """Filter fast5 files based on a quality threshold and if there is an alignment
    :param alignment_file: bam aligment file
    :param readdb: readdb or sequence summary file
    :param read_dirs: list of directories
    :param quality_threshold: phred quality score min threshold for passing
    :param recursive: search directories recursively for more fast5 dirs
    :param trim: number of bases to analyze
    """
    assert alignment_file.endswith("bam"), "Alignment file must be in BAM format: {}".format(alignment_file)
    # grab aligned segment
    if trim:
        assert isinstance(trim, int), "Trim needs to be an integer: {}".format(trim)
    else:
        trim = np.inf
    n_bases = 0
    with closing(pysam.AlignmentFile(alignment_file, 'rb')) as bamfile:
        name_indexed = pysam.IndexedReads(bamfile)
        name_indexed.build()
        for name, fast5 in parse_read_name_map_file(readdb, read_dirs, recursive=recursive):
            try:
                iterator = name_indexed.find(name)
                for aligned_segment in iterator:
                    if aligned_segment.is_secondary or aligned_segment.is_unmapped \
                            or aligned_segment.is_supplementary or aligned_segment.has_tag("SA"):
                        continue
                    # get data and sanity check
                    if aligned_segment.query_qualities is not None:
                        if np.mean(aligned_segment.query_qualities) < quality_threshold:
                            continue
                    if trim < n_bases:
                        break
                    n_bases += aligned_segment.query_length
                    yield fast5, aligned_segment

            except KeyError:
                print("Found no alignments for {}".format(fast5))


def make_readdb(alignment_file, read_dirs, readdb_out):
    """Create a readdb generator from an alignment file"""
    assert out_path.endswith("readdb"), "Readdb output file nees to end with readdb. {}".format(readdb_out)
    assert alignment_file.endswith("bam"), "Alignment file must be in BAM format: {}".format(alignment_file)
    # grab aligned segment
    for read_dir in read_dirs:
        fast5s = list_dir(read_dir, ext="fast5")
        for fast5_path in fast5s:
            assert os.path.exists(fast5_path), "fast5 path does not exist: {}".format(fast5_path)
            f5h = NanoporeRead(fast5_path)
            f5h._initialize_metadata()
            read_name = f5h.read_label
            fast5_dict[read_name] = fast5_path

    with closing(pysam.AlignmentFile(alignment_file, 'rb' if alignment_file.endswith("bam") else 'r')) as aln:
        for aligned_segment in aln.fetch():
            try:
                read_name = aligned_segment.qname.split("_")[0]
                fast5_path = fast5_dict[read_name]
                yield aligned_segment.qname, fast5_path
            except KeyError:
                print("{} not found in fast5s".format(read_name))


def filter_reads_to_string_wrapper(filter_reads_generator):
    """Wrap filter reads in order to convert the aligned segment into a string so it can be pickled"""
    for fast5, aligned_segment in filter_reads_generator:
        yield fast5, aligned_segment.to_string()


def main():
    args = parse_args()

    fast5s = list_dir(args.fast5_dir, ext='fast5')
    assert len(fast5s) > 0, "Check fast5_dir. No files with fast5 extension found: {}".format(args.fast5_dir)
    # Make output dir if it doesn't exist
    args.pass_output_dir = os.path.abspath(args.pass_output_dir)
    if not os.path.isdir(args.pass_output_dir):
        os.mkdir(args.pass_output_dir)

    if copy_dir_structure:
        for sub_dir in get_all_sub_directories(args.fast5_dir):
            # make dir if it doesnt exist
            out_dir = os.path.join(args.pass_output_dir, os.path.basename(sub_dir))
            if not os.path.isdir(out_dir):
                os.mkdir(out_dir)

            best_files = filter_reads(args.alignment_file, args.readdb, [sub_dir], args.quality_threshold, trim=args.trim)
            assert os.path.isdir(out_dir), "out_dir does not exist or get created: {}".format(out_dir)
            for path, _ in best_files:
                shutil.move(path, os.path.join(out_dir, os.path.basename(path)))

    else:
        best_files = filter_reads(args.alignment_file, args.readdb, [args.fast5_dir], args.quality_threshold,
                                  recursive=args.recursive, trim=args.trim)
        # move passed files
        assert os.path.isdir(args.pass_output_dir), "pass_output_dir does not exist or get created: {}".format(args.pass_output_dir)
        for path, _ in best_files:
            shutil.move(path, os.path.join(args.pass_output_dir, os.path.basename(path)))


if __name__ == "__main__":
    main()
    raise SystemExit
