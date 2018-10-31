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
from py3helpers.utils import list_dir


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--alignment_file', action='store',
                        dest='alignment_file', required=True, type=str, default=None,
                        help="Sam/Bam file with all alignment data")

    parser.add_argument('--fast5_dir', action='store', default=None,
                        dest='fast5_dir', required=True, type=str,
                        help="Directory of all fast5 files")

    parser.add_argument('--pass_output_dir', action='store', default=None,
                        dest='pass_output_dir', required=True, type=str,
                        help="Location where all pass reads will be moved")

    parser.add_argument('--quality_threshold', action='store', default=7,
                        dest='quality_threshold', required=False, type=float,
                        help="Minimum average base quality threshold. Default = 7")

    args = parser.parse_args()
    return args


def filter_reads(fast5s, alignment_file, quality_threshold=7):
    """Filter fast5 files based on a quality threhsold and if there is an alignment"""
    # loop through fast5s
    fast5_dict = defaultdict()
    # loop through fast5s
    for fast5_path in fast5s:
        assert os.path.exists(fast5_path), "fast5 path does not exist: {}".format(fast5_path)
        f5h = NanoporeRead(fast5_path)
        f5h._initialize_metadata()
        read_name = f5h.read_label
        fast5_dict[read_name] = fast5_path
    print("Created read_id to fast5_path mapping")
    # grab aligned segment
    with closing(pysam.AlignmentFile(alignment_file, 'rb' if alignment_file.endswith("bam") else 'r')) as aln:
        for aligned_segment in aln.fetch():
            if aligned_segment.is_secondary or aligned_segment.is_unmapped \
                    or aligned_segment.is_supplementary or aligned_segment.has_tag("SA"):
                continue
            read_name = aligned_segment.qname.split("_")[0]
            fast5_path = fast5_dict[read_name]
            # get data and sanity check
            if aligned_segment.query_qualities is not None:
                if np.mean(aligned_segment.query_qualities) > quality_threshold:
                    yield fast5_path
            else:
                yield fast5_path



def main():
    args = parse_args()

    fast5s = list_dir(args.fast5_dir, ext='fast5')
    assert len(fast5s) > 0, "Check fast5_dir. No files with fast5 extension found: {}".format(args.fast5_dir)
    # Make output dir if it doesn't exist
    args.pass_output_dir = os.path.abspath(args.pass_output_dir)
    if not os.path.isdir(args.pass_output_dir):
        os.mkdir(args.pass_output_dir)

    best_files = filter_reads(fast5s, args.alignment_file, args.quality_threshold)

    # move passed files
    assert os.path.isdir(args.pass_output_dir), "pass_output_dir does not exist or get created: {}".format(args.pass_output_dir)
    for path in best_files:
        shutil.move(path, os.path.join(args.pass_output_dir, os.path.basename(path)))


if __name__ == "__main__":
    main()
    raise SystemExit
