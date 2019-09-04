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
from signalalign.utils import multithread
from signalalign.nanoporeRead import NanoporeRead
from signalalign.utils.sequenceTools import get_full_nucleotide_read_from_alignment
from py3helpers.utils import list_dir, get_all_sub_directories, merge_dicts, merge_lists


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

    parser.add_argument('--jobs', action='store', default=1,
                        dest='jobs', required=False, type=int,
                        help="Number of jobs to start if copy_dir_structure is set")

    parser.add_argument('--debug', action='store_true', default=False,
                        dest='debug', required=False,
                        help="Will run copy_dir_structure with only one job and fail if errors arise")

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
    """Parse either a seq summary file or a readdb file
    :param read_map: either a readdb file or sequencing summary file
    :param directories: check all the directories for the fast5 path
    :param recursive: boolean option to check the sub directories of input directories
    """
    if read_map.endswith("readdb"):
        name_index = 0
        path_index = 1
    else:
        name_index = 1
        path_index = 0
    for dir_path in directories:
        assert os.path.isdir(dir_path), "Path provided does not exist or is not a directory: {}".format(dir_path)

    with open(read_map, 'r') as fh:
        for line in fh:
            split_line = line.split()
            if len(split_line) < 2:
                continue
            for dir_path in directories:
                if recursive:
                    if os.path.exists(split_line[path_index]):
                        yield split_line[name_index], os.path.abspath(split_line[path_index])
                    else:
                        directories2 = get_all_sub_directories(dir_path)
                        for dir_path2 in directories2:
                            full_path = os.path.join(dir_path2, split_line[path_index])
                            if os.path.exists(full_path):
                                yield split_line[name_index], os.path.abspath(full_path)
                else:
                    if os.path.exists(split_line[path_index]):
                        yield split_line[name_index], os.path.abspath(split_line[path_index])
                    else:
                        full_path = os.path.join(dir_path, split_line[path_index])
                        if os.path.exists(full_path):
                            yield split_line[name_index], os.path.abspath(full_path)


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
    n_files = 0
    with closing(pysam.AlignmentFile(alignment_file, 'rb')) as bamfile:
        name_indexed = pysam.IndexedReads(bamfile)
        name_indexed.build()
        for name, fast5 in parse_read_name_map_file(readdb, read_dirs, recursive=recursive):
            try:
                if trim < n_bases:
                    print("Filtered {} files for {} bases".format(n_files, n_bases))
                    break

                iterator = name_indexed.find(name)
                for aligned_segment in iterator:
                    if aligned_segment.is_secondary or aligned_segment.is_unmapped \
                            or aligned_segment.is_supplementary or aligned_segment.has_tag("SA"):
                        continue
                    # get data and sanity check
                    if aligned_segment.query_qualities is not None:
                        if np.mean(aligned_segment.query_qualities) < quality_threshold:
                            continue
                    n_files += 1
                    n_bases += aligned_segment.query_length
                    yield fast5, aligned_segment

            except KeyError:
                print("Found no alignments for {}".format(fast5))


def filter_reads_to_string_wrapper(filter_reads_generator):
    """Wrap filter reads in order to convert the aligned segment into a string so it can be pickled
    :param filter_reads_generator: a filter_reads generator object
    """
    for fast5, aligned_segment in filter_reads_generator:
        yield fast5, aligned_segment.to_string()


def filter_read_service2(work_queue, done_queue, service_name="filter_reads_service2"):
    """
    Service used by the multithread module in signal align to filter reads from a large number of directories
    :param work_queue: arguments to be done
    :param done_queue: errors and returns to be put
    :param service_name: name of the service
    """
    # prep
    total_handled = 0
    failure_count = 0
    mem_usages = list()

    # catch overall exceptions
    try:
        for f in iter(work_queue.get, 'STOP'):
            # catch exceptions on each element
            try:
                reads = filter_read_wrapper(**f)
                done_queue.put(reads)
            except Exception as e:
                # get error and log it
                message = "{}:{}".format(type(e), str(e))
                error = "{} '{}' failed with: {}".format(service_name, multithread.current_process().name, message)
                print("[{}] ".format(service_name) + error)
                done_queue.put(error)
                failure_count += 1

            # increment total handling
            total_handled += 1

    except Exception as e:
        # get error and log it
        message = "{}:{}".format(type(e), str(e))
        error = "{} '{}' critically failed with: {}".format(service_name, multithread.current_process().name, message)
        print("[{}] ".format(service_name) + error)
        done_queue.put(error)

    finally:
        # logging and final reporting
        print("[%s] '%s' completed %d calls with %d failures"
              % (service_name, multithread.current_process().name, total_handled, failure_count))
        done_queue.put("{}:{}".format(multithread.TOTAL_KEY, total_handled))
        done_queue.put("{}:{}".format(multithread.FAILURE_KEY, failure_count))
        if len(mem_usages) > 0:
            done_queue.put("{}:{}".format(multithread.MEM_USAGE_KEY, ",".join(map(str, mem_usages))))


def filter_read_service(work_queue, done_queue, service_name="filter_reads"):
    """
    Service used by the multithread module in signal align to filter reads from a large number of directories
    :param work_queue: arguments to be done
    :param done_queue: errors and returns to be put
    :param service_name: name of the service
    """
    # prep
    total_handled = 0
    failure_count = 0
    mem_usages = list()

    # catch overall exceptions
    try:
        for f in iter(work_queue.get, 'STOP'):
            # catch exceptions on each element
            try:
                filter_read_wrapper_for_making_dir(**f)
            except Exception as e:
                # get error and log it
                message = "{}:{}".format(type(e), str(e))
                error = "{} '{}' failed with: {}".format(service_name, multithread.current_process().name, message)
                print("[{}] ".format(service_name) + error)
                done_queue.put(error)
                failure_count += 1

            # increment total handling
            total_handled += 1

    except Exception as e:
        # get error and log it
        message = "{}:{}".format(type(e), str(e))
        error = "{} '{}' critically failed with: {}".format(service_name, multithread.current_process().name, message)
        print("[{}] ".format(service_name) + error)
        done_queue.put(error)

    finally:
        # logging and final reporting
        print("[%s] '%s' completed %d calls with %d failures"
              % (service_name, multithread.current_process().name, total_handled, failure_count))
        done_queue.put("{}:{}".format(multithread.TOTAL_KEY, total_handled))
        done_queue.put("{}:{}".format(multithread.FAILURE_KEY, failure_count))
        if len(mem_usages) > 0:
            done_queue.put("{}:{}".format(multithread.MEM_USAGE_KEY, ",".join(map(str, mem_usages))))


def filter_read_wrapper_for_making_dir(in_dir, out_dir, alignment_file, readdb, quality_threshold=7, trim=False):
    """Helper function for multiprocessing all the different sub directories
    :param in_dir: input directory to search for fast5s
    :param out_dir: output directory to place passing fast5s
    :param alignment_file: alignment file to look into
    :param readdb: readdb or sequencing summary file
    :param quality_threshold: q score threshold
    :param trim: option to trim for x number of bases"""
    best_files = filter_reads(alignment_file, readdb, [in_dir], quality_threshold, trim=trim)
    assert os.path.isdir(out_dir), "out_dir does not exist or get created: {}".format(out_dir)
    for path, _ in best_files:
        shutil.move(path, os.path.join(out_dir, os.path.basename(path)))


def filter_read_wrapper(in_dir, alignment_file, readdb, quality_threshold=7, trim=False):
    """Helper function for multiprocessing all the different sub directories
    :param in_dir: input directory to search for fast5s
    :param out_dir: output directory to place passing fast5s
    :param alignment_file: alignment file to look into
    :param readdb: readdb or sequencing summary file
    :param quality_threshold: q score threshold
    :param trim: option to trim for x number of bases"""
    best_files = [(fast5, cigar_string) for fast5, cigar_string in
                  filter_reads_to_string_wrapper(filter_reads(alignment_file, readdb, [in_dir],
                                                              quality_threshold, trim=trim))]
    return best_files


def create_new_directories_for_filter_reads(in_dir, out_dir):
    """Copy directory structure and return the input directory and output directory for each interal dir
    :param in_dir: top input directory to
    :param out_dir: top of the output directories
    """
    for sub_in_dir in get_all_sub_directories(in_dir):
        # make dir if it doesnt exist
        sub_out_dir = os.path.join(out_dir, os.path.basename(sub_in_dir))
        if not os.path.isdir(sub_out_dir):
            os.mkdir(sub_out_dir)
        yield sub_in_dir, sub_out_dir


def multiprocess_move_and_filter_reads(in_dir, out_dir, alignment_file, readdb, trim=False,
                                       quality_threshold=7, worker_count=1, debug=False):
    """Multiprocess for filtering reads
    :param in_dir: input directory with subdirectories assumed to have fast5s in them
    :param out_dir: head output directory
    :param alignment_file: bam file
    :param readdb: readdb or sequence summary file
    :param trim: option to trim for x number of bases
    :param quality_threshold: quality threshold
    :param worker_count: number of workers to use
    :param debug: boolean option which will only use one process in order to fail if an error arises
    :return: True
    """
    assert alignment_file.endswith("bam"), "Alignment file must be in BAM format: {}".format(alignment_file)
    # grab aligned segment
    if debug:
        for sub_in_dir, sub_out_dir in create_new_directories_for_filter_reads(in_dir, out_dir):
            best_files = filter_reads(alignment_file, readdb, [sub_in_dir],
                                      quality_threshold=quality_threshold, trim=trim)
            assert os.path.isdir(sub_out_dir), "sub_out_dir does not exist or get created: {}".format(sub_out_dir)
            for path, _ in best_files:
                shutil.move(path, os.path.join(sub_out_dir, os.path.basename(path)))
    else:
        filter_reads_args = {"readdb": readdb, "alignment_file": alignment_file,
                             "quality_threshold": quality_threshold, "trim": trim}
        total, failure, messages, output = multithread.run_service2(
            filter_read_service, create_new_directories_for_filter_reads(in_dir, out_dir),
            filter_reads_args, ["in_dir", "out_dir"], worker_count)
    return True


def multiprocess_filter_reads(in_dir, alignment_file, readdb, trim=False,
                              quality_threshold=7, worker_count=1, debug=False):
    """Multiprocess for filtering reads but dont move the files
    :param in_dir: input directory with subdirectories assumed to have fast5s in them
    :param alignment_file: bam file
    :param readdb: readdb or sequence summary file
    :param trim: option to trim for x number of bases
    :param quality_threshold: quality threshold
    :param worker_count: number of workers to use
    :param debug: boolean option which will only use one process in order to fail if an error arises
    :return: True
    """
    assert alignment_file.endswith("bam"), "Alignment file must be in BAM format: {}".format(alignment_file)
    # grab aligned segment
    if debug:
        best_files = []
        for sub_in_dir in get_all_sub_directories(in_dir):
            best_files.extend(filter_reads(alignment_file, readdb, [sub_in_dir],
                                           quality_threshold=quality_threshold, trim=trim))
    else:
        filter_reads_args = {"readdb": readdb, "alignment_file": alignment_file,
                             "quality_threshold": quality_threshold, "trim": trim}
        total, failure, messages, output = multithread.run_service2(
            filter_read_service2, get_all_sub_directories(in_dir),
            filter_reads_args, ["in_dir"], worker_count)
        best_files = merge_lists(output)
    return best_files


def find_fast5s_from_ids_readdb(readdb, read_ids, read_dirs, recursive=False):
    """Find the corresponding fast5 files given readids"""
    for name, fast5 in parse_read_name_map_file(readdb, read_dirs, recursive=recursive):
        if name.split("_")[0] in read_ids:
            yield name, fast5


def write_readdb(list_of_names_and_fast5s, out_path):
    """Write a readdb file given a list of pairs of names and fast5 files"""
    with open(out_path, "w") as fh:
        for pair in list_of_names_and_fast5s:
            fh.write(pair[0]+"\t"+os.path.basename(pair[1])+"\n")
    return 0


def copy_files_from_readdb(readdb, directories, new_dir, recursive=False):
    """Copy files from a readdb file to another location
    :param readdb: readdb file
    :param directories: directories to search for fast5 files
    :param new_dir: new directory to copy files into
    :param recursive: recursive option for looking for files
    """
    assert os.path.isdir(new_dir), "New directory must exist already"
    n_files_copied = 0
    for name, path in parse_read_name_map_file(readdb, directories, recursive=recursive):
        new_path = os.path.join(new_dir, os.path.basename(path))
        shutil.copy(path, new_path)
        n_files_copied +=1
    return n_files_copied


def main():
    args = parse_args()
    # Make output dir if it doesn't exist
    args.pass_output_dir = os.path.abspath(args.pass_output_dir)
    if not os.path.isdir(args.pass_output_dir):
        os.mkdir(args.pass_output_dir)

    if args.copy_dir_structure:
        multiprocess_move_and_filter_reads(args.fast5_dir, args.pass_output_dir, args.alignment_file, args.readdb,
                                           trim=False,
                                           quality_threshold=args.quality_threshold, worker_count=args.jobs,
                                           debug=args.debug)
        # for sub_dir in get_all_sub_directories(args.fast5_dir):
        #     # make dir if it doesnt exist
        #     out_dir = os.path.join(args.pass_output_dir, os.path.basename(sub_dir))
        #     if not os.path.isdir(out_dir):
        #         os.mkdir(out_dir)
        #
        #     best_files = filter_reads(args.alignment_file, args.readdb, [sub_dir], args.quality_threshold, trim=args.trim)
        #     assert os.path.isdir(out_dir), "out_dir does not exist or get created: {}".format(out_dir)
        #     for path, _ in best_files:
        #         shutil.move(path, os.path.join(out_dir, os.path.basename(path)))

    else:
        best_files = filter_reads(args.alignment_file, args.readdb, [args.fast5_dir], args.quality_threshold,
                                  recursive=args.recursive, trim=args.trim)
        # move passed files
        assert os.path.isdir(args.pass_output_dir), "pass_output_dir does not exist or get created: {}".format(
            args.pass_output_dir)
        for path, _ in best_files:
            shutil.move(path, os.path.join(args.pass_output_dir, os.path.basename(path)))


if __name__ == "__main__":
    main()
    raise SystemExit
