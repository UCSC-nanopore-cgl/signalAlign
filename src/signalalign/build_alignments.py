#!/usr/bin/env python
"""Class and methods to deal with alignments and generate alignment tables for HDP training"""
########################################################################
# File: build_alignments.py
#  executable: build_alignments.py
#
# Author: Andrew Bailey
# History: Created 05/02/19
########################################################################

import os
import heapq
import queue
import shutil
from py3helpers.utils import all_string_permutations, list_dir
from py3helpers.multiprocess import *
from signalalign.hiddenMarkovModel import HmmModel


def add_to_queue(queue):
    for x in range(10):
        queue.put(x)
    queue.put('STOP')


def get_from_queue(queue, n_threads):
    all_data = []
    n_stops = 0
    while True:
        x = queue.get()
        if x == "STOP":
            n_stops += 1
            if n_stops == n_threads:
                return all_data
        else:
            all_data.append(x)


def get_nlargest_dataframe(dataframe, topn=100):
    data = [(x[12], x) for x in dataframe.itertuples(index=False, name=None)]
    return heapq.nlargest(topn, data)


def get_nlargest_queue(my_queue, topn=100):
    data = []
    try:
        while True:
            data.append(my_queue.get(False))
    except queue.Empty:
        pass
    return heapq.nlargest(topn, data)


def get_nlargest_alignment_queue(my_queue, topn=100):
    data = []
    try:
        while True:
            data.append(my_queue.get(False))
    except queue.Empty:
        pass
    return heapq.nlargest(topn, data, key=lambda e: e[3])


def get_assignment_kmer_tables(file_path, kmers, min_probability, full):
    if full:
        data = parse_alignment_file(file_path)
    else:
        data = parse_assignment_file(file_path)
    data = data.loc[data['prob'] >= min_probability]
    all_kmers = []
    for k in kmers:
        all_kmers.append(data.loc[data.kmer == k])
    return all_kmers


def alignment_file_to_queues(assignment_file_path, my_queues, min_prob=0.0):
    """Put each kmer into a queue which """
    with open(assignment_file_path, "r") as fh:
        for line in fh:
            split_line = line.split()
            if float(split_line[12]) >= min_prob:
                outline = [split_line[9], split_line[4], float(split_line[13]), float(split_line[12])]
                try:
                    if split_line[4] == "c":
                        my_queues[1].put(outline, False)
                    else:
                        my_queues[0].put(outline, False)
                except queue.Full:
                    pass


def assignment_file_to_queues(assignment_file_path, my_queues, min_prob=0.0):
    """Put each kmer into a queue which """
    with open(assignment_file_path, "r") as fh:
        for line in fh:
            split_line = line.split()
            if float(split_line[3]) >= min_prob:
                outline = [split_line[0], split_line[1], float(split_line[2]), float(split_line[3])]
                try:
                    if split_line[1] == "c":
                        my_queues[1].put(outline, False)
                    else:
                        my_queues[0].put(outline, False)
                except queue.Full:
                    pass


def make_kmer_directories(dir_path, alphabet, kmer_length, complement=False):
    """Make the kmer directories where all the kmers will be written
    :param dir_path: path to directory
    :param alphabet: kmer alphabet
    :param kmer_length: length of kmer
    :param complement: boolean option to create complement kmers
    """
    assert os.path.isdir(dir_path), "dir_path is not a directory. {}".format(dir_path)
    template_dirs = []
    complment_dirs = []

    for kmer in all_string_permutations(alphabet, length=kmer_length):
        template_path = os.path.join(dir_path, kmer)
        os.mkdir(template_path)
        template_dirs.append(template_path)
        if complement:
            complement_path = os.path.join(dir_path, kmer+"_c")
            os.mkdir(complement_path)
            complment_dirs.append(complement_path)
    template_dirs.extend(complment_dirs)
    return template_dirs


def split_assignment_file(assignment_file_path, my_dirs, alphabet, kmer_length, alphabet_size, min_prob=0.0,
                          alignment=False, remove=False):
    """Split kmers and write to new files
    :param assignment_file_path: path to assignment file
    :param my_dirs: list of directories to write
    :param alphabet: kmer alphabet
    :param kmer_length: kmer length
    :param alphabet_size: size of alphabet
    :param min_prob: minimum probability
    :param alignment: if set will select columns from an alignment file with 14 columns
    """
    basename = os.path.basename(assignment_file_path)
    data = [[] for _ in range((alphabet_size**kmer_length)*2)]
    if alignment:
        prob_index = 12
        kmer_index = 9
        strand_index = 4
        mean_index = 13
    else:
        kmer_index = 0
        strand_index = 1
        mean_index = 2
        prob_index = 3

    with open(assignment_file_path, "r") as fh:
        for line in fh:
            split_line = line.split()
            if float(split_line[prob_index]) >= min_prob:
                outline = [split_line[kmer_index], split_line[strand_index],
                           split_line[mean_index], split_line[prob_index]]
                k_index = HmmModel._get_kmer_index(split_line[kmer_index], alphabet, kmer_length, alphabet_size)

                if split_line[strand_index] == "c":
                    data[k_index+(alphabet_size**kmer_length)].append(outline)
                else:
                    data[k_index].append(outline)

    for directory, kmer_data in zip(my_dirs, data):
        if len(kmer_data) > 0:
            with open(os.path.join(directory, basename), "w") as fh2:
                for line in kmer_data:
                    fh2.write("\t".join(line)+"\n")
    if remove:
        os.remove(assignment_file_path)
    return True


def multiprocess_split_sa_tsv_file(list_of_assignment_paths, my_dirs, alphabet, kmer_length,
                                       min_prob=0.0, alignment=False, remove=False, worker_count=1):

    # Multiprocess reading in files
    extra_args = {"my_dirs": my_dirs,
                  "alphabet": alphabet,
                  "kmer_length": kmer_length,
                  "alphabet_size": len(alphabet),
                  "min_prob": min_prob,
                  "alignment": alignment,
                  "remove": remove}
    service = BasicService(split_assignment_file, service_name="multiprocess_split_sa_tsv_file")
    total, failure, messages, output = run_service(service.run, list_of_assignment_paths,
                                                   extra_args, ["assignment_file_path"], worker_count)
    return True


def concatenate_files(file_paths, output_file_path, remove_files=False):
    """
    Concatenate files (efficiently) given a list of file paths to a single file
    :param file_paths: List of file path
    :param output_file_path: Output file path name
    :param remove_files: remove files that were concatenated together
    :return: None
    From Ryan Lorig-Roach
    """
    with open(output_file_path, 'wb') as out_file:
        for file_path in file_paths:
            with open(file_path, 'rb') as in_file:
                # 100MB per writing chunk to avoid reading big file into memory.
                shutil.copyfileobj(in_file, out_file, 1024*1024*100)
            if remove_files:
                os.remove(file_path)


def get_top_kmers_from_directory(kmer_dir, output_dir, n, random=False):
    """Get the top n kmers from a directory of kmer tables
    :param kmer_dir: path to directory with kmer files
    :param output_dir: path to output dir
    :param n: number of kmers to collect
    :param random: boolean option to select random kmers instead of top n
    """
    output_file = os.path.join(kmer_dir, "all_kmers.tsv")
    kmer = os.path.basename(kmer_dir)
    concatenate_files(list_dir(kmer_dir, ext="tsv"), output_file)
    with open(output_file, "r") as fh:
        data = [x.split() for x in fh.readlines()]

    # sort data
    if len(data) < n:
        print("Not enough events for kmer {}. {}".format(kmer, len(data)))
        largest = data
    else:
        print("Getting events for kmer {}. {}".format(kmer, n))

        largest = heapq.nlargest(n, data, key=lambda e: e[3])

    # write data
    output_file = os.path.join(output_dir, kmer+".tsv")
    with open(output_file, "w") as fh2:
        for line in largest:
            fh2.write("\t".join(line)+"\n")

    shutil.rmtree(kmer_dir)


def multiprocess_get_top_kmers_from_directory(kmer_dirs, output_dir, n, random=False, worker_count=1):
    """Multiprocess get_top_kmers_from_directory"""
    # Multiprocess reading in files
    extra_args = {"output_dir": output_dir,
                  "n": n,
                  "random": random}
    service = BasicService(get_top_kmers_from_directory, service_name="multiprocess_get_top_kmers_from_directory")
    total, failure, messages, output = run_service(service.run, kmer_dirs,
                                                   extra_args, ["kmer_dir"], worker_count)
    return True


def generate_top_n_kmers_from_sa_output(assignment_files, working_dir, output_file, n, alphabet="ACGT",
                                        kmer_len=5, min_prob=0.8, worker_count=1, random=False, complement=False,
                                        remove=False, alignment=False):
    """Get the most probable n events for each kmer"""
    kmer_dirs = make_kmer_directories(working_dir, alphabet, kmer_len, complement=complement)
    multiprocess_split_sa_tsv_file(assignment_files, kmer_dirs, alphabet, kmer_len, min_prob=min_prob, remove=remove,
                                   worker_count=worker_count, alignment=alignment)

    multiprocess_get_top_kmers_from_directory(kmer_dirs, working_dir, n, random=random,
                                              worker_count=worker_count)
    concatenate_files(list_dir(working_dir, ext="tsv"), output_file, remove_files=True)
    return output_file