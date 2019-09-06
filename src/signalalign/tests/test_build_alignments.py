#!/usr/bin/env python
"""Test trainModels.py"""
########################################################################
# File: test_trainModels.py
#  executable: test_trainModels.py
#
# Author: Andrew Bailey
# History: 5/21/18 Created
########################################################################


import unittest
import tempfile
import numpy as np
from signalalign.build_alignments import *
from signalalign.train.trainModels import get_kmers, multiprocess_make_kmer_assignment_tables
from py3helpers.utils import time_it


class TrainSignalAlignTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TrainSignalAlignTest, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])
        cls.reference = os.path.join(cls.HOME, "tests/test_sequences/pUC19_SspI_Zymo.fa")
        cls.ecoli_reference = os.path.join(cls.HOME, "tests/test_sequences/E.coli_K12.fasta")
        cls.ecoli_bam = os.path.join(cls.HOME, "tests/minion_test_reads/canonical_ecoli_R9/canonical_ecoli.bam")
        cls.ecoli_readdb = os.path.join(cls.HOME, "tests/minion_test_reads/canonical_ecoli_R9/canonical_ecoli.readdb")
        cls.fast5_dir = os.path.join(cls.HOME, "tests/minion_test_reads/canonical_ecoli_R9")
        cls.files = [
            "miten_PC_20160820_FNFAD20259_MN17223_mux_scan_AMS_158_R9_WGA_Ecoli_08_20_16_83098_ch138_read23_"
            "strand.fast5",
            "miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch101_"
            "read456_strand.fast5",
            "miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch101_"
            "read544_strand1.fast5",
            "miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch103_"
            "read333_strand1.fast5"]
        cls.fast5_paths = [os.path.join(cls.fast5_dir, f) for f in os.listdir(cls.fast5_dir)
                           if os.path.isfile(os.path.join(cls.fast5_dir, f))]
        cls.alignments_dir = os.path.join(cls.HOME, "tests/test_alignments/pUC_6mer_tempFiles_alignment")
        cls.alignments_path = os.path.join(cls.HOME,
                                           "tests/test_alignments/pUC_6mer_tempFiles_alignment/"
                                           "07b1cae8-9b48-426f-a2a5-b07437d5a58e.sm.backward.tsv")
        cls.assignment_path = os.path.join(cls.HOME,
                                           "tests/test_assignment_files/"
                                           "d6160b0b-a35e-43b5-947f-adaa1abade28.sm.assignments.tsv")

    def test_add_to_queue(self):
        work_queue = Manager().Queue()
        worker_count = 2
        for w in range(worker_count):
            p = Process(target=add_to_queue, args=(work_queue,), daemon=True)
            p.start()
        data = get_from_queue(work_queue, worker_count)
        # print(data)
        self.assertSequenceEqual(sorted(list(range(10)) + list(range(10))), sorted(data))

    def test_alignment_file_to_queues(self):
        max_size = 10
        work_queues = [Manager().Queue(max_size) for _ in range(2)]

        alignment_file_to_queues(self.alignments_path, work_queues, min_prob=0.9)
        first_one = work_queues[0].get()
        self.assertSequenceEqual(first_one, ["TGAAAA", "t", 75.375476, 0.999576])

    def test_assignment_file_to_queues(self):
        max_size = 10
        work_queues = [Manager().Queue(max_size) for _ in range(2)]

        assignment_file_to_queues(self.assignment_path, work_queues, min_prob=0.9)

        first_one = work_queues[0].get()
        self.assertSequenceEqual(first_one, ["GCCTTA", "t", 83.709275, 1.000000])

    def test_get_nlargest_queue(self):
        work_queue = Manager().Queue()
        all_data = []
        for x in np.random.randint(100, size=100):
            work_queue.put(x)
            all_data.append(x)

        data = get_nlargest_queue(work_queue, topn=100)
        self.assertSequenceEqual(data, sorted(all_data)[::-1])

    def test_get_nlargest_alignment_queue(self):
        work_queue = Manager().Queue()
        all_data = []
        for x in np.random.randint(100, size=100):
            data = ["GCCTTA", "t", "83.709275", x]
            work_queue.put(data)
            all_data.append(data)

        data = get_nlargest_alignment_queue(work_queue, topn=100)
        self.assertSequenceEqual(data, sorted(all_data)[::-1])

    def test_make_kmer_directories(self):
        with tempfile.TemporaryDirectory() as temdir:
            dirs = make_kmer_directories(temdir, "ACGT", 6, complement=False)
            for new_dir in dirs:
                self.assertTrue(os.path.isdir(new_dir))

    def test_split_assignment_file(self):
        with tempfile.TemporaryDirectory() as temdir:
            dirs = make_kmer_directories(temdir, "ACGT", 6, complement=False)
            split_assignment_file(self.assignment_path, dirs, "ACGT", 6, 4, min_prob=0.0)
            self.assertTrue(
                os.path.exists(os.path.join(os.path.join(temdir, "GCCTTA"), os.path.basename(self.assignment_path))))
        with tempfile.TemporaryDirectory() as temdir:
            dirs = make_kmer_directories(temdir, "ACGT", 6, complement=False)
            split_assignment_file(self.alignments_path, dirs, "ACGT", 6, 4, min_prob=0.0, alignment=True)
            self.assertTrue(
                os.path.exists(os.path.join(os.path.join(temdir, "TGAAAA"), os.path.basename(self.alignments_path))))

    def test_multiprocess_split_assignment_file(self):
        with tempfile.TemporaryDirectory() as temdir:
            dirs = make_kmer_directories(temdir, "ACGT", 6, complement=False)
            multiprocess_split_sa_tsv_file([self.assignment_path], dirs, "ACGT", 6, min_prob=0.0, worker_count=1)
            self.assertTrue(
                os.path.exists(os.path.join(os.path.join(temdir, "GCCTTA"), os.path.basename(self.assignment_path))))

        with tempfile.TemporaryDirectory() as temdir:
            dirs = make_kmer_directories(temdir, "ACGT", 6, complement=True)
            multiprocess_split_sa_tsv_file([self.alignments_path], dirs, "ACGT", 6, min_prob=0.0,
                                           worker_count=1, alignment=True)
            self.assertTrue(
                os.path.exists(os.path.join(os.path.join(temdir, "TGAAAA"), os.path.basename(self.alignments_path))))

    def test_get_top_kmers_from_directory(self):
        with tempfile.TemporaryDirectory() as temdir:
            dirs = make_kmer_directories(temdir, "ACGT", 6, complement=False)
            split_assignment_file(self.assignment_path, dirs, "ACGT", 6, 4, min_prob=0.0)
            get_top_kmers_from_directory(dirs[2428], temdir, 10, random=False)
            data = []
            with open(os.path.join(temdir, "GCCTTA.tsv"), "r") as fh:
                line = fh.readline()
                data.append(line.split()[3])

            self.assertSequenceEqual(sorted(data)[::-1], data)

    def test_multiprocess_get_top_kmers_from_directory(self):
        with tempfile.TemporaryDirectory() as temdir:
            dirs = make_kmer_directories(temdir, "ACGT", 6, complement=False)
            split_assignment_file(self.assignment_path, dirs, "ACGT", 6, 4, min_prob=0.0)
            multiprocess_get_top_kmers_from_directory(dirs, temdir, 10, random=False)
            data = []
            with open(os.path.join(temdir, "GCCTTA.tsv"), "r") as fh:
                line = fh.readline()
                data.append(line.split()[3])

            self.assertSequenceEqual(sorted(data)[::-1], data)

    def test_generate_top_n_kmers_from_sa_output(self):
        with tempfile.TemporaryDirectory() as temdir:
            output_file = os.path.join(temdir, "built_alignment.tsv")
            generate_top_n_kmers_from_sa_output([self.alignments_path], temdir, output_file, 10,
                                                kmer_len=6, min_prob=0.8, worker_count=1, random=False)

    def test_generate_buildAlignments4(self):
        kmers = get_kmers(6, alphabet="ATGC")
        data_files = [self.alignments_path]
        data, time1 = time_it(multiprocess_make_kmer_assignment_tables,
                              data_files, kmers,
                              {"t", "c"}, 0.0, True, True, 10, 2)

        with tempfile.TemporaryDirectory() as temdir:
            output_file = os.path.join(temdir, "built_alignment.tsv")
            data2, time2 = time_it(generate_top_n_kmers_from_sa_output,
                                   data_files, temdir, output_file, 10, "ACGT", 6, 0.0, 8, False, True, False, True)

            # get kmers associated with each sample
            num_lines = len(list(open(output_file)))
        print(time2, time1)
        self.assertEqual(len(data.index), num_lines)
        self.assertLess(time2, time1)


if __name__ == '__main__':
    unittest.main()
