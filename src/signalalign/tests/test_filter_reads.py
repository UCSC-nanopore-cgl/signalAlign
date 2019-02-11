#!/usr/bin/env python
"""
    Place unit tests for filter_reads.py
"""
########################################################################
# File: test_filter_reads.py
#  executable: test_filter_reads.py
# Purpose: test filter_reads.py
#
# Author: Andrew Bailey
# History: 11/1/2018 Created
########################################################################

import unittest
import tempfile
import shutil
import pysam
from signalalign.filter_reads import *


class FilterReadsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super(FilterReadsTest, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-4])

        cls.tmp_directory = tempfile.mkdtemp()
        cls.test_dir = os.path.join(cls.tmp_directory, "test")

        cls.dna_dir = os.path.join(cls.HOME, "tests/minion_test_reads/1D/")
        # copy file to tmp directory
        shutil.copytree(cls.dna_dir, cls.test_dir)
        cls.readdb = os.path.join(cls.HOME, "tests/minion_test_reads/oneD.fastq.index.readdb")
        cls.bam = os.path.join(cls.HOME, "tests/minion_test_reads/oneD.bam")

    def test_write_readdb(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            tmp_file = os.path.join(temp_dir, "test.readdb")
            write_this = [[name, path] for name, path in parse_readdb(self.readdb, [self.test_dir])]
            write_readdb(write_this, tmp_file)
            read_this = [[name, path] for name, path in parse_readdb(tmp_file, [self.test_dir])]
            self.assertSequenceEqual(read_this, write_this)

    def test_parse_readdb(self):
        for name, path in parse_readdb(self.readdb, [self.test_dir]):
            self.assertTrue(name.endswith("template"))
            self.assertTrue(os.path.exists(os.path.join(self.test_dir, path)))

    def test_filter_reads(self):
        first_dir = {x: y for x, y in filter_reads(self.bam, self.readdb, [self.test_dir], quality_threshold=7)}
        self.assertEqual(len(first_dir), 3)
        first_dir = {x: y for x, y in filter_reads(self.bam, self.readdb, [self.test_dir], quality_threshold=12)}
        self.assertEqual(len(first_dir), 0)
        first_dir = {x: y for x, y in filter_reads(self.bam, self.readdb, [os.path.split(self.test_dir)[0]],
                                                   quality_threshold=7, recursive=True)}
        self.assertEqual(len(first_dir), 3)

    def test_filter_reads_to_string_wrapper(self):
        a = filter_reads(self.bam, self.readdb, [self.test_dir], quality_threshold=7)
        for x, y in a:
            self.assertTrue(type(x) is str)
            self.assertTrue(type(y) is pysam.AlignedSegment)

        a = filter_reads(self.bam, self.readdb, [self.test_dir], quality_threshold=7)
        b = filter_reads_to_string_wrapper(a)
        for x, y in b:
            self.assertTrue(type(x) is str)
            self.assertTrue(type(y) is str)

    def test_multiprocess_move_and_filter_reads(self):
        with tempfile.TemporaryDirectory() as tempdir:
            multiprocess_move_and_filter_reads(self.tmp_directory, tempdir, self.bam, self.readdb, trim=False,
                                               quality_threshold=False, worker_count=1, debug=True)
            self.assertEqual(3, len(os.listdir(os.path.join(tempdir, "test"))))
            shutil.rmtree(self.test_dir)
            shutil.copytree(self.dna_dir, self.test_dir)

        with tempfile.TemporaryDirectory() as tempdir:
            multiprocess_move_and_filter_reads(self.tmp_directory, tempdir, self.bam, self.readdb, trim=10,
                                               quality_threshold=False, worker_count=1, debug=True)
            self.assertEqual(1, len(os.listdir(os.path.join(tempdir, "test"))))
            shutil.rmtree(self.test_dir)
            shutil.copytree(self.dna_dir, self.test_dir)
        with tempfile.TemporaryDirectory() as tempdir:
            multiprocess_move_and_filter_reads(self.tmp_directory, tempdir, self.bam, self.readdb, trim=False,
                                               quality_threshold=False, worker_count=1, debug=False)
            self.assertEqual(3, len(os.listdir(os.path.join(tempdir, "test"))))
            shutil.rmtree(self.test_dir)
            shutil.copytree(self.dna_dir, self.test_dir)

    def test_multiprocess_filter_reads(self):
        reads = multiprocess_filter_reads(self.tmp_directory, self.bam, self.readdb, trim=False,
                                          quality_threshold=False, worker_count=1, debug=True)
        self.assertEqual(len(reads), len(list_dir(self.test_dir, ext="fast5")))
        counter = 0
        for read, sam_str in multiprocess_filter_reads(self.tmp_directory, self.bam, self.readdb, trim=False,
                                                       quality_threshold=False, worker_count=1, debug=False):
            self.assertTrue(os.path.exists(read))
            counter += 1
        self.assertEqual(counter, len(list_dir(self.test_dir, ext="fast5")))

    def test_find_fast5s_from_ids_readdb(self):
        read_ids = ["5cc86bac-79fd-4897-8631-8f1c55954a45"]
        path = [y for x, y in find_fast5s_from_ids_readdb(self.readdb, read_ids, [self.dna_dir])]
        self.assertEqual(os.path.basename(path[0]), "LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch92_read1108_strand.fast5")

    def test_copy_files_from_readdb(self):
        with tempfile.TemporaryDirectory() as tempdir:
            n_files = copy_files_from_readdb(self.readdb, [self.dna_dir], tempdir, recursive=False)
            self.assertEqual(n_files, 3)
            self.assertEqual(len(list_dir(tempdir)), 3)

    @classmethod
    def tearDownClass(cls):
        """Remove test fast5 file"""
        shutil.rmtree(cls.tmp_directory)


if __name__ == '__main__':
    unittest.main()
