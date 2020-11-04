#!/usr/bin/env python3

import glob
import os
import shutil
import tempfile
import unittest
from contextlib import closing
from subprocess import call, PIPE

import numpy as np
import pandas as pd
from py3helpers.utils import captured_output
from signalalign.nanoporeRead import NanoporeRead
from signalalign.utils.sequenceTools import getFastaDictionary

SIGNALALIGN_ROOT = '/'.join(os.path.abspath(__file__).split("/")[:-4])
ZYMO_C_READS = os.path.join(SIGNALALIGN_ROOT, "tests/minion_test_reads/C/")
ZYMO_REFERENCE = os.path.join(SIGNALALIGN_ROOT, "tests/test_sequences/zymo_sequence.fasta")
BIN_PATH = os.path.join(SIGNALALIGN_ROOT, "bin")

scripts_PATH = os.path.join(SIGNALALIGN_ROOT, "src/signalalign/scripts/")


def parse_alignment_full(alignment_file):
    data = pd.read_csv(alignment_file,
                       usecols=(1, 2, 4, 5, 9, 12, 13),
                       sep="\t",
                       dtype={'ref_pos': np.int64,
                              'ref_kmer': np.str,
                              'strand': np.str,
                              'event_index': np.int64,
                              'kmer': np.str,
                              'posterior_prob': np.float64,
                              'event_mean': np.float64},
                       header=None,
                       names=['ref_pos', 'ref_kmer', 'strand', 'event_index', 'kmer', 'posterior_prob', 'event_mean'])
    return data


class LibTest(unittest.TestCase):
    def test_signalAlign_library(self):
        current_wd = os.getcwd()
        command = os.path.join(BIN_PATH, "signalAlignLibTests")
        os.chdir(BIN_PATH)
        result = call(command, shell=True, bufsize=-1, stdout=PIPE, stderr=PIPE)
        self.assertTrue(result == 0, "signalAlign Library Tests Fail")
        os.chdir(current_wd)


class SignalAlignLibTests(unittest.TestCase):
    def setUp(self):
        self.work_dir = "./signalAlign_pylibTest/"
        os.makedirs(self.work_dir)

    def tearDown(self):
        shutil.rmtree(self.work_dir)


class SignalAlignAlignmentTest(unittest.TestCase):
    def setUp(self):
        self.current_wd = os.getcwd()
        os.chdir(BIN_PATH)
        if os.path.exists("./signalAlign_unittest/"):
            shutil.rmtree("./signalAlign_unittest/")
        os.makedirs("./signalAlign_unittest/")

    def tearDown(self):
        # if os.path.exists("./signalAlign_unittest/"):
        #     shutil.rmtree("./signalAlign_unittest/")
        os.chdir(self.current_wd)

    def check_alignments(self, reads, reference, kmer_length, contig_name, extra_args=None, rna=False):

        def get_kmer(start):
            kmer = reference_sequence[start:start + kmer_length]
            if type(kmer) is str:
                return kmer
            else:
                return bytes.decode(kmer)

        input_fast5s = glob.glob(os.path.join(reads, "*.fast5"))
        assert len(input_fast5s) > 0, "Didn't find test MinION reads"
        assert os.path.isfile(reference), "Didn't find reference sequence"

        # it's this or rewrite all the relative locations of the files
        os.chdir(BIN_PATH)

        # prep command
        run_signal_align = os.path.join(scripts_PATH, "runSignalAlign.py")
        # removed: --debug
        alignment_command = "python {runsignalalign} run2 -d={reads} --bwa_reference={ref} " \
                            "-smt=threeState -o={testDir} " \
                            "--keep_tmp_folder ".format(runsignalalign=run_signal_align,
                                                        reads=reads,
                                                        ref=reference,
                                                        testDir="./signalAlign_unittest/")
        if extra_args is not None:
            alignment_command += extra_args

        # run signalAlign
        result = call(alignment_command, shell=True, bufsize=-1, stdout=PIPE, stderr=PIPE)

        self.assertTrue(result == 0, "Error running signalAlign. Command was {}".format(alignment_command))

        # get alignments
        test_alignments = glob.glob("./signalAlign_unittest/tempFiles_alignment/*.tsv")
        self.assertTrue(len(test_alignments) == len(input_fast5s),
                        "Didn't make all alignments got {got} should be {should}".format(
                            got=len(test_alignments), should=len(input_fast5s)))

        # prep for verification
        reference_sequence = getFastaDictionary(reference)[contig_name]
        alignment2events = dict()
        for fast5 in input_fast5s:
            with closing(NanoporeRead(fast5, initialize=True)) as read:
                event_count = len(read.get_template_events())
                read_id = read.read_label
                self.assertTrue(event_count > 0, "Got no events for fast5 {} with read_id {}".format(fast5, read_id))
                for alignment in test_alignments:
                    if os.path.basename(alignment).startswith(read_id):
                        self.assertTrue(alignment not in alignment2events,
                                        "Fast5 {} matched read_id {} with multiple output alignments".format(
                                            fast5, read_id))
                        alignment2events[alignment] = event_count

        for alignment in test_alignments:
            # expected = parse_alignment_full(os.path.join(true_alignments, alignment_file))
            obs = parse_alignment_full(alignment)
            for row in obs.itertuples():
                ref_pos = row[1]
                obs_kmer = row[2]
                strand = row[3]
                exp_kmer = get_kmer(ref_pos)
                if rna:
                    exp_kmer = exp_kmer[::-1]
                self.assertEqual(obs_kmer, exp_kmer, msg="kmer at index {idx} on strand {strand} is {obs} "
                                                         "should be {exp}, file {f}".format(idx=ref_pos,
                                                                                            strand=strand,
                                                                                            obs=obs_kmer,
                                                                                            exp=exp_kmer,
                                                                                            f=alignment))
            signal_align_event_count = len(obs)
            initial_event_count = alignment2events[alignment]
            self.assertTrue(signal_align_event_count >= initial_event_count,
                            "SignalAlign produced {} events, less than initial count {}".format(
                                signal_align_event_count, initial_event_count))
            # this is a magic number
            self.assertTrue(signal_align_event_count <= initial_event_count * 3,
                            "SignalAlign produced {} events, more than 3x the initial count {}".format(
                                signal_align_event_count, initial_event_count))

    def test_pUC_r9_reads_5mer(self):
        with captured_output() as (_, _):

            puc_reads = os.path.join(SIGNALALIGN_ROOT, "tests/minion_test_reads/pUC/")
            puc_reference = os.path.join(SIGNALALIGN_ROOT, "tests/test_sequences/pUC19_SspI.fa")
            self.check_alignments(reads=puc_reads,
                                  reference=puc_reference,
                                  kmer_length=5,
                                  contig_name="pUC19",
                                  extra_args="-T=../models/testModelR9_5mer_acegot_template.model "
                                             "-C=../models/testModelR9_5mer_acegot_complement.model "
                                             "--2d ")

    def test_pUC_r9_reads_6mer(self):
        with captured_output() as (_, _):

            puc_reads = os.path.join(SIGNALALIGN_ROOT, "tests/minion_test_reads/pUC/")
            puc_reference = os.path.join(SIGNALALIGN_ROOT, "tests/test_sequences/pUC19_SspI.fa")
            self.check_alignments(reads=puc_reads,
                                  reference=puc_reference,
                                  kmer_length=6,
                                  contig_name="pUC19",
                                  extra_args="--2d ")

    def test_Ecoli1D_reads_5mer(self):
        with captured_output() as (_, _):

            reads = os.path.join(SIGNALALIGN_ROOT, "tests/minion_test_reads/1D/")
            reference = os.path.join(SIGNALALIGN_ROOT, "tests/test_sequences/E.coli_K12.fasta")
            self.check_alignments(reads=reads,
                                  reference=reference,
                                  kmer_length=5,
                                  contig_name="gi_ecoli",
                                  extra_args="-T=../models/testModelR9p4_5mer_acegt_template.model ")

    def test_RNA_edge_alignments_reads_5mer(self):
        with captured_output() as (_, _):

            edge_case_reads = os.path.join(SIGNALALIGN_ROOT, "tests/minion_test_reads/RNA_edge_cases/")
            with tempfile.TemporaryDirectory() as tempdir:
                new_dir = os.path.join(tempdir, "new_dir")
                shutil.copytree(edge_case_reads, new_dir)
                edge_case_reference = os.path.join(SIGNALALIGN_ROOT, "tests/test_sequences/fake_rna_ref.fa")
                rna_alignments = os.path.join(SIGNALALIGN_ROOT, "tests/minion_test_reads/RNA_edge_cases/rna_reads.sam")
                self.check_alignments(reads=new_dir,
                                      reference=edge_case_reference,
                                      kmer_length=5,
                                      contig_name="rna_fake",
                                      extra_args="-T=../models/testModelR9p4_5mer_acgt_RNA.model "
                                                 "--alignment_file {}".format(rna_alignments), rna=True)

    def test_signal_files_without_events(self):
        """Test if signalAlign can handle signal files without event information"""
        with captured_output() as (_, _):

            signal_file_reads = os.path.join(SIGNALALIGN_ROOT, "tests/minion_test_reads/no_event_data_1D_ecoli")

            with tempfile.TemporaryDirectory() as tempdir:
                new_dir = os.path.join(tempdir, "new_dir")
                shutil.copytree(signal_file_reads, new_dir)
                ecoli_reference = os.path.join(SIGNALALIGN_ROOT, "tests/test_sequences/E.coli_K12.fasta")
                signal_file_guide_alignment = os.path.join(SIGNALALIGN_ROOT, "tests/minion_test_reads/oneD.bam")
                self.check_alignments(reads=new_dir,
                                      reference=ecoli_reference,
                                      kmer_length=5,
                                      contig_name="gi_ecoli",
                                      extra_args="-T=../models/testModelR9p4_5mer_acegt_template.model "
                                                 "--alignment_file {}".format(signal_file_guide_alignment))


if __name__ == '__main__':
    unittest.main()
