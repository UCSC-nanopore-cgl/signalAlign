
from __future__ import print_function
import sys
import os
import subprocess
from contextlib import closing
import pysam
import numpy as np

from signalalign import _parseCigar
from Bio import SeqIO


class GuideAlignment(object):
    def __init__(self, cigar=False, strand=False, reference_name=False):
        self.cigar = cigar
        self.strand = strand
        self.reference_name = reference_name

    def validate(self, reference_names=None):
        if self.cigar == "" or self.cigar is False or self.cigar is None:
            print("[GuideAlignment]Illegal cigar %s" % self.cigar)
            return False
        if self.strand not in ("+", "-"):
            print("[GuideAlignment]Illegal strand %s" % self.strand)
            return False
        if reference_names is not None:
            if self.reference_name not in reference_names:
                print("[GuideAlignment]Reference %s not in reference names %s: " % (self.reference_name, reference_names))
                return False
        return True


class TargetRegions(object):
    def __init__(self, tsv, already_sorted=False):
        assert(os.stat(tsv).st_size != 0), "Empty regions file"

        self.region_array = np.loadtxt(tsv, usecols=(0, 1), dtype=np.int32)

        if len(self.region_array.shape) == 1:
            a = np.empty([1, 2], dtype=np.int32)
            a[0] = self.region_array
            self.region_array = a

        if not already_sorted:
            self.region_array = np.sort(self.region_array, axis=1)

    def check_aligned_region(self, left, right):
        if right < left:
            left, right = right, left
        for region in self.region_array:
            if (region[0] >= left) and (region[1] <= right):
                return True
            else:
                continue
        return False


class Bwa(object):
    """Wrapper of BWA aligner, requires bwa to be in path.
    Citation:
        Program: bwa (alignment via Burrows-Wheeler transformation)
        Contact: Heng Li <lh3@sanger.ac.uk>
    """
    def __init__(self, target):
        self.target = target
        self.db_handle = ''

    def build_index(self, destination, output=None, log=None):
        self.db_handle = os.path.join(destination, self.target)
        # is this a directory and are all bwa files present? we can return early
        if False not in set(map(os.path.isfile, ["{}{}".format(self.db_handle, suffix) for suffix in self.suffixes()])):
            return self.db_handle

        if log:
            print("[{}] creating BWA index for {}".format(log, self.db_handle))

        cmd = "bwa index -p {0} {1}".format(self.db_handle, self.target)
        if output is None:
            output = open(os.devnull, 'w')
        else:
            output = open(output, 'w')
        try:
            subprocess.check_call(cmd.split(), stdout=output, stderr=output)
            return self.db_handle
        except subprocess.CalledProcessError:
            return None
        finally:
            output.close()

    @staticmethod
    def suffixes():
        return [".amb", ".ann", ".bwt", ".pac", ".sa"]

    @staticmethod
    def align(reference_fasta, query, output_sam_path, outerr=None):
        """Generate alignment to a reference sequence.

        :param reference_fasta: path reference fasta
        :param query: fasta file to align to reference
        :param output_sam_path: path to place sam file
        :param outerr: another option to write stderr from subprocess command to get alignment
        """
        for suff in Bwa.suffixes():
            assert os.path.exists(reference_fasta + suff),\
                "[Bwa:align] Didn't find index files {}".format(reference_fasta + suff)
        assert os.path.exists(query), "[Bwa::align] Didn't find query file {}".format(query)
        cmd = "bwa mem -x ont2d {idx} {query}".format(idx=reference_fasta, query=query)
        if outerr is None:
            outerr = open(os.devnull, 'w')
        else:
            outerr = open(outerr, 'w')
        try:
            with open(output_sam_path, 'w') as fH:
                # print(bytes.decode(subprocess.check_output(cmd.split(), stderr=outerr)))
                fH.write(bytes.decode(subprocess.check_output(cmd.split(), stderr=outerr)))
            outerr.close()
            return True
        except subprocess.CalledProcessError:
            outerr.close()
            return False


def buildBwaIndex(reference, dest, output=None, log=None):
    """Create BWA index for fasta file

    :param reference: path to fasta reference file
    :param dest: directory to place bwt file
    :param output: file handle to write stdout and stderr. If None output is not captured
    :param log: notation for logging
    """
    assert os.path.isfile(reference), "Reference sequence does not exist: {}".format(reference)
    # check to see if bwt already exists in original ref dir
    if False in set(map(os.path.isfile, ["{}{}".format(reference, suffix) for suffix in Bwa.suffixes()])):
        bwa = Bwa(reference)
        reference = bwa.build_index(dest, output=output, log=log)
    return reference


def getGuideAlignmentFromAlignmentFile(alignment_location, read_name=None, target_regions=None):
    # get reads from alignment file (sam or bam)
    aligned_segment, n_aligned_segments, reference_name = get_aligned_segment_from_alignment_file(alignment_location, read_name)
    assert aligned_segment is not None, "[generateGuideAlignment] Read {} has no passing alignment".format(read_name)
    query_name = aligned_segment.qname
    flag = aligned_segment.flag
    reference_pos = aligned_segment.pos + 1  # pysam gives the 0-based leftmost start
    sam_cigar = aligned_segment.cigarstring

    # couldn't find anything
    if n_aligned_segments == 0:
        print("[generateGuideAlignment] Found no aligned segments" +
              ("" if read_name is None else " for read {}".format(read_name)))
        return None

    if n_aligned_segments > 1:
        print("[generateGuideAlignment] WARNING more than 1 mapping, taking the first one heuristically")

    # get strand
    strand = ""
    assert(flag is not None), "[generateGuideAlignment] ERROR flag is None"

    if int(flag) == 16:
        strand = "-"
    if int(flag) == 0:
        strand = "+"
    elif int(flag) != 0 and int(flag) != 16:
        print("[generateGuideAlignment] ERROR unexpected alignment flag {flag}, not continuing with signal alignment"
              " for {query}".format(flag=flag, query=query_name), file=sys.stderr)
        return None

    # get cigar info
    try:
        query_start, query_end, reference_start, reference_end, cigar_string = _parseCigar(sam_cigar, reference_pos,
                                                                                           forward=strand == "+")
    except AssertionError as e:
        print("[generateGuideAlignment] ERROR %s" % e)
        return None

    # account for strand
    if strand == '-':
        temp = reference_start
        reference_start = reference_end
        reference_end = temp

    # santity
    assert(reference_name is not None), "[generateGuideAlignment] ERROR reference_name is None"
    assert(query_name is not None), "[generateGuideAlignment] ERROR query_name is None"

    completeCigarString = "cigar: %s %i %i + %s %i %i %s 1 %s" % (
        query_name, query_start, query_end, reference_name, reference_start, reference_end, strand, cigar_string)

    if target_regions is not None:
        target_regions = target_regions if type(target_regions) == TargetRegions else TargetRegions(target_regions)
        keep = target_regions.check_aligned_region(reference_start, reference_end)
        if keep is False:
            print("[generateGuideAlignment] Read does not map witin the target regions, passing "
                  "on signal-level alignment", file=sys.stderr)
            return None
        else:
            pass

    return GuideAlignment(completeCigarString, strand, reference_name)


def getInfoFromCigarFile(cigar_file):
    with open(cigar_file, 'r') as cig:
        for line in cig:
            if not line.startswith("cigar:"): continue
            parts = line.split()
            assert len(parts) >= 11, "[generateGuideAlignment] Malformed cigar file {}".format(cigar_file)
            strand = parts[8]
            assert strand == '+' or strand == "-", "[generateGuideAlignment] Unexpected strand '{}' in cigar line {}".format(strand, line)
            reference_name = parts[5]
            return strand, reference_name
    assert False, "[generateGuideAlignment] No cigar line found in {}".format(cigar_file)


def generateGuideAlignment(reference_fasta, query, temp_sam_path, target_regions=None):
    # type: (string, string, string, TargetRegions) -> GuideAlignment
    """Aligns the read sequnece with BWA to get the guide alignment,
    returns the CIGAR (in exonerate format), the strand (plus or minus) and the
    contig mapped to if the read aligned. Returns (False, False, False) if there
    is a problem with any of the steps or if the read maps to a region not included
    within TargetRegions
    """

    # align with bwa
    ok = Bwa.align(reference_fasta=reference_fasta, query=query, output_sam_path=temp_sam_path)
    if not ok:  # Bwa alignment fail
        return None

    return getGuideAlignmentFromAlignmentFile(temp_sam_path, target_regions=target_regions)


def get_aligned_segment_from_alignment_file(alignment_location, read_name):
    """Find an AlignedSegment from an alignment file via the pysam.AlignmentFile interface

    :param alignment_location: sam or bam file
    :param read_name: name of the read you want the alignment information
    :return: the pysam AlignedSegment object
    """
    # get reads from alignment file (sam or bam)
    n_aligned_segments = 0
    correct_segment = None
    reference_name = None
    with closing(pysam.AlignmentFile(alignment_location, 'rb' if alignment_location.endswith("bam") else 'r')) as aln:
        for aligned_segment in aln.fetch():
            if aligned_segment.is_secondary or aligned_segment.is_unmapped \
                    or aligned_segment.is_supplementary or aligned_segment.has_tag("SA"):
                continue
            if read_name is not None and aligned_segment.qname != read_name and read_name not in aligned_segment.qname:
                continue
            n_aligned_segments += 1
            correct_segment = aligned_segment
        if correct_segment is not None:
            reference_name = aln.getrname(correct_segment.rname)

    return correct_segment, n_aligned_segments, reference_name


def check_number_of_passing_reads(sam, fastq):
    """Calculate fraction of passing reads given original fastq and sam file

    note: assumes phred-33 encoding

    :param sam: sam/bam file containing alignment information
    :param fastq: compiled fastq from nanopore reads
    """
    assert os.path.exists(fastq), "Fastq file does not exist: {}".format(fastq)
    assert os.path.exists(sam), "Sam file does not exist: {}".format(sam)

    total = 0
    passed = 0
    not_primary = 0
    av_phred_score_count = 0

    for record in SeqIO.parse(fastq, "fastq"):
        av_phred_score = np.mean(record.letter_annotations["phred_quality"])
        correct_segment, n_aligned_segments, reference_name = get_aligned_segment_from_alignment_file(sam, record.id)
        if correct_segment is None and n_aligned_segments != 1:
            not_primary += 1
            if av_phred_score < 7:
                av_phred_score_count += 1
        elif av_phred_score < 7:
            av_phred_score_count += 1
        else:
            passed += 1

        total += 1

    print("total = {}\npassed = {}\nno_suitable_alignment = {}\nnum av_phred_score < 7 = {}".format(total,
                                                                                                passed,
                                                                                                not_primary,
                                                                                                av_phred_score_count))

    return total, passed
