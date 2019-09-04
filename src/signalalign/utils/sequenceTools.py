#!/usr/bin/env python3
"""Functions to create a positions file."""
########################################################################
# File: makePositionsFiles.py
#  executable: makePositionsFiles.py
#
# Author: Andrew Bailey
# History: 5/21/18 Created
########################################################################

import os
import sys
import string
import array
import subprocess
import numpy as np
import pandas as pd
import pysam
from contextlib import closing
from collections import Counter
from signalalign.utils.parsers import read_fasta
from py3helpers.utils import find_substring_indices, all_string_permutations


AMBIG_BASES = {
    "AG": "R",
    "CT": "Y",
    "CG": "S",
    "AT": "W",
    "GT": "K",
    "AC": "M",
    "CGT": "B",
    "AGT": "D",
    "ACT": "H",
    "ACG": "V",
    "ACGT": "X",
    "CEO": "L",
    "CE": "P",
    "AI": "Q",
    "AF": "f",
    "ACEGOT": "U",
    "JT": "Z",
    "A": "A",
    "T": "T",
    "C": "C",
    "G": "G",
    "E": "E",
    "O": "O",
    "F": "F",
    "J": "J",
    "I": "I",
    "R": "R",
    "Y": "Y",
    "S": "S",
    "W": "W",
    "K": "K",
    "M": "M",
    "B": "B",
    "D": "D",
    "H": "H",
    "V": "V",
    "X": "X",
    "L": "L",
    "P": "P",
    "Q": "Q",
    "f": "f",
    "U": "U",
    "Z": "Z",
    "p": "p",
    "b": "b",
    "d": "d",
    "e": "e",
    "h": "h",
    "i": "i",
    "j": "j",
    "k": "k",
    "l": "l",
    "m": "m",
    "n": "n",
    "o": "o",
    "Tp": "j",
    "Gb": "k",
    "Gd": "l",
    "Ce": "m",
    "Th": "n",
    "Ai": "o"
}


def find_gatc_motifs(sequence):
    """Generate index of 'A' within the 'GATC' motifs in a nucleotide sequence

    :param sequence: since GATC motif is in DNA, expecting a DNA nucleotide sequence
    :return: generator yielding index of 'A' within the 'GATC'
    """
    return find_substring_indices(sequence.upper(), "GATC", offset=1)


def find_different_char_index(start_string, edit_string):
    """Compares standard and modified string and identifies the index of the modified character
        ex: find_char_difference_index("CCAGG","CFAGG") = 1

    :param start_string: starting string
    :param edit_string: string with a single character edit from start_string
    :return: index of string difference
    """
    assert len(start_string) == len(edit_string), ""
    pos = [i for i in range(len(start_string)) if start_string[i] != edit_string[i]]
    assert len(pos) == 1, "Only one character difference allowed. " \
                          "start_string={}, edit_string={}".format(start_string, edit_string)
    return pos[0]


def find_modification_index_and_character(canonical_motif, replacement_motif):
    """Compares canonical and modified motif and identifies the
            index of the modified nucleotide, canonical character, and replacement character.

    note. everything is converted to uppercase

    ex. find_modification_index_and_character("ATGC", "ETGC") = 0, "A", "E"

    :param canonical_motif: canonical nucleotide bases
    :param replacement_motif: replacement motif
    :return: mod_index, old_char, new_char
    """
    canonical_motif = canonical_motif.upper()
    replacement_motif = replacement_motif.upper()
    assert canonical_motif != replacement_motif, "Canonical motif cannot be the same as replacement motif"
    assert set(canonical_motif) <= set("ATGC"), "Canonical motif must only have canonical nucleotides"
    pos = find_different_char_index(canonical_motif, replacement_motif)
    old_char = canonical_motif[pos]
    new_char = replacement_motif[pos]
    return pos, old_char, new_char


def make_positions_file(reference, output_path, motifs, overlap=False):
    """Creates a tsv file with the following format ("contig", "position", "strand", "change_from", "change_to").
    Given a reference sequence and sets of sequence motif changes we report the location of each change.

    NOTE: the motifs cannot create new characters on the opposite strand!

    :param reference: path to reference sequence
    :param output_path: output path of positions file
    :param motifs: list of lists of find replace motifs ex: [("CCAGG","CFAGG"), ("CCTGG","CFTGG")]
    :param overlap: if the motif can overlap with its self, find index of overlap if set to true
    """
    rev_motifs = []
    for motif in motifs:
        rev_motifs.append([x[::-1] for x in motif])

    with open(output_path, "w") as outfile:
        for header, comment, sequence in read_fasta(reference):
            fwd_seq = sequence
            bwd_seq = reverse_complement(fwd_seq, reverse=False, complement=True).upper()
            for index, old_char, substitution_char in find_motifs_sequence_positions(fwd_seq, motifs, overlap=overlap):
                outfile.write(header + "\t" + np.str(index) + "\t" + "+" + "\t"
                              + old_char + "\t" + substitution_char + "\n")

            for index, old_char, substitution_char in find_motifs_sequence_positions(bwd_seq, rev_motifs, overlap=overlap):
                outfile.write(header + "\t" + np.str(index) + "\t" + "-" + "\t"
                              + old_char + "\t" + substitution_char + "\n")

    return output_path


def replace_motifs_sequence_positions(sequence, motifs, overlap=False):
    """Edit nucleotide sequence using find and replace motifs

    note: we convert sequence to uppercase

    :param sequence: nucleotide sequence
    :param motifs: list of motif's which need to be replaced: eg [[find, replace]], [["CCAGG", "CEAGG"]]
    :param overlap: boolean option to look for motif overlaps
    """
    new_sequence = list(sequence)
    for index, old_char, substitution_char in find_motifs_sequence_positions(sequence, motifs, overlap=overlap):
        new_sequence[index] = substitution_char
    subst_sequence = ''.join(new_sequence).upper()
    return subst_sequence


def find_motifs_sequence_positions(sequence, motifs, overlap=False):
    """Find locations of edited nucleotide nucleotide sequence using find and replace motifs

    note: we convert sequence to uppercase

    :param sequence: nucleotide sequence
    :param motifs: list of motif's which need to be replaced: eg [[find, replace]], [["CCAGG", "CEAGG"]]
    :param overlap: boolean option to look for motif overlaps
    """
    already_repaced_indexes = set()
    # gather motifs
    for motif_pair in motifs:
        assert len(motif_pair) == 2 and type(
            motif_pair) is list, "Motifs must be structured as list of lists, even for one motif find and replace"
        # find edit character and offset
        offset, old_char, substitution_char = find_modification_index_and_character(motif_pair[0], motif_pair[1])
        for index in find_substring_indices(sequence.upper(), motif_pair[0].upper(), offset=offset, overlap=overlap):
            # make sure that there is no overlapping assignments of characters
            assert index not in already_repaced_indexes, "Motifs has two different edits to a single nucleotide " \
                                                         "location. Check motifs {}".format(motifs)
            already_repaced_indexes.add(index)

            yield index, old_char, substitution_char


def replace_periodic_sequence_positions(sequence, step_size, offset, substitution_char):
    """Edit nucleotide sequence using by replacing every 'step_size' nucleotides with an offset

    note: we convert sequence to uppercase

    eg: replace_periodic_sequence_positions("ATGCATGC", 3, 1, "F") = "AFGCFTGF"

    :param sequence: nucleotide sequence
    :param step_size: every 'step_size' locations the offset position is changed
    :param offset: the
    :param substitution_char: replacement character
    """
    assert offset < step_size, "Offset has to be less than step size"
    sequence = list(sequence)
    for i in range(offset, len(sequence), step_size):
        sequence[i] = substitution_char
    subst_sequence = ''.join(sequence).upper()

    return subst_sequence


def replace_periodic_reference_positions(reference_location, sub_fasta_path, step, offset, substitution_char='X'):
    """Edit and write a reference sequence to a specified path by replacing periodic characters

    note: if sub_fasta_path exists it will return the path without creating a new file

    :param reference_location: input reference
    :param sub_fasta_path: location of edited reference
    :param step: size of gap between substitution characters
    :param offset: offset of when to start creating substiutions
    :param substitution_char: character to replace original character
    """
    if os.path.isfile(sub_fasta_path):
        print("[substitute_reference_positions] Substituted reference fasta file exists: {}".format(
            sub_fasta_path))
        return sub_fasta_path
    else:
        print("[substitute_reference_positions] Creating substituted reference fasta file: {}".format(
            sub_fasta_path))
        # write
        with open(sub_fasta_path, 'w') as outfasta:
            for header, comment, sequence in read_fasta(reference_location):
                subst_sequence = replace_periodic_sequence_positions(sequence, step, offset, substitution_char)
                print(
                    ">%s %s\n%s" % (header, "substituted:{},step:{},offset:{}".format(substitution_char, step, offset),
                                    subst_sequence), file=outfasta)

    return sub_fasta_path


def replace_motif_reference_positions(reference_location, sub_fasta_path, motifs, overlap=False):
    """Replace motif  reference sequence to a specific path

    :param reference_location: input reference
    :param sub_fasta_path: location of edited reference
    :param motifs: list of motif's which need to be replaced: eg [[find, replace]], [["CCAGG", "CEAGG"]]
    :param overlap: of overlap is possible, replace with overlap: eg [["AAA", "AAT"]] :  AAAA -> AATT
    """
    if os.path.isfile(sub_fasta_path):
        print("[substitute_reference_positions] Substituted reference fasta file exists: {}".format(
            sub_fasta_path))
        return sub_fasta_path
    else:
        print("[substitute_reference_positions] Creating substituted reference fasta file: {}".format(
            sub_fasta_path))
        # write
        with open(sub_fasta_path, 'w') as outfasta:
            for header, comment, sequence in read_fasta(reference_location):
                subst_sequence = replace_motifs_sequence_positions(sequence, motifs, overlap)
                print(">%s %s\n%s" % (header, "substituted:{}".format(motifs),
                                      subst_sequence), file=outfasta)
    return sub_fasta_path


def samtools_faidx_fasta(fasta_path, log=None):
    """Index fasta using samtools faidx

    note: samtools must be in PATH

    :param fasta_path: path to fasta file
    """
    # produce faidx file
    assert os.path.isfile(fasta_path), "Path to fasta file does not exist"
    index_path = "{}.fai".format(fasta_path)
    if not os.path.exists(index_path):
        if log:
            print("[{}] indexing reference {}".format(log, fasta_path))
        args = ["samtools", "faidx", fasta_path]
        subprocess.check_call(args)
    assert os.path.isfile(index_path), "Error creating FAIDX file for: {}".format(fasta_path)
    return index_path


def count_all_sequence_kmers(seq, k=5, rev_comp=False):
    """Count all the 5'-3' kmers of a nucleotide sequence, rev_comp counts rev_comp seq IN ADDITION to given sequence

    :param seq: nucleotide sequence
    :param k: size of kmer
    :param rev_comp: boolean option to count reverse complement kmers as well
    :return: dictionary of kmers with counts as values
    """
    # loop through kmers
    kmers = Counter()
    for kmer in kmer_iterator(seq, k):
        kmers[kmer] += 1
    if rev_comp:
        # loop through rev_comp kmers
        seq1 = reverse_complement(seq, reverse=True, complement=True)
        for kmer in kmer_iterator(seq1, k):
            kmers[kmer] += 1

    return kmers


def get_sequence_kmers(seq, k=5, rev_comp=False):
    """Get the set of all kmers from a sequence.

    :param seq: nucleotide sequence
    :param k: size of kmer
    :param rev_comp: boolean option to count reverse complement kmers as well
    :return: set of kmers
    """
    return set(count_all_sequence_kmers(seq, k=k, rev_comp=rev_comp).keys())


def get_motif_kmers(motif_pair, k, alphabet="ATGC"):
    """Given a motif pair, create a list of all kmers which contain modification


    """
    assert len(motif_pair) == 2, "Motif pair must be a list of length 2. len(motif_pair) = {}".format(len(motif_pair))
    canonical = motif_pair[0]
    modified = motif_pair[1]
    motif_len = len(canonical)
    # get mod index and chars
    mod_index, old_char, new_char = find_modification_index_and_character(canonical, modified)
    bases_after = motif_len - mod_index - 1

    # get overlaps for front and back of kmer
    front_overlap, back_overlap = get_front_back_kmer_overlap(k, motif_len, mod_index)
    # pre-compute kmers
    kmer_set_dict = dict()
    for i in range(1, max(front_overlap, back_overlap) + 1):
        kmer_set_dict[i] = [x for x in all_string_permutations(alphabet, i)]
    kmer_set_dict[0] = ['']

    motif_kmers = []
    for i in range(k):
        # get prepend kmers and index for front of motif
        if i >= front_overlap:
            front_index = i - front_overlap
            prepend_kmers = ['']
        else:
            prepend_kmers = kmer_set_dict[front_overlap - i]
            front_index = 0
        # get append kmers and index for back of motif
        if i > bases_after:
            append_kmers = kmer_set_dict[i - bases_after]
            back_index = motif_len
        else:
            back_index = mod_index + i + 1
            append_kmers = ['']

        kmer = modified[front_index:back_index]
        motif_kmers.extend(
            [front + kmer + back for front in prepend_kmers for back in append_kmers if front + kmer + back is not ''])

    return set(motif_kmers)


def get_front_back_kmer_overlap(k, motif_len, mod_index):
    """Get the largest number of bases overlap at front and back of motif

    eg: k=3 , motif_len = 2, mod_index = 1
        motif = GE

        X G E X
        _ _ _       max front_overlap = 1
          _ _ _     max back_overlap = 1

    :param k: length of kmer
    :param motif_len: length of motif
    :param mod_index: index position of modification
    :return: largest overlap in the front and back of a generated kmer
    """
    assert k >= 1, "k cannot be less than 1. k: {}".format(k)
    front_overlap = k - mod_index - 1
    back_overlap = k - (motif_len - mod_index)
    return front_overlap, back_overlap


# TODO write these tests ya dig
def getFastaDictionary(fastaFile):
    """Returns a dictionary of the first words of fasta headers to their corresponding
    fasta sequence
    """
    namesAndSequences = [(x[0].split()[0], x[1]) for x in fastaRead(open(fastaFile, 'r'))]
    names = [x[0] for x in namesAndSequences]
    assert len(names) == len(set(names))  # Check all the names are unique
    return dict(namesAndSequences)  # Hash of names to sequences


def fastaRead(fileHandleOrFile):
    """iteratively yields a sequence for each '>' it encounters, ignores '#' lines
    """
    fileHandle = _getFileHandle(fileHandleOrFile)
    line = fileHandle.readline()
    chars_to_remove = "\n "
    valid_chars = {x for x in string.ascii_letters + "-"}
    while line != '':
        if line[0] == '>':
            name = line[1:-1]
            line = fileHandle.readline()
            seq = array.array('b')
            while line != '' and line[0] != '>':
                line = line.translate(str.maketrans('', '', chars_to_remove))
                if len(line) > 0 and line[0] != '#':
                    seq.extend(list(map(ord, line)))
                line = fileHandle.readline()
            try:
                assert all(chr(x) in valid_chars for x in seq)
            except AssertionError:
                bad_chars = {chr(x) for x in seq if chr(x) not in valid_chars}
                raise RuntimeError("Invalid FASTA character(s) see in fasta sequence: {}".format(bad_chars))
            yield name, seq.tobytes()
        else:
            line = fileHandle.readline()
    if isinstance(fileHandleOrFile, "".__class__):
        fileHandle.close()


def _getFileHandle(fileHandleOrFile, mode="r"):
    if isinstance(fileHandleOrFile, "".__class__):
        return open(fileHandleOrFile, mode)
    else:
        return fileHandleOrFile


def fastaWrite(fileHandleOrFile, name, seq, mode="w"):
    """Writes out fasta file
    """
    fileHandle = _getFileHandle(fileHandleOrFile, mode)
    valid_chars = {x for x in string.ascii_letters + "-"}
    try:
        assert any([isinstance(seq, str), isinstance(seq, str)])
    except AssertionError:
        raise RuntimeError("Sequence is not unicode or string")
    try:
        assert all(x in valid_chars for x in seq)
    except AssertionError:
        bad_chars = {x for x in seq if x not in valid_chars}
        raise RuntimeError("Invalid FASTA character(s) see in fasta sequence: {}".format(bad_chars))
    fileHandle.write(">%s\n" % name)
    chunkSize = 100
    for i in range(0, len(seq), chunkSize):
        fileHandle.write("%s\n" % seq[i:i + chunkSize])
    if isinstance(fileHandleOrFile, "".__class__):
        fileHandle.close()


def getFastaDictionary(fastaFile):
    """Returns a dictionary of the first words of fasta headers to their corresponding
    fasta sequence
    """
    namesAndSequences = [(x[0].split()[0], x[1]) for x in fastaRead(open(fastaFile, 'r'))]
    names = [x[0] for x in namesAndSequences]
    assert len(names) == len(set(names))  # Check all the names are unique
    return dict(namesAndSequences)  # Hash of names to sequences


def kmer_iterator(dna, k):
    """Generates kmers of length k from a string with one step between kmers

    :param dna: string to generate kmers from
    :param k: size of kmer to generate
    """
    assert len(dna) >= 1, "You must select a substring with len(dna) >= 1: {}".format(dna)
    assert k >= 1, "You must select a main_string with k >= 1: {}".format(k)

    for i in range(len(dna)):
        kmer = dna[i:(i + k)]
        if len(kmer) == k:
            yield kmer


def reverse_complement(dna, reverse=True, complement=True):
    """
    Make the reverse complement of a DNA sequence. You can also just make the
    complement or the reverse strand (see options).

    Input: A DNA sequence containing 'ATGC' base pairs and wild card letters

    Output: DNA sequence as a string.

    Options: Specify reverse and/or complement as False to get the complement or
             reverse of the input DNA.  If both are False, input is returned.

    """

    # Make translation table
    trans_table = str.maketrans('ACGTMKRYBVDHNacgtmkrybvdhn',
                                "TGCAKMYRVBHDNtgcakmyrvbhdn")
    # Make complement to DNA
    comp_dna = dna.translate(trans_table)
    # Output all as strings
    if reverse and complement:
        return comp_dna[::-1]
    if reverse and not complement:
        return dna[::-1]
    if complement and not reverse:
        return comp_dna
    if not complement and not reverse:
        return dna


def count_kmers(dna, k):
    """Count all kmers of length k in a string

    :param dna: string to search and count kmers
    :param k: size of kmer
    """
    assert len(dna) >= 1, "You must select a substring with len(dna) >= 1: {}".format(dna)
    assert k >= 1, "You must select a main_string with k >= 1: {}".format(k)
    kmer_count = Counter()
    for i in range(len(dna)):
        kmer = dna[i:(i + k)]
        if len(kmer) == k:
            kmer_count[kmer] += 1
    return kmer_count


def parse_full_alignment_file(alignment_file):
    data = pd.read_csv(alignment_file, usecols=(1, 4, 5, 9, 12, 13), sep="\t",
                         dtype={'ref_pos': np.int64,
                                'strand': np.str,
                                'event_index': np.int64,
                                'kmer': np.str,
                                'posterior_prob': np.float64,
                                'event_mean': np.float64},
                         header=None,
                         names=['ref_pos', 'strand', 'event_index', 'kmer', 'posterior_prob', 'event_mean'])
    return data


class CustomAmbiguityPositions(object):
    def __init__(self, ambig_filepath):
        """Deal with ambiguous positions from a tsv ambiguity position file with the format of
        contig  position            strand  change_from change_to
        'name'  0 indexed position   +/-    C           E


        :param ambig_filepath: path to ambiguity position file"""

        self.ambig_df = self.parseAmbiguityFile(ambig_filepath)

    @staticmethod
    def parseAmbiguityFile(ambig_filepath):
        """Parses a 'ambiguity position file' that should have the format:
            contig  position    strand  change_from change_to

        :param ambig_filepath: path to ambiguity position file
        """
        return pd.read_csv(ambig_filepath, sep="\t",
                             usecols=(0, 1, 2, 3, 4),
                             names=["contig", "position", "strand", "change_from", "change_to"],
                             dtype={"contig": np.str,
                                    "position": np.int,
                                    "strand": np.str,
                                    "change_from": np.str,
                                    "change_to": np.str})

    def getForwardSequence(self, contig, raw_sequence):
        """Edit 'raw_sequence' given a ambiguity positions file. Assumes raw_sequence is forward direction( 5'-3')
        :param contig: which contig the sequence belongs (aka header)
        :param raw_sequence: raw nucleotide sequence
        :return: edited nucleotide sequence
        """
        return self._get_substituted_sequence(contig, raw_sequence, "+")

    def getBackwardSequence(self, contig, raw_sequence):
        """Edit 'raw_sequence' given a ambiguity positions file, Assumes raw_sequence is forward direction( 5'-3')
        :param contig: which contig the sequence belongs (aka header)
        :param raw_sequence: raw nucleotide sequence
        :return: edited nucleotide sequence
        """
        raw_sequence = reverse_complement(raw_sequence, reverse=False, complement=True)
        return self._get_substituted_sequence(contig, raw_sequence, "-")

    def _get_substituted_sequence(self, contig, raw_sequence, strand):
        """Change the given raw nucleotide sequence using the edits defined in the positions file

        :param contig: name of contig to find
        :param raw_sequence: nucleotide sequence (note: this is note edited in this function)
        :param strand: '+' or '-' to indicate strand
        """
        contif_df = self._get_contig_positions(contig, strand)
        raw_sequence = list(raw_sequence)
        for _, row in contif_df.iterrows():
            if raw_sequence[row["position"]] != row["change_from"]:
                raise RuntimeError(
                    "[CustomAmbiguityPositions._get_substituted_sequence]Illegal substitution requesting "
                    "change from %s to %s, row: %s" % (raw_sequence[row["position"]], row["change_to"], row))
            raw_sequence[row["position"]] = AMBIG_BASES["".join(sorted(row["change_to"]))]
        return "".join(raw_sequence)

    def _get_contig_positions(self, contig, strand):
        """Get all unique locations within the positions file

        :param contig: name of contig to find
        :param strand: '+' or '-' to indicate strand
        """
        df = self.ambig_df.loc[
            (self.ambig_df["contig"] == contig) & (self.ambig_df["strand"] == strand)].drop_duplicates()
        assert len(df['position']) == len(set(df['position'])), "Multiple different changes for a single position. {}" \
            .format(df['position'])
        return df


def processReferenceFasta(fasta, work_folder, name, motifs=None, positions_file=None):
    """loops over all of the contigs in the reference file, writes the forward and backward sequences
    as flat files (no headers or anything) for signalMachine, returns a dict that has the sequence
    names as keys and the paths to the processed sequence as keys

    :param fasta: path to un-edited fasta file
    :param work_folder: FolderHandler object
    :param motifs: list of tuple pairs for motif edits. ex [["CCAGG", "CEAGG"]]
    :param positions_file: ambiguous positions file which can be processed via CustomAmbiguityPositions
    :return: paths to possibly edited forward reference sequence and backward reference sequence
    """
    positions = None
    # if no processing needs to happen
    if positions_file is None and motifs is None:
        return fasta, None
    # Cant pass positions file and motifs
    if positions_file is not None and motifs is not None:
        raise RuntimeError("[processReferenceFasta] Cannot specify motif key and ambiguity position file")
    # get positions object (if appropriate)
    if positions_file:
        if not os.path.exists(positions_file):
            raise RuntimeError("[processReferenceFasta] Did not find ambiguity position file here: %s" %
                               positions_file)
        positions = CustomAmbiguityPositions(positions_file)

    # process fasta
    fw_fasta_path = work_folder.add_file_path("forward.{}.{}".format(name, os.path.basename(fasta)))
    bw_fasta_path = work_folder.add_file_path("backward.{}.{}".format(name, os.path.basename(fasta)))
    print("[SignalAlignment.run] NOTICE: Creating forward and backward fasta files.")
    with open(bw_fasta_path, 'w') as bw_outfasta, open(fw_fasta_path, 'w') as fw_outfasta:
        for header, comment, sequence in read_fasta(fasta):
            # signalAlign likes uppercase
            if positions is not None:
                fw_sequence = positions.getForwardSequence(contig=header, raw_sequence=sequence.upper())
                bw_sequence = positions.getBackwardSequence(contig=header, raw_sequence=sequence.upper())
            else:
                fw_sequence = sequence.upper()
                bw_sequence = reverse_complement(fw_sequence, reverse=False, complement=True).upper()
                if motifs:
                    fw_sequence = replace_motifs_sequence_positions(fw_sequence, motifs, True)
                    bw_sequence = replace_motifs_sequence_positions(bw_sequence, motifs, True)

            print(">%s %s\n%s" % (header, "backward", bw_sequence), file=bw_outfasta)
            print(">%s %s\n%s" % (header, "forward", fw_sequence), file=fw_outfasta)

    return fw_fasta_path, bw_fasta_path


def get_full_nucleotide_read_from_alignment(alignment_location, read_name, hardclip_character=None):
    sequence, qualities, hardclipped_start, hardclipped_end = None, None, 0, 0
    with closing(pysam.AlignmentFile(alignment_location, 'rb' if alignment_location.endswith("bam") else 'r')) as aln:
        for aligned_segment in aln.fetch(until_eof=True):
            if read_name not in aligned_segment.qname:
                continue
            BAM_CHARD_CLIP = 5

            # get data and sanity check
            sequence = aligned_segment.query_sequence.upper()
            qualities = aligned_segment.qual
            cigar_tuples = aligned_segment.cigartuples
            if cigar_tuples is None or len(cigar_tuples) == 0:
                print("[get_full_nucleotide_read_from_alignment] no alignment found for {} in {}".format(
                    read_name, alignment_location), file=sys.stderr)
                break

            # check for hard clipping
            if cigar_tuples[0][0] == BAM_CHARD_CLIP:
                hardclipped_start = cigar_tuples[0][1]
                if hardclip_character is not None:
                    sequence = (hardclip_character * hardclipped_start) + sequence
                    if qualities is not None and len(qualities) != 0:
                        qualities = ("!" * hardclipped_start) + qualities
            if cigar_tuples[-1][0] == BAM_CHARD_CLIP:
                hardclipped_end = cigar_tuples[-1][1]
                if hardclip_character is not None:
                    sequence = sequence + (hardclip_character * hardclipped_end)
                    if qualities is not None and len(qualities) != 0:
                        qualities = qualities + ("!" * hardclipped_end)

            # check for reverse mapping
            if aligned_segment.is_reverse:
                sequence = reverse_complement(sequence, reverse=True, complement=True)
                if qualities is not None and len(qualities) != 0:
                    qualities = ''.join(reversed(list(qualities)))
                tmp = hardclipped_end
                hardclipped_end = hardclipped_start
                hardclipped_start = tmp

            # stop looking (assuming only one alignment per read in file)
            break

    return sequence, qualities, hardclipped_start, hardclipped_end, aligned_segment
