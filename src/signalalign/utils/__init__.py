
from __future__ import print_function
import os

from collections import Counter

import numpy as np
import pandas as pd

from signalalign.motif import getMotif
from signalalign.utils.parsers import read_fasta


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
    data = pd.read_table(alignment_file, usecols=(1, 4, 5, 9, 12, 13),
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
        return pd.read_table(ambig_filepath,
                             usecols=(0, 1, 2, 3, 4),
                             names=["contig", "position", "strand", "change_from", "change_to"],
                             dtype={"contig"      : np.str,
                                    "position"    : np.int,
                                    "strand"      : np.str,
                                    "change_from" : np.str,
                                    "change_to"   : np.str})

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
                raise RuntimeError("[CustomAmbiguityPositions._get_substituted_sequence]Illegal substitution requesting "
                                   "change from %s to %s, row: %s" % (raw_sequence[row["position"]], row["change_to"], row))
            raw_sequence[row["position"]] = row["change_to"]
        return "".join(raw_sequence)

    def _get_contig_positions(self, contig, strand):
        """Get all unique locations within the positions file

        :param contig: name of contig to find
        :param strand: '+' or '-' to indicate strand
        """
        df = self.ambig_df.loc[(self.ambig_df["contig"] == contig) & (self.ambig_df["strand"] == strand)].drop_duplicates()
        assert len(df['position']) == len(set(df['position'])), "Multiple different changes for a single position. {}"\
            .format(df['position'])
        return df


def processReferenceFasta(fasta, work_folder, motif_key=None, sub_char=None, positions_file=None):
    """loops over all of the contigs in the reference file, writes the forward and backward sequences
    as flat files (no headers or anything) for signalMachine, returns a dict that has the sequence
    names as keys and the paths to the processed sequence as keys

    :param fasta: path to un-edited fasta file
    :param work_folder: FolderHandler object
    :param motif_key:
    :param sub_char:
    :param positions_file: ambiguous positions file which can be processed via CustomAmbiguityPositions
    :return: paths to possibly edited forward reference sequence and backward reference sequence
    """
    # argument sanitization
    if positions_file is not None and motif_key is not None:
        raise RuntimeError("[processReferenceFasta] Cannot specify motif key and ambiguity position file")
    if positions_file is not None and sub_char is not None:
        raise RuntimeError("[processReferenceFasta] Cannot specify a substitution character and an ambiguity position file")

    # get positions object (if appropriate)
    if positions_file is not None:
        if not os.path.exists(positions_file):
            raise RuntimeError("[processReferenceFasta] Did not find ambiguity position file here: %s" %
                               positions_file)
        positions = CustomAmbiguityPositions(positions_file)
    else:
        positions = None

    if positions_file is None and motif_key is None and sub_char is None:
        return fasta, None
    # process fasta
    fw_fasta_path = work_folder.add_file_path("forward.{}".format(os.path.basename(fasta)))
    bw_fasta_path = work_folder.add_file_path("backward.{}".format(os.path.basename(fasta)))
    print("[SignalAlignment.run] NOTICE: Creating forward and backward fasta files.")
    with open(bw_fasta_path, 'w') as bw_outfasta, open(fw_fasta_path, 'w') as fw_outfasta:
        for header, comment, sequence in read_fasta(fasta):
            # the motif label allows us to make multiple copies of the reference with unique file names
            # motif_lab = "" if motif_key is None else "%s." % motif_key
            # these are the paths to the flat files that have the references
            # signalAlign likes uppercase
            # TODO make motifs work the way I built them in SequenceTools.py
            if motif_key is not None:
                motif, ok = getMotif(motif_key, sequence)
                if not ok:
                    raise RuntimeError("[processReferenceFasta]Illegal motif key %s" % motif_key)
                fw_sequence = motif.forwardSubstitutedSequence(sub_char).upper()
                bw_sequence = motif.complementSubstitutedSequence(sub_char).upper()
            elif positions is not None:
                fw_sequence = positions.getForwardSequence(contig=header, raw_sequence=sequence.upper())
                bw_sequence = positions.getBackwardSequence(contig=header, raw_sequence=sequence.upper())
            else:
                fw_sequence = sequence.upper()
                bw_sequence = reverse_complement(fw_sequence, reverse=False, complement=True).upper()

            print(">%s %s\n%s" % (header, "backward", bw_sequence), file=bw_outfasta)
            print(">%s %s\n%s" % (header, "forward", fw_sequence), file=fw_outfasta)

    return fw_fasta_path, bw_fasta_path
