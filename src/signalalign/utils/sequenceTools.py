#!/usr/bin/env python3

import string
from collections import Counter
# from sonLib.bioio import fastaRead
import array


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
    trans_table = string.maketrans('ATGCatgc', 'TACGtacg')
    
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
    """count the kmers of length k in a string"""
    kmer_count = Counter()
    for i in range(len(dna)):
        kmer = dna[i:(i+k)]
        if len(kmer) == k:
            kmer_count[kmer] += 1
    return kmer_count


def getFastaDictionary(fastaFile):
    """Returns a dictionary of the first words of fasta headers to their corresponding
    fasta sequence
    """
    namesAndSequences = [(x[0].split()[0], x[1]) for x in fastaRead(open(fastaFile, 'r'))]
    names = [x[0] for x in namesAndSequences]
    assert len(names) == len(set(names)) #Check all the names are unique
    return dict(namesAndSequences) #Hash of names to sequences


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
                    seq.extend(map(ord, line))
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
        fileHandle.write("%s\n" % seq[i:i+chunkSize])
    if isinstance(fileHandleOrFile, "".__class__):
        fileHandle.close()
