#!/usr/bin/env python

__author__ = "DELEVOYE Guillaume"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "DELEVOYE Guillaume"
__email__ = "delevoye@ens.fr"
__status__ = "Developpment"

import os

def load_fasta(fastafile):
    """Returns a python dict { id : sequence } for the given .fasta file"""
    with open(os.path.realpath(fastafile), 'r') as filin:
        fasta = filin.read()
        fasta = fasta.split('>')[1:]
        outputdict = {x.split('\n')[0].strip(): "".join(x.split('\n')[1:]) for x in fasta}
    return outputdict

def load_fasta_special(fastafile):
    """Returns a python dict { id : sequence } for the given .fasta file
    BUT as the name warns it s a SPECIAL way: Only the fist part of the identifier will be taken """
    with open(os.path.realpath(fastafile), 'r') as filin:
        fasta = filin.read()
        fasta = fasta.split('>')[1:]
        outputdict = { x.split('\n')[0].split()[0] : "".join(x.split('\n')[1:]) for x in fasta}
    return outputdict

def reverse_cpl(seq):
    """
    >>> reverse_cpl("GATC")
    'GATC'
    >>> reverse_cpl("AT")
    'AT'
    >>> reverse_cpl("GNGCTAGN")
    'NCTAGCNC'
    """

    cpl = {"A": "T", "C": "G", "T": "A", "G": "C", "N": "N"}
    return "".join([cpl[x] for x in seq[::-1]])

def get_snipet(seq, position, strand, n=12):
    """Returns the context around +/- n bases (and adds padding if at the extremity of the sequence)
    'position' msut be 0-based (seq is considered as a python string)
    n must be even

    >>> testseq = "GATCGCT"
    >>> get_snipet(testseq,0,0,2)
    'NNGAT'
    >>> get_snipet(testseq,position=1,strand=1,n=3)
    'CGATCNN'
    >>> get_snipet(testseq,position=6,strand=1,n=3)
    'NNNAGCG'
    >>> get_snipet(testseq,position=3,strand=0,n=2)
    'ATCGC'
    >>> testseq = "NNGTAGCNGTCAATCG"
    >>> get_snipet(testseq,position=3,strand=0,n=12)
    'NNNNNNNNNNNGTAGCNGTCAATCG'
    >>> reverse = reverse_cpl(testseq)
    >>> reverse
    'CGATTGACNGCTACNN'
    >>> get_snipet(testseq,position=1,strand=1,n=1)
    'CNN'
    >>> get_snipet(testseq,position=1,strand=1,n=3)
    'TACNNNN'


    """
    assert n >= 1
    assert position >= 0 and position <= len(seq)
    assert seq
    assert strand in [0, 1]

    seq = seq.upper()

    left = seq[max(0, position - n):position]
    right = seq[position + 1:min(position + n, len(seq)) + 1]

    left = left.rjust(n, "N")
    right = right.ljust(n, "N")

    returnseq = left + seq[position] + right

    if strand == 0:
        return returnseq
    else:
        return reverse_cpl(returnseq)


def dict_to_fasta(myDict, fastafile):
    """ from a python dict {Identifier:sequence}, writes a fastafile
    If "specify_HoleID is true, then Identifier is an int and therefore is
    preceeded by a "Hole_ID" in the .fasta. This is meant to avoid having
    int-starting identifiers, which is forbidden in many .fasta conventions"""
    def divide_seq(seq, length = 60, sep = '\n'):
        """From a full-length sequence, divides it in chunks of 60 letters (most frequent.fasta format)"""
        return( str(sep).join([seq[x:x+int(length)] for x in range(0, len(seq), length)]))
        # Tricky one liner --> Splits a sequence into chunks of size "length"
    with open(os.path.realpath(fastafile), "w") as filout:
        for key in myDict:
            filout.write('>'+str(key)+'\n'+str(divide_seq(myDict[key]))+'\n')

# # Test the code above
# test1 = "ATGCTAGCTTTTTCGTGAGGCTGATCGATATATATATAGCGGCGCTAGCTGATCGATTAGCTA"
# test2 = "GCGCGATCGATCGATGCTGTATATTAGCGATGCTAGCTCGATTCTCTAGAGATCGATCGATTACGGTAGC"
# fastadict = {"scaffoldtest1":test1,"scaffoldtest2":test2}
# dict_to_fasta(fastadict,"test.fasta")
# myfa = load_fasta("test.fasta")
# get_snipet(seq=myfa["scaffoldtest1"],
#            position=2,
#            n=6,
#            strand=0)


def compute_chunk_infos(real_start, real_end, fasta_this_scaffold,distance=100): #TODO #FIXME Maybe one day replace it by get_snippet that is way more generalistic
    """Creates a chunk to build a false ref at +100 / -100 of the begin/end where the CCS mapped"""
    if real_start - distance < 0:
        chunk_start = 0
    else:
        chunk_start = real_start - distance
    if real_end + distance > len(fasta_this_scaffold):
        chunk_end = len(fasta_this_scaffold)
    else:
        chunk_end = real_end + distance

    chunk_size = chunk_end - chunk_start
    offset = real_start - distance
    sequence = fasta_this_scaffold[chunk_start:chunk_end]

    return chunk_start, chunk_end, chunk_size, offset, sequence