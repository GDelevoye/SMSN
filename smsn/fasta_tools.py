#!/usr/bin/env python

__author__ = "DELEVOYE Guillaume"
__credits__ = ["BAHIN Mathieu", "MEYER Eric"]
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
    cpl = {"A": "T", "C": "G", "T": "A", "G": "C", "N": "N"}
    return "".join([cpl[x] for x in seq[::-1]])

def get_snipet(seq, position, n, strand, padding_eachside_upto=12):
    """Returns the context around n bases (and adds padding if at the extremity of the sequence)
    'position' msut be 0-based (seq is considered as a python string)"""
    assert n >= 1
    assert position >= 0 and position <= len(seq)
    assert seq
    assert padding_eachside_upto > 0
    assert strand in [0, 1]

    left = seq[max(0, position - n):position]
    right = seq[position + 1:min(position + n, len(seq)) + 1]

    left = left.rjust(padding_eachside_upto, "N")
    right = right.ljust(padding_eachside_upto, "N")

    returnseq = left + seq[position] + right

    if strand == 0:
        return returnseq
    else:
        return reverse_cpl(returnseq)


def dict_to_fasta(myDict, fastafile, specify_HoleID = False):
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
            if specify_HoleID:
                filout.write('>'+'HoleID_'+str(key)+'\n'+str(divide_seq(myDict[key]))+'\n')
            else:
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


def compute_chunk_infos(real_start, real_end, fasta_this_scaffold):
    """Creates a chunk to build a false ref at +100 / -100 of the begin/end where the CCS mapped"""
    if real_start - 100 < 0:
        chunk_start = 0
    else:
        chunk_start = real_start - 100
    if real_end + 100 > len(fasta_this_scaffold):
        chunk_end = len(fasta_this_scaffold)
    else:
        chunk_end = real_end + 100

    chunk_size = chunk_end - chunk_start
    offset = real_start - 100
    sequence = fasta_this_scaffold[chunk_start:chunk_end]

    return chunk_start, chunk_end, chunk_size, offset, sequence