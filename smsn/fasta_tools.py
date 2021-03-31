#!/usr/bin/env python

__author__ = "DELEVOYE Guillaume"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "DELEVOYE Guillaume"
__email__ = "delevoye@ens.fr"
__status__ = "Developpment"

import os
import pandas as pd
import logging

def parse_gff(gfffile): # Not really a .fasta tool but there's no better place than here

    logging.debug('[DEBUG] Parsing gff {}'.format(gfffile))
    gfffile = os.path.realpath(gfffile)
    gfffile = open(gfffile, "r")

    proto_df = []
    for line in gfffile:
        if line[0] == "#":
            continue

        line_dict = {}
        line = line.split('\t')
        i = 0

        for elt in ["scaffold", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]:
            if elt == "attributes":
                miniline = line[-1].strip()
                for subelt in miniline.split(';'):
                    key = subelt.split('=')[0]
                    value= subelt.split('=')[1]
                    line_dict[key] = value
                if "identificationQv" not in miniline:
                    key = "identificationQv"
                    value = "NaN"
                    line_dict[key] = value
            else:
                line_dict[elt] = line[i]
                i += 1
        if line_dict:
            proto_df.append(line_dict)

    if len(proto_df)>0:
        return pd.DataFrame(proto_df)[
            ["scaffold", "source", "type", "start", "end", "score", "strand", "phase","coverage","context","IPDRatio","identificationQv"]]
    else:
        logging.error('[ERROR] gfffile {} was apparently empty'.format(gfffile))
        return -1


def load_fasta(fastafile):
    """Returns a python dict { id : sequence } for the given .fasta file"""
    logging.debug("[DEBUG] Loading fastafile {}".format(fastafile))
    with open(os.path.realpath(fastafile), 'r') as filin:
        fasta = filin.read()
        fasta = fasta.split('>')[1:]
        outputdict = {x.split('\n')[0].strip(): "".join(x.split('\n')[1:]) for x in fasta}
    return outputdict


def load_fasta_special(fastafile):
    """Returns a python dict { id : sequence } for the given .fasta file
    BUT as the name warns it s a SPECIAL way: it returns identifiers in the form "_".join(identifier.split(' '))
    This is needed since PacBio systematically ignores the " " characters"""

    def handle_characters(name):
        for elt in [" ","#",".",",",",",",","*","-","&","|","[","]","{","}","\"","\'","~","%","+",">","<","!","?","/","\\","é","è",":"]:
            name = name.replace(elt,"_")
        return name

    return_dict = {}

    with open(os.path.realpath(fastafile), 'r') as filin:
        fasta = filin.read()
        fasta = fasta.split('>')[1:] #If the fasta starts with an ">" like it should, then the first split is empty ''
        # so we take [1:] instead of all the splitted list

        for line in fasta:
            scaffoldname = line.split('\n')[0].strip()
            scaffoldname = handle_characters(scaffoldname)

            fastaseq = "".join(line.split('\n')[1:])
            return_dict[scaffoldname] = fastaseq

    return return_dict


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
    logging.debug('[DEBUG] (dict to fasta), fastafile = {}'.format(fastafile))
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


def compute_chunk_infos(real_start, real_end, fasta_this_scaffold, distance=100):
    """Creates a chunk to build a false ref at +100 / -100 of the begin/end where the CCS mapped"""

    logging.debug("[DEBUG] (compute_chunk_infos) real_start = {}, real_end = {}, distance = {}".format(real_start,real_end,distance))

    if real_start - distance < 0:
        chunk_start = 0
    else:
        chunk_start = real_start - distance
    if real_end + distance > len(fasta_this_scaffold):
        chunk_end = len(fasta_this_scaffold)
    else:
        chunk_end = real_end + distance


    chunk_size = chunk_end - chunk_start
    offset = max(0,real_start - distance)
    sequence = fasta_this_scaffold[chunk_start:chunk_end]

    logging.debug('[DEBUG] (compute_chunk_infos) returning chunk_start = {},chunk_end = {}, chunk_size = {},offset = {}, sequence = {}'.format(chunk_start,chunk_end,chunk_size,offset,sequence))

    return chunk_start, chunk_end, chunk_size, offset, sequence