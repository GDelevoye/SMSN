"""Tools to handle the .fasta files and the .gff annotations"""

__author__ = "DELEVOYE Guillaume"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "DELEVOYE Guillaume"
__email__ = "delevoye@ens.fr"
__status__ = "Production"

import os
import pandas as pd
import logging
import typing


def parse_gff(gfffile: os.path.realpath, parse_attributes: bool = True) -> pd.DataFrame:
    """Parses a .gff file and returns a pandas dataframe.
    In parse_attributes is set to True, separated columns will be returned for each attribute
    (in addition to the "attributes" column itself)"""

    # Not really a .fasta tool but there's no better place than here

    logging.debug('[DEBUG] Parsing gff {}'.format(gfffile))
    gfffile = os.path.realpath(gfffile)
    gfffile = open(gfffile, "r")

    list_interest = ["scaffold", "source", "type", "start",
                     "end", "score", "strand", "phase", "attributes"]  # The base files

    list_attributes = set()

    proto_df = []
    for line in gfffile:
        if line[0] == "#":
            continue

        line_dict = {}
        line = line.strip().split()
        i = 0

        for elt in list_interest:
            line_dict[elt] = line[i]
            i += 1

            if elt == "attributes" and parse_attributes:
                miniline = line[-1].strip()

                for subelt in miniline.split(';'):  # Parse the attribute column, which is ";" separated
                    key = subelt.split('=')[0]
                    value = subelt.split('=')[1]
                    line_dict[key] = value
                    list_attributes.add(key)

        if line_dict:
            proto_df.append(line_dict)

    if parse_attributes:
        list_interest = list_interest + sorted([x for x in list_attributes])

    if len(proto_df) > 0:
        return pd.DataFrame(proto_df)[list_interest]
    else:
        return pd.DataFrame(columns=list_interest)


def load_fasta(fastafile: os.path.realpath) -> typing.Dict[str, str]:
    """Returns a python dict { id : sequence } for the given .fasta file"""
    logging.debug("[DEBUG] Loading fastafile {}".format(fastafile))
    with open(os.path.realpath(fastafile), 'r') as filin:
        fasta = filin.read()
        fasta = fasta.split('>')[1:]
        outputdict = {
            x.split('\n')[0].strip(): "".join(x.split('\n')[1:]) for x in fasta}
    return outputdict


def load_fasta_special(fastafile: os.path.realpath) -> typing.Dict[str,str]:
    """Returns a python dict { id : sequence } for the given .fasta file
    BUT as the name warns it s a somewhat special way:
    It replaces whitespaces and separator special characters by "_"
    This is needed since PacBio systematically ignores the whitespace characters and tabulations might break
    the .bam format."""

    def handle_characters(name: str) -> str:
        """Spots special characters in the .fasta file that might trigger an error in the PacBio tools.
        In case it does, the string is returned corrected"""
        for elt in ["\\t"," ", ".", ",", ",", ",", "-", "&", "|", "[", "]",
                    "{", "}", "/"]:
            name = name.replace(elt, "_")
        return name

    return_dict = {}

    with open(os.path.realpath(fastafile), 'r') as filin:
        fasta = filin.read()
        # If the fasta starts with an ">" like it should, then the first split
        # is empty ''
        fasta = fasta.split('>')[1:]
        # so we take [1:] instead of all the splitted list

        for line in fasta:
            scaffoldname = line.split('\n')[0].strip()
            scaffoldname = handle_characters(scaffoldname)

            fastaseq = "".join(line.split('\n')[1:]).upper()
            return_dict[scaffoldname] = fastaseq

    return return_dict


def reverse_cpl(seq: str) -> str:
    """
    From a DNA sequence (iterable : string, list, etc...) returns the reverse complementary

    >>> reverse_cpl("GATC")
    'GATC'
    >>> reverse_cpl("AT")
    'AT'
    >>> reverse_cpl("GNGCTAGN")
    'NCTAGCNC'
    """

    cpl = {"A": "T", "C": "G", "T": "A", "G": "C", "N": "N"}
    return "".join([cpl[x] for x in seq.upper()[::-1]])


def get_snipet(seq: str, position: int, strand: int = 0, n: int = 12) -> str:
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
    assert 0 <= position <= len(seq)
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


def dict_to_fasta(myDict: typing.Dict[str, str], fastafile: os.path.realpath) -> None:
    """ from a python dict {Identifier:sequence}, writes a fastafile
    If "specify_HoleID is true, then Identifier is an int and therefore is
    preceeded by a "Hole_ID" in the .fasta. This is meant to avoid having
    int-starting identifiers, which is forbidden in many .fasta conventions"""
    logging.debug('[DEBUG] (dict to fasta), fastafile = {}'.format(fastafile))

    def divide_seq(seq, length=60, sep='\n'):
        """From a full-length sequence, divides it in chunks of 60 letters (most frequent.fasta format)"""
        return str(sep).join([seq[x:x + int(length)] for x in range(0, len(seq), length)])
        # Tricky one liner --> Splits a sequence into chunks of size "length"

    with open(os.path.realpath(fastafile), "w") as filout:
        first = True
        for key in myDict:
            if not first:
                filout.write('\n')
            else:
                first = False
            filout.write('>' + str(key) + '\n' +
                         str(divide_seq(myDict[key])))


def compute_chunk_infos(real_start: int,
                        real_end: int,
                        fasta_this_scaffold: os.path.realpath,
                        distance: int = 100) -> typing.Tuple[int, int, int, int, str]:
    """Creates a chunk to build a false ref at +100 / -100 of the begin/end where the CCS mapped
    Coordinate system is 0-based and is half-open [start;end)
    This means that if you consider the "real_end" as the last position where a
    nucleotide mapped, you must pass real_end+1 to the function."""

    logging.debug(
        "[DEBUG] (compute_chunk_infos) real_start = {}, real_end = {}, distance = {}".format(
            real_start,
            real_end,
            distance))

    assert (real_end - real_start) <= len(fasta_this_scaffold) + 1
    assert real_end <= len(fasta_this_scaffold) + 1

    if real_start - distance < 0:
        chunk_start = 0
    else:
        chunk_start = real_start - distance
    if real_end + distance > len(fasta_this_scaffold):
        chunk_end = len(fasta_this_scaffold)
    else:
        chunk_end = real_end + distance

    chunk_size = chunk_end - chunk_start
    offset = max(0, real_start - distance)
    sequence = fasta_this_scaffold[chunk_start:chunk_end]

    logging.debug(
        '[DEBUG] (compute_chunk_infos) returning chunk_start = {},chunk_end = {}, chunk_size = {},'
        'offset = {}, sequence = {}'.format(
            chunk_start,
            chunk_end,
            chunk_size,
            offset,
            sequence))

    return chunk_start, chunk_end, chunk_size, offset, sequence
