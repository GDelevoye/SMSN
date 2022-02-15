"""Tools to handle the PacBio .bam and .sam formats"""

__author__ = "DELEVOYE Guillaume"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "DELEVOYE Guillaume"
__email__ = "delevoye@ens.fr"
__status__ = "Production"

import numpy as np
import logging
import subprocess
import os
import pandas as pd
import smsn
import subprocess

import typing


def check_sorted(bamfilePath: os.path.realpath) -> bool:
    """Checks that in a PacBio subreads .bam file, the subreads' holeID are not mixed"""
    encountered = set()
    is_sorted = True
    first = True

    for line in smsn.bam_toolbox.read_bam(bamfilePath):
        current_holeID = int(line.split()[0].split('/')[1])

        if first:
            encountered.add(current_holeID)
            first = False
            previous_hole_id = current_holeID

        elif current_holeID not in encountered and previous_hole_id != current_holeID:
            encountered.add(current_holeID)
            previous_hole_id = current_holeID

        elif (previous_hole_id != current_holeID) and (current_holeID in encountered):
            is_sorted = False
            return is_sorted

    return is_sorted


def read_bam(bamfilePath: os.path.realpath) -> typing.Generator[str, None, None]:
    """Reads a .bam line by line, as (stripped) .sam text lines (ignores header)"""
    logging.debug(
        '[DEBUG] (read_bam) Launching samtools in subprocess to yield the lines of the bamfile {}'.format(bamfilePath))
    # Samtools must be installed !
    cmd = 'samtools view ' + os.path.realpath(bamfilePath)
    logging.debug('[DEBUG] (read_bam) cmd = {}'.format(cmd))

    with subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE) as proc:

        logging.debug(
            '[DEBUG] (read_bam) Starting to read output from samtools for {}'.format(bamfilePath))

        while True:
            # We will read the output of samtools line by line. Decoding is needed
            # because popen returns bytes, not strings
            line = proc.stdout.readline().decode('utf-8')
            if line:  # While there is something to read (else it will break)
                yield line.strip()
            else:  # There is nothing left in the stdout of samtools
                logging.debug(
                    '[DEBUG] (read_bam) reading bamfile {} seems over'.format(bamfilePath))
                break

        proc.stdout.close()
        return


def read_header(bamfilePath: os.path.realpath) -> str:
    """ returns the header as a python-string .sam"""

    logging.debug("[DEBUG] Reading headers of {}".format(bamfilePath))

    samtxt = []
    cmd = 'samtools view -H ' + str(bamfilePath)
    completed = subprocess.run(["samtools","view","-H",str(bamfilePath)],
                               capture_output=True,
                               check=True,
                               encoding="utf-8")
    if completed.returncode != 0:
        logging.error("[ERROR] {}".format(completed.stderr))
    else:
        return completed.stdout.strip()




def get_hole_id(samseq: str) -> int:
    """Returns the holeID of a .sam string, as an integer) -- Considers that all reads/subreads have the same id"""
    holeID = "".join([x for x in samseq.split('\n') if '@' not in x]).split()[0].split('/')[1]
    return int(holeID)


def sort_by_name(bamfilepath: os.path.realpath,
                 newpath: os.path.realpath,
                 nbcore: int = 1) -> None:
    """ Will sort the bam by read name (it will do it numerically for the numbers,
    not alphabetically, which is convenient for us)"""
    bampath = os.path.realpath(bamfilepath)
    bamdir = os.path.dirname(bampath)
    newpath = os.path.realpath(newpath)
    # HERE = os.getcwd()
    # See https://stackoverflow.com/questions/36506513/pytest-conftest-py-incorrectly-shows-the-directory-os-getcwd-in-pycharm
    HERE = os.path.dirname(os.path.abspath(__file__))
    os.system('mkdir -p ' + str(bamdir))
    os.chdir(bamdir)
    cmd = 'samtools sort -n ' + \
          str(os.path.basename(bamfilepath)) + ' ' + \
          '-o ' + str(newpath) + ' -@ ' + str(nbcore)
    os.system(cmd)
    os.chdir(HERE)


def inmemory_asbam(samstring: str) -> typing.ByteString:
    """ Returns a python byte-encoded .bam from a .sam string.
    Everything is done in memory, not on hard drive.
    Can crash if the RAM is insufficient"""

    p = subprocess.run('samtools view -bS -h', shell=True,
                       stdout=subprocess.PIPE, input=samstring.encode('utf-8'))
    return (p.stdout)


def yield_hole_by_hole(subreads_bamfile: os.path.realpath,
                       include_header: bool = True,
                       restrictionlist=None,
                       min_subreads: int = 0) -> typing.Generator[str, None, None]:
    """ Yields the subreads HoleID by HoleID as a .sam textstring subreads_bamfile must be already
    sorted by HoleID. "restricted_to" can be the exhaustive set of HoleID (int) to include.
    If include_header is set to True, will return the header each time for each hole
    Can also exclude the holes with a critically low number of subreads (0 = No threshold)
    """

    # #TODO This is not a real todo but just a remark if you want to change anything:
    # --> DO NOT read the header each time it's killing everything !

    if not restrictionlist:
        restrictionlist = set()  # Just a trick to avoid a default argument that is mutable
        # https://stackoverflow.com/questions/41686829/warning-about-mutable-default-argument-in-pycharm

    reader = smsn.bam_toolbox.read_bam(subreads_bamfile)  # An iterator
    header = read_header(subreads_bamfile)
    setEncountered = set()
    current_yield = []
    first = True  # Simple flag to treat the first hole specifically

    while True:

        try:
            current_line = next(reader).strip()
            current_holeID = int(current_line.split()[0].split('/')[1])

            if (current_holeID not in restrictionlist) and (len(restrictionlist) > 0):  # We don't care about these
                continue

            elif first:
                setEncountered.add(current_holeID)
                current_yield.append(current_line)
                first = False

            else:
                if current_holeID not in setEncountered:
                    if len(current_yield) > min_subreads:
                        if include_header:
                            yield header.strip() + "\n" + "\n".join(current_yield)
                        else:
                            yield "\n".join(current_yield)
                        current_yield = []
                    else:  # What we stocked until there appears to be insufficiently covered
                        current_yield = []

                    setEncountered.add(current_holeID)
                current_yield.append(current_line)


        except StopIteration:  # Return the last hole if there is one
            if len(current_yield) > min_subreads:
                if include_header:
                    yield header.strip() + "\n" + "\n".join(current_yield)
                else:
                    yield "\n".join(current_yield)
            break

    logging.debug('[INFO] Yielding .bam file {} hole : over'.format(subreads_bamfile))
    return


def parse_header(bamfile: os.path.realpath) -> typing.Dict[str,typing.Union[str,float]]:
    """Returns as a dict the crucial informations that interest us in the header
     (ID, PL, PU, PM, READTYPE, Ipd:CodecV1, FRAMERATEHZ, etc"""

    logging.debug("[DEBUG] Parsing informations from the .bam input")

    splitted_header = read_header(bamfile).split('\n')
    lineRG = [x for x in splitted_header if x[0:3] == "@RG"][0]
    predict = lineRG.split()[1:]

    outputdict = {}
    for elt in predict:
        if elt[0:2] != "DS":
            key = elt.split(':')[0]
            value = elt.split(':')[1]
            outputdict[key] = value

    for elt in predict:
        if elt[0:2] == "DS":
            DSline = elt[3:]

    for elt in DSline.split(';'):
        key = elt.split('=')[0]
        value = elt.split('=')[1]
        outputdict[key] = value

    if "Ipd:Frames" in outputdict:
        logging.warning("[WARNING] Your data seems to contain raw frames of inter-pulse durations (IPDs), rather than "
                        "the usual lossy-encoding (CodecV1). See : "
                        "https://pacbiofileformats.readthedocs.io/en/3.0/BAM.html. The program will continue anyway, "
                        "and might work properly. Be carefull, however, when analyzing your data because the pipeline "
                        "was never tested without CodecV1. If your data was generated from old h5 files, "
                        "it is possible to re-use bax2bam to encode the IPDs with the CodecV1.")

    outputdict["FRAMERATEHZ"] = float(outputdict["FRAMERATEHZ"])
    assert outputdict["FRAMERATEHZ"] > 0

    return outputdict


# def parse_line(line, include_ipds=False):
#     line = line.split('\t')
#
#     returndict = {}
#     returndict["alignment_flag"] = int(line[1])
#     returndict["scaffold"] = line[2]
#     returndict["start"] = line[3]
#     # A CIGAR object will be created later to avoid saturating the RAM
#     returndict["CIGAR_string"] = line[5]
#     returndict["HoleID"] = int(line[0].split('/')[1])
#     returndict["sequence"] = line[9]
#     returndict["moviename"] = line[0].split('/')[0]
#
#     if include_ipds:
#         returndict["encoded_ipd"] = [x for x in line if x[:3] == "ip:"][0]
#         # This should then be decoded, but using np.uint8 ensures that not too
#         # much place is taken in RAM
#         returndict["encoded_ipd"] = np.array(
#             returndict["encoded_ipd"].split(',')[
#                 1:], dtype=np.uint8)
#     return returndict


def detailed_get_pd_dataframe_alignment(bamfile: os.path.realpath) -> pd.DataFrame:
    """ From a pacbio .bam alignment file produced by blasr, and gives a pandas dataframe corresponding to the
    alignment"""

    list_alignments = []  # We will return a dataframe from a list of dicts
    # In each dict, the key is the name of the feature while the value is
    # its... value (Wow !)
    for line in read_bam(
            os.path.realpath(bamfile)):  # Reads the .bam line after line
        if line:
            line = line.split('\t')
            moviename = line[0].split('/')[0]
            scaffold = line[2]
            start = int(line[3])
            CIGAR_string = line[5]

            sequence = line[9]

            if CIGAR_string != "*":
                CIGAR_obj = CIGAR(
                    CIGAR_string)
            else:  # Happens for unaligned reads
                CIGAR_string = str(len(sequence)) + "S"
                CIGAR_obj = CIGAR(CIGAR_string)

            identity_score = CIGAR_obj.identity
            # Reflen is the size of reference deduced from the CIGAR string
            end = start + CIGAR_obj.reflen
            # Because the alignment starting position doesn't take the soft-clipping in account,
            # We don't take it in account either (reflen doesn't take clipped
            # bases in account)

            # But clipped bases are counted as unaligned for %identity
            # computation !

            reflen = end - start

            HoleID = int(line[0].split('/')[1])
            alignment_flag = int(line[1])

            mapQV = int(line[4])
            clipped_bases = CIGAR_obj.number_of_clipped
            matching_bases = CIGAR_obj.number_of_matches

            list_interest = ["moviename", "scaffold", "start", "CIGAR_string", "identity_score", "end", "HoleID",
                             "alignment_flag", "reflen", "mapQV", "matching_bases", "clipped_bases", "sequence"]

            alignment_dict = {}
            for elt in list_interest:
                # We add all these values to a dict
                alignment_dict[elt] = eval(elt)

            list_alignments.append(alignment_dict)
            del line
    # We return a list of dict transformed in a pandas DataFrame
    return pd.DataFrame(list_alignments)


class CIGAR:
    """Class to parse the alignment CIGAR

    >>> cigarstring = "3=1I1D1X1="
    >>> cigObject = CIGAR(cigarstring)
    >>> cigObject.number_of_deletions
    1
    >>> cigObject.number_of_clipped
    0
    >>> cigObject.number_of_insertions
    1
    >>> cigObject.number_of_matches
    4
    >>> cigObject.seqlen
    6
    >>> cigObject.reflen
    6
    >>> cigObject.identity
    0.6666666666666666
    >>> cigObject2 = CIGAR("1=1X1D1I3=")
    >>> cigObject2 != cigObject
    True
    >>> smsn.bam_toolbox.CIGAR("ecoli")
    Traceback (most recent call last):
     ...
    ValueError: invalid literal for int() with base 10: ''
    """

    def __init__(self, strCIGAR: str, scaffold=None):
        self.string = str(strCIGAR)

        # List of tuples [ (32 , X), (7, =) ... ]
        self.listgroups = self._compute_listgroups()
        # Python list of True and False, corresponding do Match/Unmatch
        self.binarlist = self._binarize()
        self.matchArray = self._transform_asNumpyMatches()  # A numpy array
        # I think it's a name explicit enough ?
        self.number_of_matches = np.sum(self.matchArray)

        # Python list of True (Deletion)
        self.deletionbinarlist = self._binarize_deletions()
        self.deletionArray = self._transform_deletions_asNumpy()  # A numpy array
        self.number_of_deletions = np.sum(self.deletionArray)

        self.insertionbinarlist = self._binarize_insertions()
        self.insertionArray = self._transform_insertions_asNumpy()
        self.number_of_insertions = np.sum(self.insertionArray)

        self.substitutionbinarlist = self._binarize_substitutions()
        self.substitutionsArray = self._transform_substitutions_asNumpy()
        self.number_of_substitutions = np.sum(self.substitutionsArray)

        self.number_of_clipped = self.get_clip()
        self.scaffold = scaffold

        self.reflen = self.number_of_matches + \
                      self.number_of_substitutions + self.number_of_deletions
        self.seqlen = self.number_of_matches + self.number_of_clipped + \
                      self.number_of_substitutions + self.number_of_insertions
        self.identity = self.number_of_matches * \
                        2 / (self.reflen + self.seqlen)

    def _compute_listgroups(self):
        listgroups = []
        for group in self._enumerate_groups():
            listgroups.append(group)

        return listgroups

    def _enumerate_groups(self):
        tmp_string = ''

        for elt in self.string:
            if elt.isnumeric():
                tmp_string = tmp_string + str(elt)
            else:
                yield (int(tmp_string), str(elt))
                tmp_string = ''

    def _binarize(self):
        binarlist = []
        for elt in self.listgroups:
            if elt[1] == '=':
                for i in range(0, elt[0]):
                    binarlist.append(True)
            else:
                for i in range(0, elt[0]):
                    binarlist.append(False)
            del i

        return binarlist

    def _transform_asNumpyMatches(self):
        """ From the CIGAR string, returns a numpy array of True (Match) and False (Mismatch/Indel etc) base/base """
        return np.array(self.binarlist, dtype=bool)

    def _binarize_deletions(self):
        binarlist = []
        for elt in self.listgroups:
            if elt[1] == 'D':
                for i in range(0, elt[0]):
                    binarlist.append(True)
                del i
            else:
                for i in range(0, elt[0]):
                    binarlist.append(False)
                del i

        return binarlist

    def _transform_deletions_asNumpy(self):
        return np.array(self.deletionbinarlist)

    def _binarize_insertions(self):
        binarlist = []
        for elt in self.listgroups:
            if elt[1] == 'I':
                for i in range(0, elt[0]):
                    binarlist.append(True)
                del i
            else:
                for i in range(0, elt[0]):
                    binarlist.append(False)
                del i

        return binarlist

    def _transform_insertions_asNumpy(self):
        return np.array(self.insertionbinarlist)

    def _binarize_substitutions(self):
        binarlist = []
        for elt in self.listgroups:
            if elt[1] == 'X':
                for i in range(0, elt[0]):
                    binarlist.append(True)
                del i
            else:
                for i in range(0, elt[0]):
                    binarlist.append(False)
                del i

        return binarlist

    def _transform_substitutions_asNumpy(self):
        return np.array(self.substitutionbinarlist)

    def _enumerate_clipping(self):
        tmp_string = ''

        for elt in self.string:
            if elt.isnumeric():
                tmp_string = tmp_string + str(elt)
            else:
                if elt == 'S':
                    yield (int(tmp_string), str(elt))
                    tmp_string = ''
                else:
                    tmp_string = ''

    def get_clip(self):
        """ Returns the number of soft-clipped bases (mainly for debugging purposes)"""
        count = 0
        for clipped_bases in self._enumerate_clipping():
            if clipped_bases[1] == 'S':
                count += clipped_bases[0]

        return count


def listholes(bamfilePath: os.path.realpath) -> typing.Set[int]:
    """Returns the set of all the HoleID contained in the bamfile (used in tests)"""

    bamfilePath = os.path.realpath(bamfilePath)
    logging.debug(
        '[DEBUG] Listing holes present in bamfile {}'.format(bamfilePath))
    setHoles = set()
    for line in read_bam(bamfilePath):
        holeID = int(line.split()[0].split('/')[1])
        setHoles.add(holeID)
    return setHoles
