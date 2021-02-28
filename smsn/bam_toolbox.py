import numpy as np
import logging
import subprocess
import os
import pandas as pd

def read_bam(bamfilePath):
    """Reads a .bam line by line"""
    logging.debug(
        '[DEBUG] (read_bam) Launching samtools in subprocess to yield the lines of the bamfile {}'.format(bamfilePath))
    # Samtools must be installed !
    cmd = 'samtools view '+os.path.realpath(bamfilePath)
    logging.debug('[DEBUG] (read_bam) cmd = {}'.format(cmd))
    proc = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE)  # No .sam dropping

    list_txt_sam = ''  # Will contain a .sam string
    first = True

    logging.debug(
        '[DEBUG] (read_bam) Starting to read output from samtools for {}'.format(bamfilePath))

    while True:
        # We will read the output of samtools line by line. Decoding is needed because popen returns bytes, not strings
        line = proc.stdout.readline().decode('utf-8')
        if line:  # While there is something to read (else it will break)
            yield line
        else:  # There is nothing left in the stdout of samtools
            logging.debug(
                '[DEBUG] (read_bam) reading bamfile {} seems over'.format(bamfilePath))
            break

    return

def read_header(bamfilePath):
    """ returns the header as a python-string .sam"""

    samtxt = ''
    cmd = 'samtools view -h -S '+str(bamfilePath)
    bam = None
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                            stdin=subprocess.PIPE)  # We launch samtools view

    while True:
        output = proc.stdout.readline().decode('utf-8')
        if output[0] == '@':
            samtxt = samtxt + output
        else:
            break

    return samtxt

def parse_line(line,framerate=80):
    line = line.split('\t')

    returndict = {}
    returndict["scaffold"] = line[2]
    returndict["start"] = line[3]
    returndict["CIGAR_string"] = line[5]
    returndict["HoleID"] = int(line[0].split('/')[1])
    returndict["alignment_flag"] = line[1]
    returndict["sequence"] = line[9]
    returndict["ipd"] = np.array([np.int(x) for x in line[12]])

    return returndict



def get_pd_dataframe_alignment_custom(bamfile, framerate=80, genome = "None"):

    list_alignments = [] # We will return a list of dicts
    # In each dict, the key is the name of the feature while the value is its... value (Wow !)
    for line in read_bam(os.path.realpath(bamfile)): # Reads the .bam line after line
        if line:
            alignment_dict = parse_line(line,framerate)
            list_alignments.append(alignment_dict.copy())
    return pd.DataFrame(list_alignments) # We return a list of dict transformed in a pandas DataFrame

class CIGAR:

    def __init__(self, strCIGAR, scaffold=None):
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
        self.len = self.getlen()
        self.reflen = self._compute_reflen()
        self.identity = self._compute_identity()

    def __repr__(self):
        return self.string

    def __str__(self):
        return self.string

    def __eq__(self, other):
        if self.string == other.string and self.scaffold == other.scaffold:
            return True
        else:
            return False

    def __ne__(self, other):
        if self.string == other.string and other.scaffold == other.scaffold:
            return False
        else:
            return True

    def _compute_identity(self):
        # We will consider that clip != Identity
        return self.number_of_matches / (self.reflen + self.number_of_clipped + self.number_of_insertions)
        # Clipping will cause penalty...

    def _compute_reflen(self):
        binarlist = []
        for elt in self.listgroups:
            if elt[1] != 'I':
                for i in range(0, elt[0]):
                    binarlist.append(True)
            else:
                for i in range(0, elt[0]):
                    binarlist.append(False)
            del i

        # We have no idea if there are insertions in the clipped part
        reflen = np.sum(np.array(binarlist)) - self.number_of_clipped
        return reflen

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
        return np.array(self.binarlist, dtype=np.bool)

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

    def getlen(self):
        return len(self.binarlist)

def get_samflag(number):
    """According to https://broadinstitute.github.io/picard/explain-flags.html"""

    list_significations = ["read paired",
                           "read mapped in proper pair",
                           "read_unmapped",
                           "mate_unmapped",
                           "read reverse strand",
                           "mate reverse strand",
                           "first in pair",
                           "second in pair",
                           "not primary alignment",
                           "read failed platform/vendor quality checks",
                           "read is PCR or optical duplicate",
                           "supplementary alignment"
                           ]
    list_significations = list_significations[::-1]
    binarized = list(bin(number))[2:]

    while(len(binarized) < len(list_significations)):
        binarized.insert(0,"0")

    flags = []
    for i in binarized:
        if i == "1":
            flags.append(True)
        else:
            flags.append(False)

    returndict = {x:y for (x,y) in zip(list_significations,flags)}
    return returndict