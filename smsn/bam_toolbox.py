import numpy as np
import logging
import subprocess
import os
import pandas as pd
import smsn

def read_bam(bamfilePath):
    """Reads a .bam line by line"""
    logging.debug(
        '[DEBUG] (read_bam) Launching samtools in subprocess to yield the lines of the bamfile {}'.format(bamfilePath))
    # Samtools must be installed !
    cmd = 'samtools view '+os.path.realpath(bamfilePath)
    logging.debug('[DEBUG] (read_bam) cmd = {}'.format(cmd))
    proc = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE)  # No .sam dropping

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
    logging.debug("[DEBUG] Reading headers of {}".format(bamfilePath))

    samtxt = ''
    cmd = 'samtools view -h -S '+str(bamfilePath)
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                            stdin=subprocess.PIPE)  # We launch samtools view

    while True:
        output = proc.stdout.readline().decode('utf-8')
        if output[0] == '@':
            samtxt = samtxt + output
        else:
            break

    return samtxt

def get_hole_id(samseq):
    try:
        holeID = "".join([x for x in samseq.split('\n') if '@' not in x]).split()[0].split('/')[1]
    except:
        holeID = -1
    return holeID

def sort_by_name(bamfilepath, newpath, nbcore=1):
    """ Will sort the bam by read name (it will do it numerically for the numbers, not alphabetically, which is convenient for us)"""
    bampath = os.path.realpath(bamfilepath)
    bamdir = os.path.dirname(bampath)
    newpath = os.path.realpath(newpath)
    HERE = os.getcwd()
    os.system('mkdir -p '+str(bamdir))
    os.chdir(bamdir)
    cmd = 'samtools sort -n ' + \
        str(os.path.basename(bamfilepath))+' ' + \
        '-o '+str(newpath)+' -@ '+str(nbcore)
    os.system(cmd)
    os.chdir(HERE)


def inmemory_asbam(samstring):
    """ Returns a python byte-encoded .bam from a .sam string. Everything is done in memory, not on hard drive """

    p = subprocess.run('samtools view -bS -h', shell=True,
                       stdout=subprocess.PIPE, input=samstring.encode('utf-8'))
    return(p.stdout)


def yield_hole_by_hole(subreads_bamfile, header=None, restricted_to=None, min_subreads=0):
    """ Yields the subreads HoleID by HoleID as .bam
    subreads_bamfile must be already sorted by HoleID.
    "restricted_to" can be the exhaustive set of HoleID (int) to include.

    Can also exclude the holes with a critically low number of subreads

    Yields a tuple (HoleID,ignored,reason,samtxt)

    """

    def yield_hole(samtxt, header, restricted_to, min_subreads):  # subfunction
        nbsubreads = len(samtxt)
        holeID = int(samtxt[0].split("/")[1])

        samtxt = "\n".join([x.strip() for x in samtxt])

        if header:
            samtxt = header + samtxt


        logging.debug("[DEBUG] (yield_hole) {} subreads for HoleID {}".format(nbsubreads, holeID))

        if restricted_to:
            logging.debug("[DEBUG] (yield_hole) There's a restriction list")
            isin_shortlist = holeID in restricted_to

            if not isin_shortlist:
                logging.debug("[DEBUG] (yield_hole) holeID {} is not in the restriction list".format(holeID))
                logging.debug("[DEBUG] (yield_hole) HoleID {} Not in shortlist: Skipped".format(holeID))
                return (holeID, True, "Not in shortlist", None)
            else:
                logging.debug('[DEBUG] (yield_hole) holeID {} is in the shortlist'.format(isin_shortlist))

        # In case it is in the shortlist
        if nbsubreads < min_subreads: # Check if thats ok for the subreads
            logging.debug(
                '[DEBUG] (yield_hole) HoleID {} not enough subreads ({} VS {} min): Skipped'.format(holeID,
                                                                                                    nbsubreads,
                                                                                                    min_subreads))
            return (holeID, True, "Not enough subreads", None)
        else:
            logging.debug("[DEBUG] (yield_hole) HoleID {} has sufficient subreads ({} VS {})".format(holeID,
                                                                                                     nbsubreads,
                                                                                                     min_subreads))

        logging.debug('[DEBUG] (yield_hole) HoleID {} yielded successfully'.format(holeID))
        return (holeID, False, None, samtxt)

    #############################

    if restricted_to:
        logging.debug("[DEBUG] (yield hole by hole) Size of restriction list = {}".format(len(restricted_to)))
        logging.debug('[DEBUG] (yield hole by hole). List restriction = {}'.format(restricted_to))

    logging.debug('[DEBUG] (yield by hole). min_subreads = {}'.format(min_subreads))

    reader = smsn.bam_toolbox.read_bam(subreads_bamfile)
    setEncountered = set()
    current_yield = []

    first = True

    for current_line in reader:
        try:
            current_holeID = int(current_line.split()[0].split('/')[1])
        except Exception as e:
            logging.error("[ERROR] (yield hole by hole) Error encountered with the following exception :\n{}".format(e))
            logging.error("[ERROR] (yield hole by hole) Exiting -1")
            exit(-1)

        if current_holeID not in setEncountered:
            logging.debug("[DEBUG] Never encountered {}".format(current_holeID))
            if first == False:
                yield yield_hole(samtxt=current_yield, header=header,
                                 restricted_to=restricted_to, min_subreads=min_subreads)
                current_yield = []
                setEncountered.add(current_holeID)
                current_yield.append(current_line)
            else:
                first = False
                logging.debug('[DEBUG] Not first anymore')
                setEncountered.add(current_holeID)
                current_yield.append(current_line)
        else:
            current_yield.append(current_line)

    # Don't forget the last line
    if current_yield:
        logging.debug("[DEBUG] lastline")
        logging.debug('[DEBUG] yielding {}'.format(current_holeID))
        yield yield_hole(samtxt=current_yield, header=header,
                         restricted_to=restricted_to, min_subreads=min_subreads)

    logging.debug('[DEBUG] Yielding bamfile hole by hole is over')
    return

def parse_header(bamfile):
    """Returns as a dict the crucial informations that interest us in the header (ID, PL, PU, PM, READTYPE, Ipd:CodecV1, FRAMERATEHZ, etc"""

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


def parse_line(line,include_ipds=False):
    line = line.split('\t')

    returndict = {}
    returndict["alignment_flag"] = int(line[1])
    returndict["scaffold"] = line[2]
    returndict["start"] = line[3]
    returndict["CIGAR_string"] = line[5] # A CIGAR object will be created later to avoid saturating the RAM
    returndict["HoleID"] = int(line[0].split('/')[1])
    returndict["sequence"] = line[9]
    returndict["moviename"] = line[0].split('/')[0]

    if include_ipds:
        returndict["encoded_ipd"] = [x for x in line if x[:3]=="ip:"][0]
        returndict["encoded_ipd"] = np.array(returndict["encoded_ipd"].split(',')[1:],dtype=np.uint8) # This should then be decoded, but using np.uint8 ensures that not too much place is taken in RAM

    return returndict


def detailed_get_pd_dataframe_alignment(bamfile, include_ipds=False):
    """ From a .bam alignment file, gives a pandas dataframe corresponding to the alignment"""

    list_alignments = []  # We will return a list of dicts
    # In each dict, the key is the name of the feature while the value is its... value (Wow !)
    for line in read_bam(os.path.realpath(bamfile)):  # Reads the .bam line after line
        if line:
            line = line.split('\t')
            moviename = line[0].split('/')[0]
            scaffold = line[2]
            start = int(line[3])
            CIGAR_string = line[5]
            try:
                CIGAR_obj = CIGAR(
                    CIGAR_string)  # CIGAR is a very usefull object I've developped to parse CIGAR strings...
            # [...] and get many informations about them
            except Exception as e:
                logging.error("[ERROR] A problem occured on line : \n{} \n\n\t***** Problem = CIGAR {}".format(line, CIGAR_string))
                logging.error('[ERROR] Exception recieved = {}'.format(e))

            identity_score = CIGAR_obj.identity
            end = start + CIGAR_obj.reflen  # Reflen is the size of reference deduced from the CIGAR string
            # Because the alignment starting position doesn't take the soft-clipping in account,
            # We don't take it in account either (reflen doesn't take clipped bases in account)

            # But clipped bases are counted as unaligned for %identity computation !

            reflen = end - start

            HoleID = int(line[0].split('/')[1])
            alignment_flag = int(line[1])

            mapQV = int(line[4])
            clipped_bases = CIGAR_obj.number_of_clipped
            matching_bases = CIGAR_obj.number_of_matches

            sequence = line[9]

            list_interest = ["moviename", "scaffold", "start", "CIGAR_string", "identity_score", "end", "HoleID",
             "alignment_flag", "reflen", "mapQV", "matching_bases", "clipped_bases", "sequence"]

            if include_ipds:
                encoded_ipds = [x for x in line if x[:3] == "ip:"][0]
                encoded_ipds = np.array(encoded_ipds.split(',')[1:],
                                                     dtype=np.uint8)  # This should then be decoded, but using np.uint8 ensures that not too much place is taken in RAM
                list_interest.append("encoded_ipds")

            alignment_dict = {}
            for elt in list_interest:
                # We add all these values to a dict
                alignment_dict[elt] = eval(elt)

            list_alignments.append(alignment_dict)
            del line
    return pd.DataFrame(list_alignments)  # We return a list of dict transformed in a pandas DataFrame

def listholes(bamfilePath):
    """Returns the set of all the HoleID contained in the bamfile"""

    bamfilePath = os.path.realpath(bamfilePath)
    logging.debug('[DEBUG] Listing holes present in bamfile {}'.format(bamfilePath))
    setHoles = set()
    for line in read_bam(bamfilePath):
        holeID = int(line.split()[0].split-('/')[1])
        setHoles.add(holeID)
    return setHoles

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
