#!/usr/bin/env python

__author__ = "DELEVOYE Guillaume"
__credits__ = ["BAHIN Mathieu", "MEYER Eric"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "DELEVOYE Guillaume"
__email__ = "delevoye@ens.fr"
__status__ = "Developpment"

import os
import logging
import sys
import subprocess

def parse_hacked_output(line):
    """The hack of ipdSummary works as follows: Printing on stdout the formatted results at the right moment.
    This function allows to parse the lines producted by the hacked ipdSummary to determine wether it is
    a stderr from ipdSummary or our hacked_results."""

    current_dict = {}

    if 'type_capping' not in line: # Trick to distinguish PacBio logging lines from actual outputs that we need
    # It is mandatory because workers that end send a message, something that we absolutely don't care about
        # We will still log the output in case it would be important...
        logging.debug('[DEBUG] (parse_hacked_output) Unusual line catch: {}'.format(line))
        current_dict['type_capping'] = 'ENDED' # This will serve as a flag for the calling function
        current_dict['message'] = line
        # to indicate that it shouldn't take care of it
        return current_dict


    # The following lines will deal about converting the printed output of hacked ipdSummary into their corresponding real python objects
    for elt in line.strip().split(';'): # ; separated every key:value
    #Like this:
        # key1:value1;key2:value2;key3:value3\n [etc]
        try:
            key = str(elt.split(':')[0]).strip() # This line and the following one should be easier to understand by reading the explanation just above
            value = str(elt.split(':')[1]).strip()
        except IndexError:
            logging.warning('[WARNING] (parse_hacked_output) The line contained type_capping but wasn\'t splitable: {}'.format(line))
        try:
            current_dict[str(key)] = int(value)# Try to get it as a python integer
        except:
            try:
                current_dict[str(key)] = float(value) # else a float...
            except:
                # if key == 'capped_values' or key == 'rawData': # And if you can't, then probably the output is actually the __repr__ of a python list
                #     value = value.replace('[', '').replace(']', '')
                #     value = "[" + ",".join(value.split()) + "]"
                #     # In the two lines above I'm just transforming the __repr__ list into a real python list
                #     current_dict[str(key)] = eval(value) # tadaa !
                # else:
                current_dict[str(key)] = str(value) # If it's not an integer nor a float nor a list then let's accept it as a str...

    return current_dict

def get_ipdSummary_details(aligned_subreads, reference, holeID, args):
    """Returns a list of dicts corresponding to the hacked output of ipdSummary"""

    HERE = os.getcwd()
    alignmentFile = os.path.realpath(aligned_subreads)
    reference = os.path.realpath(reference)

    logging.debug('[DEBUG] (get_ipdSummary_details) Pbindexing the bam in order to analyze it with ipdSummary')
    os.system('pbindex '+str(os.path.realpath(alignmentFile)))

    list_output = [] # A list of dict (one dict = one nucleotide on 1 strand)
    logging.debug('[DEBUG] (get_ipdSummary_details) Working on a single mode approach for hole {}'.format(str(holeID)))

    workdir = os.path.dirname(aligned_subreads)
    cmd = 'ipdSummary '+os.path.join(workdir,'aligned_on_restrictedscaffold_'+str(holeID)+'.bam')+' --reference '+reference+' --methylFraction --identify m4C,m6A,m5C --gff '+str(holeID)+'.gff'

    logging.debug('[DEBUG] (simple_ipdSummary_stats) Launching the following cmd: {}'.format(cmd))
    logging.info('[INFO] Launching ipdSummary (methylation analysis) --> {}'.format(cmd))
    proc = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = sys.stderr)

    for line in proc.stdout:
        output = parse_hacked_output(line.decode('utf-8'))
        if output:
            try:
                list_output.append(output)
            except KeyError:
                logging.warning('[WARNING] (get_ipdSummary_details) Skipping unexpected line: {}'.format(line))

    os.chdir(HERE)
    return list_output # This is a list of dicts !