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
from pandas import read_csv
from smsn.pipeline import call_process


def launch_ipdSummary(aligned_subreads, reference, holeID, args):
    """Returns a list of dicts corresponding to the hacked output of ipdSummary"""

    # reference = os.path.realpath(reference)
    # alignmentFile = os.path.realpath(aligned_subreads)
    # Déjà fait ?
    # logging.debug('[DEBUG] (get_ipdSummary_details) Pbindexing the bam in order to analyze it with ipdSummary')
    #
    # cmd = 'pbindex '+str(os.path.realpath(alignmentFile))
    # call_process(cmd)
    # logging.debug('[DEBUG] (get_ipdSummary_details) Working on a single mode approach for hole {}'.format(str(holeID)))

    workdir = os.path.dirname(aligned_subreads)
    csvoutput = os.path.join(os.path.join(workdir,str(holeID)),str(holeID)+".csv")

    cmd = 'ipdSummary '+aligned_subreads+' --reference '+reference+' --pvalue 1 --identify m4C,m6A,m5C --csv '+csvoutput+' --log-level '+args["verbosity"]+' --identifyMinCov 3'
    logging.debug('[DEBUG] (simple_ipdSummary_stats) Launching the following cmd: {}'.format(cmd))
    logging.debug('[DEBUG] Launching ipdSummary (methylation analysis) --> {}'.format(cmd))
    call_process(cmd)

    df = read_csv(csvoutput)
    return df # This a pandas DataFrame