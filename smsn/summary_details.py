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
    logging.debug('[DEBUG] (launch_ipdSummary) recieved aligned_subreads = {} and reference = {} and holeID = {}'.format(aligned_subreads,reference,holeID))

    workdir = os.path.dirname(aligned_subreads)
    csvoutput = os.path.join(workdir,str(holeID)+".csv")
    gffoutput = os.path.join(workdir,'output.gff')

    logging.debug('[DEBUG] (launch_ipdSummary) workdir = {}, csvoutput = {}'.format(workdir,csvoutput))
    cmd = 'ipdSummary '+aligned_subreads+' --reference '+reference+' --pvalue 1 --identify m4C,m6A --csv '+csvoutput+' --log-level DEBUG '+' --minCoverage 0 --identifyMinCov 0 --identifyMinCov 0 --gff '+gffoutput

    if args["model"].upper() != "AUTO":
        cmd = cmd + " --ipdModel "+args["pathmodel"]

    logging.debug('[DEBUG] (simple_ipdSummary_stats) Launching the following cmd: {}'.format(cmd))
    logging.debug('[DEBUG] Launching ipdSummary (methylation analysis) --> {}'.format(cmd))
    call_process(cmd)

    df = read_csv(csvoutput)

    # Handle cases where it's not computed because.. Well... kineticsTools refuses to output it (mostly because of coverage)
    if len(df) < 1:
        logging.warning('[WARNING] (launch_ipdSummary) Failed to get any data from ipdSummary on hole {}. This might be caused by a coverage that is too low on this holeID.'.format(holeID))
        artificial_tmpdict = {}
        for column in list(df):
            artificial_tmpdict[column] = "ERROR: CSV WAS EMPTY"
        artificial_tmpdict["HoleID"] = holeID
        artificial_tmpdict["isboundary"] = "ERROR: CSV WAS EMPTY"
        artificial_tmpdict["context"] = "ERROR: CSV WAS EMPTY"
        artificial_tmpdict["idQvs"] = "ERROR: CSV WAS EMPTY"
        artificial_tmpdict["tpl"] = "ERROR: CSV WAS EMPTY"
        logging.debug('[DEBUG] (launch_ipdSummary) Retunring an error for holeID {}'.format(holeID))
        df = df.append(artificial_tmpdict,ignore_index=True)


    return df # This a pandas DataFrame