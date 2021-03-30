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

<<<<<<< HEAD
    cmd = 'ipdSummary '+aligned_subreads+' --reference '+reference+' --pvalue 1 --identify m4C,m6A,m5C --csv '+csvoutput+' --log-level '+args["verbosity"]+' --identifyMinCov 3 --identify m6A,m4C --identifyMinCov 3 --gff '+gffoutput+" --ipdModel "+args["pathmodel"]
=======
    cmd = 'ipdSummary '+aligned_subreads+' --reference '+reference+' --pvalue 1 --identify m4C,m6A,m5C --csv '+csvoutput+' --log-level '+args["verbosity"]+' --identifyMinCov 3 --ipdModel /home/guillaume/conda3/envs/smsn/lib/python3.7/site-packages/kineticsTools/resources/SP2-C2.npz.gz'
>>>>>>> main
    logging.debug('[DEBUG] (simple_ipdSummary_stats) Launching the following cmd: {}'.format(cmd))
    logging.debug('[DEBUG] Launching ipdSummary (methylation analysis) --> {}'.format(cmd))
    call_process(cmd)

    df = read_csv(csvoutput)
    return df # This a pandas DataFrame