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
import numpy as np
import pandas as pd
from pandas import read_csv
from smsn.pipeline import call_process
from smsn.fasta_tools import get_snipet, parse_gff, load_fasta


def launch_ipdSummary(aligned_subreads, reference, holeID,
                      offset, real_start, real_end, args):
    logging.debug(
        '[DEBUG] (launch_ipdSummary) recieved aligned_subreads = {} and reference = {} and holeID = {}'.format(
            aligned_subreads,
            reference,
            holeID))

    workdir = os.path.dirname(aligned_subreads)
    csvoutput = os.path.join(workdir, str(holeID) + ".csv")
    gffoutput = os.path.join(workdir, 'output.gff')

    logging.debug(
        '[DEBUG] (launch_ipdSummary) workdir = {}, csvoutput = {}'.format(
            workdir, csvoutput))
    cmd = 'ipdSummary ' + aligned_subreads + ' --reference ' + reference + ' --pvalue 1 --identify m4C,m6A --csv ' + \
        csvoutput + ' --log-level DEBUG ' + \
        ' --minCoverage 0 --identifyMinCov 0 --identifyMinCov 0 --gff ' + gffoutput
    cmd = cmd + " --ipdModel " + args["pathmodel"]

    logging.debug(
        '[DEBUG] (simple_ipdSummary_stats) Launching the following cmd: {}'.format(cmd))
    logging.debug(
        '[DEBUG] Launching ipdSummary (methylation analysis) --> {}'.format(cmd))

    stdout, stderr = call_process(cmd)  # Launch ipdSummary

    # Now let's get the results...

    list_interest = ["HoleID", "scaffold", "tpl", "strand", "base", "score", "tMean", "tErr", "modelPrediction",
                     "ipdRatio", "coverage", "isboundary"]  # identificationQv, context
    if args["idQvs"]:
        list_interest.append("identificationQv")
    if args["add_context"]:
        list_interest.append("context")

    error = False

    try:
        df = read_csv(csvoutput)
    except:
        logging.error('[ERROR] Exception while reading output of holeID {}'.format(holeID))
        error = True

    gff = parse_gff('output.gff')

    if len(gff) > 0 and "identificationQv" in gff.columns:
        gff["strand"] = [0 if x == "+" else 1 for x in gff["strand"]]

        gff["uniqueID"] = [str(x) + "_" + str(y) for (x, y) in zip(
            gff["start"],
            gff["strand"])]

        df["uniqueID"] = [str(x) + "_" + str(y) for (x, y) in zip(
            df["tpl"],
            df["strand"])]

        df = pd.merge(left=df,
                      right=gff[["uniqueID",
                                 "identificationQv"]],
                      how="outer",
                      on=["uniqueID"])

            # This operation is necessarly a bit complex since idQv can be computed only for adenines and cytosines
        # This being said, it's not impossible to understand it: The idQvs are outputed in the .gff only
        # So we get them from the .gff, and we put them back into the methylationout.csv file, using the fact
        # that only 1 line is outputed in both file for each physical
        # nucleotide that was sequenced
    else:
        df["identificationQv"] = [np.nan for x in range(len(df))]

    if not error:

        # Returning the results into the proper coordinates # Warning:
        # KineticsTools is already 1-indexed
        df['tpl'] = df['tpl'] + offset
        df["isboundary"] = [min(x - real_start, real_end - x)
                            < 10 for x in df["tpl"].astype(int)]

        df["HoleID"] = [holeID] * len(df)
        df["scaffold"] = df["refName"]
        del df["refName"]
        if args["add_context"]:
            fasta_chunked_ref = load_fasta(reference)
            df["context"] = df.apply(
                lambda x: get_snipet(seq=list(fasta_chunked_ref.values())[0],
                                     position=x["tpl"] - 2 - offset,
                                     strand=x["strand"],
                                     n=12),
                axis=1)

    else: # An error has occured
        logging.warning(
            '[WARNING] (launch_ipdSummary) Failed to get any data from ipdSummary on hole {}. This might be caused by '
            'a coverage that is too low on this holeID.'.format(holeID))
        artificial_tmpdict = {}  # So that we still return something...
        for column in list_interest:
            if column == "HoleID":
                artificial_tmpdict["HoleID"] = int(holeID)
            artificial_tmpdict[column] = np.nan

        df = pd.DataFrame([artificial_tmpdict])

    return df[list_interest]


