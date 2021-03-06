#!/usr/bin/env python

__author__ = "DELEVOYE Guillaume"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "DELEVOYE Guillaume"
__email__ = "delevoye@ens.fr"
__status__ = "Developpment"

import subprocess
import os
import smsn
import logging
import pandas as pd
from pandarallel import pandarallel
import gc

from smsn.extract_model import transform_model_name

def call_process(cmd):
    processname = cmd.split()[0]
    logging.debug('[DEBUG] (call_process) Processname = {}'.format(processname))
    logging.debug('[DEBUG] (call_process) cmd = {}'.format(cmd))

    process = subprocess.run(cmd, shell=True,capture_output=True)

    std_output = process.stdout.decode('utf-8')
    error_output = process.stderr.decode('utf-8')

    if std_output.strip():
        if "ERROR" in std_output.strip():
            logging.error('[ERROR] (stdout of process {} yielded an error: {}'.format(processname,std_output.strip()))
        else:
            logging.debug("[DEBUG] (stdout of process {} : {}\n".format(processname,std_output.strip().replace('\n','\\n')))

    if error_output.strip():
        if "ERROR" not in error_output and "error" not in error_output: # In order to not be spammed by CCS which has a bugged output
            logging.debug('[DEBUG] (stderr of process {} : {}\n'.format(processname,error_output.strip().replace('\n','\\n')))
        else:
            logging.error("[ERROR] (stderr output of process {}) : {}\n".format(processname,error_output.strip()))
# This will show you eventual errors + will force the kernel to wait the end of the process

def handle_spaces_reference(referencepath,tmpdir):
    """A preprocessing is needed because none of the PacBio tools (BLASR, CCS, kineticsTools) handles
    whitespaces or special characters in theidentifiers of the scaffolds in the .fasta reference file.

    In case the user's reference contains whitespaces in its identifiers, we'll automatically switch them by "_"
    in another .fasta file located in the tmpdir. Since this can break future analysis, user will be properly
    warned, and asked to provide whitespace-free references instead."""

    loaded_fasta = smsn.fasta_tools.load_fasta(referencepath)

    contains_some = False
    list_specials = [" ","#",".",",",",",",","*","-","&","|","[","]","{","}","\"","\'","~","%","+",">","<","!","?","/","\\","??","??",":"]
    for key in list(loaded_fasta):
        if any([character in list_specials for character in key]):
            contains_some = True
        else:
            pass

    if contains_some:
        logging.warning('[WARNING] At least one ID in the .fasta file contains whitespaces or unsupported characters.')
        referencename = os.path.basename(referencepath)
        newfasta = os.path.join(tmpdir,"whitespacefree_"+referencename)
        loaded_fasta = smsn.fasta_tools.load_fasta_special(referencepath)
        smsn.fasta_tools.dict_to_fasta(loaded_fasta,newfasta)
        logging.warning('[WARNING] Since the identifiers in the .fasta reference contain special characters or whitespaces, we built another reference at {}, where identifiers are whitespace-free without special characters. This might force you to make additionnal preprocessing before analyzing your results, since the scaffold names are now changed. Please rename your scaffolds to avoid future pitfalls.'.format(newfasta))
        return newfasta
    else:
        return referencepath

def launch_smsn(args):
    """Launchs the whole pipeline. That's where every functions are called in order.
    Args are parsed from user input through the smsn_launcher.py script"""
    moviename = smsn.bam_toolbox.parse_header(args["bam"])["PU"]
    bam_header = smsn.bam_toolbox.read_header(args["bam"])

    args["moviename"] = moviename
    args["bam_header"] = bam_header

    if args["model"].upper() != "AUTO":
        logging.warning('[WARNING] You asked to use model {}, but maybe it will not be appropriate. Please be sure of yourself'.format(args["model"]))
        args["pathmodel"] = transform_model_name(args["model"])
        logging.info('[INFO] Using model = {}'.format(args["pathmodel"]))

    logging.debug("[DEBUG] Ensuring that the .fasta doesn't contained empty spaces in the references")
    args["reference"] = handle_spaces_reference(args["reference"],args["tmpdir"]) # If there is nothing to change, the reference will stay the same
    # Otherwise, it will be changed to another ref that will be in the tmpdir

    ####################################################################################################################
    logging.info("[INFO] Step 1 : Creating the consensus")
    if args["CCS"]:
        # This step is really expensive computationnally. We'll let the user give its own consensus if they're already
        # built
        args["CCS"] = os.path.realpath(args["CCS"])
        logging.info("[INFO] CCS {} will be used - CCS creation is skipped.".format(args["CCS"]))
    else:
        recreate_CCS(args["bam"],args["nb_proc"],tmpdir=args["tmpdir"])
        args["CCS"] = os.path.join(args["tmpdir"],moviename+".CCS.bam")

    ####################################################################################################################
    logging.info("[INFO] Step 2 - 1: Aligning all the CCS against the reference using BLASR")

    args["aligned_CCS"] = os.path.join(args["tmpdir"],"alignedCCS_"+os.path.basename(args["CCS"]))
    args["notaligned_CCS"] = os.path.join(args["tmpdir"],"NOTalignedCCS_"+moviename+".fasta")

    align_CCS(args["CCS"],
              args["reference"],
              args["aligned_CCS"],
              args["notaligned_CCS"],
              args["tmpdir"],
              args["nb_proc"])

    logging.debug("[DEBUG] We make the list of all holes to be analyzed")

    ####################################################################################################################
    logging.info('[INFO] Step 2 - 2: Loading the alignment produced')
    df_aligned_CCS = smsn.bam_toolbox.detailed_get_pd_dataframe_alignment(args["aligned_CCS"],include_ipds=False)

    try: # Safety check
        assert len(df_aligned_CCS) == len(df_aligned_CCS['HoleID'].unique()) # Just checking that BLASR did everything right
    except AssertionError:
        logging.critical('[CRITICAL] There is not only one alignment per HoleID. The rest of the analysis might be '
                         'corrupted')

    output_alignmentcsv = os.path.join(args["tmpdir"],"aligned_CCS_"+args["moviename"]+".csv")
    logging.debug("[DEBUG] Saving alignment (.csv) at {}".format(output_alignmentcsv))
    df_aligned_CCS.to_csv(output_alignmentcsv,sep=";",index=False)

    ####################################################################################################################
    logging.info('[INFO] Step 2 - 3: Keeping only alignments with >= {} identity'.format(args["min_identity"]))
    df_aligned_CCS = df_aligned_CCS[df_aligned_CCS["identity_score"] >= args["min_identity"]]
    df_aligned_CCS.index = [int(x) for x in df_aligned_CCS["HoleID"]]

    set_holes_to_analyse = set(df_aligned_CCS["HoleID"].unique())

    logging.debug('[DEBUG] {} holes are kept after filtering upon identity (Before filtering we had {})'.format(len(set_holes_to_analyse),len(df_aligned_CCS)))


    ####################################################################################################################
    logging.debug("[DEBUG] Initializing pandarallel")
    if args["progress_bar"] == True: # Initialize pandarallel
        logging.debug('[DEBUG] Pandarallel will use a progress_bar')
        pandarallel.initialize(progress_bar=True,nb_workers=args["nb_proc"])
    else:
        logging.debug('[DEBUG] Pandarallel will not use a progress_bar')
        pandarallel.initialize(progress_bar=False,nb_workers=args["nb_proc"])

    ####################################################################################################################
    logging.debug("[DEBUG] Entering chunked data + Parallelism")
    logging.info("[INFO] Starting analysis")
    reader = smsn.bam_toolbox.yield_hole_by_hole(args["bam"],
                                                 restricted_to=set_holes_to_analyse,
                                                 min_subreads=args["min_subreads"])


    nb_analyzed = 0
    proto_df = []
    i = 0
    chunknb = 0
    for HoleID,ignored,reason,samtxt in reader:

            if not ignored:
                alignment_this_hole = df_aligned_CCS.loc[HoleID]
                newdict = {}
                newdict["HoleID"] = HoleID
                newdict["sam"] = args["bam_header"] + samtxt
                newdict["scaffold"] = alignment_this_hole["scaffold"]
                newdict["start"] = alignment_this_hole["start"]
                newdict["end"] = alignment_this_hole["end"]

                logging.debug('[DEBUG] (pipeline) New hole added to analysis --> HoleID {}, scaffold {}, start {}, end {}, subreads = {}'.format(newdict["HoleID"],newdict["scaffold"],newdict['start'],newdict['end'],len(samtxt.split("\n"))))
                proto_df.append(newdict.copy())

                i+=1

                if i >= args["sizechunks"]:
                    logging.info('[INFO] Analyzing chunk n?? {} of size {}'.format(chunknb,i))

                    ### Indicate the mapping positions:
                    df = pd.DataFrame(proto_df)
                    outputs = df.parallel_apply(lambda x: smsn.single_hole.analyze_singleHole(x["HoleID"],
                                                                             x["sam"],
                                                                             x["scaffold"],
                                                                             x["start"],
                                                                             x["end"],
                                                                             args),axis=1) # Simple trick to make it parallel
                    nb_analyzed += len(outputs)
                    outputs = pd.concat([x for x in outputs],ignore_index=True)

                    chunk_csvpath = os.path.join(args["tmpdir"],"tmp_analysis_chunk_"+str(chunknb)+".csv")
                    logging.debug('[DEBUG] Compiling all results for chunk n?? {} into {}'.format(chunknb,chunk_csvpath))
                    outputs.to_csv(chunk_csvpath,sep=";")

                    chunknb += 1
                    i=0

                    del df
                    del outputs
                    gc.collect()

                    proto_df = []


    # Handling the last chunk (happens if its size is < sizechunks)
    if len(proto_df)>0:
        logging.info('[INFO] Analyzing last chunk - n?? {} of size {}'.format(chunknb,i))

        ### Indicate the mapping positions:
        df = pd.DataFrame(proto_df)
        outputs = df.parallel_apply(lambda x: smsn.single_hole.analyze_singleHole(x["HoleID"],
                                                                                  x["sam"],
                                                                                  x["scaffold"],
                                                                                  x["start"],
                                                                                  x["end"],
                                                                                  args),
                                    axis=1)  # Simple trick to make it parallel

        nb_analyzed += len(outputs)
        chunk_csvpath = os.path.join(args["tmpdir"], "tmp_analysis_chunk_" + str(chunknb) + ".csv")
        logging.info('[DEBUG] Compiling all results for chunk n?? {} into {}'.format(chunknb, chunk_csvpath))
        outputs = pd.concat([x for x in outputs], ignore_index=True)
        outputs.to_csv(chunk_csvpath, sep=";")

        chunknb += 1
        i = 0

        del df
        del outputs
        gc.collect()

    ####################################################################################################################
    logging.info("[INFO] Compiling all results of all chunks into a single file")

    if nb_analyzed > 0: # Sometimes, no molecule can be analyzed and it's cleaner to leave like this
        # than to leave with a python exception traceback
        pass
    else:
        logging.error('[ERROR] Not a single molecule could be analyzed given the following criteria : {}'.format(args))
        logging.info('[INFO] SMSN ended.')
        exit(-1)

    outputs = [elt for elt in os.listdir(args["tmpdir"]) if 'tmp_analysis_chunk' in elt]

    if len(outputs) < 1:
        logging.error("[ERROR] No output .csv file")

    try: # If compilation fails for some reason, be certain that we don't destroy the temporary files
        # Computation time is expensive and boring for everyone, and also it kills polar bears
        # So, give a chance to the user to get the best out of what has already been computed
        compilation = []
        for elt in os.listdir(args["tmpdir"]):
            if "tmp_analysis_chunk" in elt:
                tmp_path = os.path.join(args["tmpdir"],elt)
                compilation.append(pd.read_csv(tmp_path,sep=";"))

        list_interest = ["HoleID","scaffold","tpl","strand","base","score","tMean","tErr","modelPrediction","ipdRatio","coverage","isboundary"] # identificationQv, context
        # list_interest = list(compilation[0])
        if args["idQvs"]:
            list_interest.append("identificationQv")

        if args["add_context"]:
            list_interest.append("context")

        output_df = pd.concat(compilation,ignore_index=True).drop_duplicates().reset_index(drop=True)
        output_df.sort_values(["HoleID","tpl","strand"],inplace=True)
        output_df[list_interest].to_csv(args["output_csv"],sep=";",index=False)
        logging.info('[INFO] Result saved at {}'.format(args["output_csv"]))
        failed = False

    except Exception as e: # Warn when operation failed (ex: I/O error, server down, whatever)
        logging.error('[ERROR] Compilation of all the output csv files failed with the following error: \n{}'.format(e))
        logging.error("[ERROR] Problem with compiling into {}".format(args["output_csv"]))
        failed = True

    if not failed: # Clean the temporary files only if the compilation suceeded
        for filetmp in os.listdir(args["tmpdir"]):
            if "tmp_analysis_chunk" in filetmp:
                os.remove(os.path.join(args["tmpdir"],filetmp))
    else:
        logging.info('[INFO] For some reason, the compilation was not possible. However, it should still possible to '
                     'retrieve the tmp .csv files at {}'.format(args["tmpdir"]))


    logging.info('[INFO] SMSN ended.')



#### STEP 1

def recreate_CCS(subread_file, nbproc,tmpdir):
    """Recreates the circular consensus (CCS) from a .subreads.bam (will be named [moviename].CCS.bam, in the same
    directory as the subreads. Very laxist parameters will be used in order to yield as many CCS as possible. Because we
    are not interested here in doing an assembly, it is the mapping will then determine whether a given CCS is a good
    one or not. Some files produced by PacBio will be interminently found outside of tmpdir"""

    logging.info(
        '[INFO] Recreating the circular consensus from subreads : {}, using {} CPU.'.format(subread_file, nbproc))

    os.system("mkdir -p "+os.path.realpath(tmpdir))

    subread_file = os.path.realpath(subread_file)
    subread_dir = os.path.dirname(subread_file)
    moviename = smsn.bam_toolbox.parse_header(subread_file)["PU"]
    outputbam = os.path.join(os.path.realpath(tmpdir), moviename + ".CCS.bam")

    # Here I use very laxist parameters so that I'm sure that almost every hole is going to give a consensus
    cmd = "ccs --minLength 0 --minIdentity 0 --minZScore NaN --minSnr 1 --polish --minReadScore 0 --richQVs --maxLength 5000000 --minPasses 0 --minPredictedAccuracy 0 --numThreads " + str(
        nbproc) + " " + subread_file + " " + outputbam
    logging.debug("[DEBUG] Launching cmd = {}".format(cmd))
    call_process(cmd)

    filin = open('ccs_report.txt', "r")
    output_report = filin.read()
    filin.close()

    logging.info("[INFO] Report generated = {}".format(output_report))

    logging.debug("[DEBUG] Renaming the report")

    filout = open(moviename + "_ccs_report.txt", "w")
    filout.write(output_report)
    filout.close()

    os.remove("ccs_report.txt")

    logging.debug("[DEBUG] Moving the CCS report into the .tmp directory")
    reportpath = os.path.join(subread_dir, moviename + "_ccs_report.txt")
    reporttmppath = os.path.join(tmpdir,moviename+"_ccs_report.txt")
    os.rename(reportpath,reporttmppath)

#### STEP 2

def align_CCS(CCS,reference,aligned_CCS,notaligned_CCS,tmpdir,nb_proc):
    cmd = "blasr " + CCS + " " + reference + " --bestn 1 --hitPolicy leftmost --bam --out " + aligned_CCS + " --clipping none --nproc " + str(nb_proc) + " --unaligned " + os.path.join(tmpdir,notaligned_CCS)
    call_process(cmd)
    logging.debug('[DEBUG] Mapping DONE for {} on {}. Output canbe found at {}'.format(CCS,reference,aligned_CCS))