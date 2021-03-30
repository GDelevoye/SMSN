#!/usr/bin/env python

__author__ = "DELEVOYE Guillaume"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "DELEVOYE Guillaume"
__email__ = "delevoye@ens.fr"
__status__ = "Developpment"

import os
import logging
from smsn.bam_toolbox import inmemory_asbam
# import pandas as pd
from smsn.summary_details import launch_ipdSummary
import copy
import shutil
from smsn.pipeline import call_process
from smsn.fasta_tools import load_fasta_special, compute_chunk_infos, dict_to_fasta, get_snipet, load_fasta

def analyze_singleHole(holeID,samseq,scaffold,real_start,real_end,args):

    fastafile = args["reference"]
    workdir = args["tmpdir"]

    fasta = load_fasta_special(os.path.realpath(fastafile))
    HERE = os.getcwd()
    os.chdir(workdir)

    os.system('mkdir -p '+str(holeID))
    os.chdir(str(holeID))

    filout = open(str(holeID)+'.bam', 'wb') # Let's drop the .bam file here
    filout.write(inmemory_asbam(samseq)) # We need a clean unaligned .bam
    filout.close()

    aligned_bam_path = os.path.join(os.getcwd(),"aligned_on_restrictedscaffold_"+str(holeID)+".bam")
    chunked_ref_path = os.path.join(os.getcwd(),"chunked_ref.fasta")
    unaligned_bam_path = os.path.join(os.getcwd(),str(holeID)+".bam")

    # Indexing the unaligned .bam
    logging.debug('[DEBUG] (analyze_singleHole) Generating index for the unaligned .bam of holeID {}'.format(holeID))
    cmd = 'pbindex '+unaligned_bam_path
    call_process(cmd)
    # cmd = 'samtools index '+unaligned_bam_path
    # call_process(cmd)

    (chunk_start, chunk_end, chunk_size, offset, sequence) = compute_chunk_infos(int(real_start), int(real_end), str(fasta[scaffold]))

    # Taking only a sub-reference (+ 100 nt / -100 nt) to re-align all the subreads precisely where the CCS mapped

    dict_to_fasta({ str(scaffold) : str(sequence) }, './chunked_ref.fasta')
    logging.debug('[DEBUG] (analyze_singleHole) Indexing the chunk_ref.fasta of hole {}'.format(holeID))
    cmd ='samtools faidx '+chunked_ref_path
    call_process(cmd)

    # Map all the subreads restrictively on the scaffold
    cmd = 'blasr ' + unaligned_bam_path +" "+chunked_ref_path +' --useccs --bestn 1 --clipping none --bam --out '+aligned_bam_path+' --unaligned '+ os.path.join(os.getcwd(),str(holeID)+'.unaligned.fasta')
    call_process(cmd)

    # Indexing the mapped .bam
    logging.debug('[DEBUG] (analyze_singleHole) Generating index for the aligned .bam on restricted scaffold for hole {}'.format(holeID))
    cmd = 'pbindex '+aligned_bam_path
    call_process(cmd)
    # cmd = 'samtools index '+aligned_bam_path
    # call_process(cmd)

    logging.debug('[DEBUG] Just before calling ipdSummary, content of directory {} is {}'.format(os.getcwd(),os.listdir(os.getcwd())))
    logging.debug("[DEBUG] aligned_bam_path = {} and chunked_ref_path = {}".format(aligned_bam_path,chunked_ref_path))

    # Perform the analysis itself
    results = launch_ipdSummary(aligned_bam_path,
                                     chunked_ref_path,
                                     holeID = holeID,
                                     args = args) # WE RECIEVE A PD.DATAFRAME

    try:
        results['tpl'] = results['tpl'] + 1 + offset# Returning the results into the proper coordinates
    except:
        results['tpl'] = 'NaN'

    results["HoleID"] = int(holeID)
    results["scaffold"] = str(scaffold)

    # Markine the nucleotides that are less than 10nt away from the extremity
    results["isboundary"] = [min(x-real_start,real_end-x) <= 10 for x in results["tpl"]]

    try:
        if args["add_context"]:
            fasta_chunked_ref = load_fasta(chunked_ref_path)
            results["context"] = results.apply(lambda x: get_snipet(seq=list(fasta_chunked_ref.values())[0],
                                             position=x["tpl"]-2-offset,
                                             strand=x["strand"],
                                             n=12),axis=1)
    except KeyError:
        logging.debug("[DEBUG] Adding the contexts is skipped")

    os.chdir(HERE)

    path_thishole_tmpdir = os.path.join(args["tmpdir"],str(holeID))
    logging.debug('[DEBUG] Deleting {}'.format(path_thishole_tmpdir))

    if not args["preserve_tmpdir"]:
        shutil.rmtree(path_thishole_tmpdir,ignore_errors=True)

    return results



def clean_sam(samseq):
    """ Some .bam have already been aligned once, which might cause problems when trying to re-align them.
    To avoid that, this function will remove any existing alignment and returning an non-aligned .sam sequence """
    outlines = []
    for line in samseq.split('\n'):
        if line[0] == '@':
            outlines.append(line)
        else:
            splitted = line.split()

            id = str(splitted[0])
            other = '\t'.join(splitted[9:])
            newline = '\t'.join([id, '4', '*', '0', '255', '*', '*', '0', '0', other])
            outlines.append(copy.deepcopy(newline))
            del newline

    return ('\n'.join(outlines))

