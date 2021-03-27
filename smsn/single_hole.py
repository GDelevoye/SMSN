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
from smsn.bam_toolbox import inmemory_asbam
# import pandas as pd
from smsn.summary_details import launch_ipdSummary
import copy
import shutil
from smsn.pipeline import call_process
from smsn.fasta_tools import load_fasta_special, compute_chunk_infos, dict_to_fasta



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

    dict_to_fasta({ str(scaffold) : str(sequence) }, './chunked_ref.fasta', specify_HoleID = True)
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



def reshape_header(samstring):
    """from a .sam python string that contains the header (with the queries), will remove all the lines of the reference in the header that are not explicitely cited in the .sam.
    writes a warning in the log if more than one reference is contained in the .sam queries
    This step is mandatory when analyzing the hole alone because otherwise ipdSummary will shout something like 'I find more references in the header than in the real reference' ..."""

    logging.debug('[DEBUG] (reshape_header) Excision of the useless header references ')

    lines = samstring.split('\n')

    # 1 - List all the references in the queries and in the header

    list_ref_header = []
    list_ref_queries = []

    retour = ''

    first = True

    for line in lines:
        if line: # Won't consider empty lines of the last line
            if line[0] != '@': # If it is not a header line
                if first:
                    holenb = line.split()[0].split('/')[1]
                    first = False
                ref = str(line.split()[2]) # Look for the reference string
                if ref not in list_ref_queries: # If this reference hasn't been added yet (NB #FIXME --> Maybe using a set instead would be interesting, but I won't change the code because it works already well as it is)
                    list_ref_queries.append(ref)
                else:
                    pass
            else:
                if line[0:3] == '@SQ':
                    ref = str(line.split()[1]).split(':')[1] # It is situated after the 'SN:'
                    list_ref_header.append(ref)
    # logging.debug('[DEBUG] (reshape_header) The following scaffolds appear in the header : {}'.format(list_ref_header))
    logging.debug('[DEBUG] (reshape_header) The following scaffolds appear in the queries (hole n° {}) : {}'.format(holenb, list_ref_queries))

     # 2 - Make the list of the references that are in both

    list_to_keep = [x for x in list_ref_header if x in list_ref_queries] # Everything that is contained in both the header and the queries

    logging.debug('[DEBUG] (read_header) The following scaffolds will be kept in the header (hole n°) : {}'.format(holenb, list_to_keep))

    if len(list_to_keep) > 1: # In our case ipdSummary works with one reference only, which might cause a problem if it isn't the case
    # In case you are having trouble with this section and don't want to let pass sequences that map on more than 1 ref
    # then I think the better solution would be to raise an exception here (rather than printing a warning) and catch it upstream in YAD_perform_analysis_one_hole
        logging.warning('[WARNING] There is more than 1 reference in the queries and this would probably cause issues if you want to use ipdSummary for holenb {} : {}'.format(holenb, list_to_keep))

    # 3 - Generate a new sam file with the truncated header (references in the header that are not present in the queries are excised)

    newsam = '' # This is what we will return at the end of the function

    del (line) # I will re-use it just after and a generator cannot be properly reset

    first = True

    for line in lines: # Reminder: Lines are the splitted .sam string, splitted over the '\n'
        if line:
            if line[0] != '@': # If it is not a header
                newsam = newsam + '\n' + line # The '\n' has been deleted a few lines above, when I splitted the lines according to their \n so this is why i'm not doing newsam + line diretly
            else:
                if line[0:3] == '@SQ': # If it is a ref_seq in the header...
                    if line.split()[1].split(':')[1] in list_to_keep: # ... Then it is where the scaffold id is supposed to be found on @SQ lines (line.split [...]), and we look in the list_to_keep if we can find it
                        newsam = newsam + '\n' + line
                    else: # If its a @SQ header but not in the list_to_keep
                        pass
                else: # It is not a header, but it doesnt start with @SQ
                    if first: # To avoid having one '\n' at the beginning of the file that would cause samtools to deadlock
                        newsam = line
                        first = False
                    else: # If it's a header but it's not @SQ and it is not the first line
                        newsam = newsam + '\n' + line

    logging.debug('[DEBUG] (read_header) Excising the annoying headers is done, returning the newsam')

    return newsam