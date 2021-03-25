#!/usr/bin/env python

__author__ = "DELEVOYE Guillaume"
__credits__ = ["BAHIN Mathieu", "MEYER Eric"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "DELEVOYE Guillaume"
__email__ = "delevoye@ens.fr"
__status__ = "Developpment"


from smsn import bam_toolbox
from smsn.pipeline import call_process
import os
import logging
import gc


def launch_smsn_condor(args):
    """DEPRECATED: Launches the pipeline with HT_condor. See smsn --help for more info"""
    moviename = bam_toolbox.parse_header(args["bam"])["PU"]
    bam_header = bam_toolbox.read_header(args["bam"])
    args["moviename"] = moviename
    args["bam_header"] = bam_header

    assert args["CCS"]
    assert args["nb_proc"]
    assert args["sizechunks"]


    ####################################################################################################################
    logging.info("[INFO] Step 1 : Creating the consensus")
    if args["CCS"]:
        # This step is really expensive computationnally. We'll let the user give its own consensus if they're already
        # built
        args["CCS"] = os.path.realpath(args["CCS"])
        logging.info("[INFO] CCS {} will be used - CCS creation is skipped.".format(args["CCS"]))
    ####################################################################################################################
    logging.info("[INFO] Step 2 - 1: Aligning all the CCS against the reference using BLASR --- SKIPPED")

    args["aligned_CCS"] = args["CCS"]

    ####################################################################################################################
    logging.info('[INFO] Step 2 - 2: Loading the alignment produced')
    df_aligned_CCS = bam_toolbox.detailed_get_pd_dataframe_alignment(args["aligned_CCS"],include_ipds=False)

    try:
        assert len(df_aligned_CCS) == len(df_aligned_CCS['HoleID'].unique()) # Just checking that BLASR did everything right
    except AssertionError:
        logging.critical('[CRITICAL] There is not only one alignment per HoleID. The rest of the analysis might be '
                         'corrupted')

    output_alignmentcsv = os.path.join(args["tmpdir"],"aligned_CCS_"+args["moviename"]+".csv")
    logging.debug("[DEBUG] Saving alignment (.csv) at {}".format(output_alignmentcsv))
    df_aligned_CCS.to_csv(output_alignmentcsv,sep=";")

    ####################################################################################################################
    logging.info('[INFO] Step 2 - 3: Keeping only alignments with >= {} identity'.format(args["min_identity"]))
    df_aligned_CCS = df_aligned_CCS[df_aligned_CCS["identity_score"] >= args["min_identity"]]
    df_aligned_CCS.index = [int(x) for x in df_aligned_CCS["HoleID"]]

    set_holes_to_analyse = set(df_aligned_CCS["HoleID"].unique())

    logging.debug('[DEBUG] {} holes are kept after filtering upon identity (Before filtering we had {})'.format(len(set_holes_to_analyse),len(df_aligned_CCS)))


    ####################################################################################################################
    logging.debug("[DEBUG] Skipping pandarallel")
    ####################################################################################################################
    logging.debug("[DEBUG] Entering chunked data + Parallelism")
    logging.info("[INFO] Starting analysis")
    logging.debug("[DEBUG] Initializing reader")
    reader = bam_toolbox.yield_hole_by_hole(args["bam"],
                                                 header=None,
                                                 restricted_to=set_holes_to_analyse,
                                                 min_subreads=args["min_subreads"])
    logging.debug("[DEBUG] Reader initialized")

    logging.info('[INFO] Cutting the input .bam files into k chunks.')

    i = 0
    chunknb = 0

    bamfile_list = []

    tmp_sam = []
    for samtxt in reader:
        # logging.debug('[DEBUG] i = {}'.format(i))
        if i <= args["sizechunks"]:
            tmp_sam.append(samtxt.strip())
            # logging.debug("[DEBUG] Added samline : {}".format(samtxt.strip()))
            i+=1
        else:
            logging.debug('[DEBUG] Cutting .bam chunk n° {}'.format(chunknb))
            tmp_sam = args["bam_header"].strip() + '\n' + "\n".join(tmp_sam) + "\n"
            bamtxt = bam_toolbox.inmemory_asbam(tmp_sam)

            filename = os.path.join(os.path.realpath(args["tmpdir"]),moviename+"_chunknb"+str(chunknb)+"_HTcondor.bam")
            bamfile_list.append(filename)
            with open(filename,"wb") as bamoutfile:
                bamoutfile.write(bamtxt)

            i = 0
            chunknb +=1
            del tmp_sam
            gc.collect()
            tmp_sam = []

    if tmp_sam:
        logging.debug('[DEBUG] Cutting .bam chunk n° {}'.format(chunknb))
        tmp_sam = args["bam_header"].strip() + '\n' + "\n".join(tmp_sam) + "\n"
        bamtxt = bam_toolbox.inmemory_asbam(tmp_sam)

        filename = os.path.join(os.path.realpath(args["tmpdir"]),
                                moviename + "_chunknb" + str(chunknb) + "_HTcondor.bam")
        bamfile_list.append(filename)
        with open(filename, "wb") as bamoutfile:
            bamoutfile.write(bamtxt)

        i = 0
        chunknb += 1
        del tmp_sam
        gc.collect()
        # tmp_sam = []

    default_condorsubmit = 'concurrency_limits = limit_smsn:{} \nexecutable	=	bash_chunkid.sh \nnice_user	=	True ' \
                           '\noutput	=	bashchunkid.stdout \nerror	=	bashchunkid.stderr \nlog    =    '\
                           'bashchunkid.log\ninitialDir	=	HERE ' \
                           '\ngetenv    =    True\nrequest_cpus    =    1\nrequest_memory    =    3G\nqueue '.strip().format(10000//args["nb_proc"])

    HERE = os.getcwd()
    chunkid= 0
    for filename in bamfile_list:
        this_tmpdir = os.path.join(args["tmpdir"], 'Condor_TMPDIR'+str(chunkid))
        os.system("mkdir -p "+this_tmpdir)
        os.chdir(this_tmpdir)
        bashcmd = 'smsn --reference ' + args["reference"] + " --bamfile " + filename + " --output_csv " + os.path.join(args["tmpdir"], str(chunkid) + "_tmpchunk.csv") + " --min_identity 0.99 --nb_proc 1 --sizechunks " + str(args["sizechunks"]) + " --verbosity DEBUG --tmpdir " + this_tmpdir
        logging.debug('[DEBUG] writing bashcmd = {}'.format(bashcmd))
        logging.debug('[DEBUG] Destination = {}'.format(os.path.realpath("./bash_"+str(chunkid)+".sh")))
        bashfile = open("./bash_"+str(chunkid)+".sh","w")
        bashfile.write(bashcmd)
        bashfile.close()

        os.system('chmod +x '+"./bash_"+str(chunkid)+".sh")

        condor_submit_txt = default_condorsubmit.replace('chunkid',str(chunkid)).replace('HERE',os.path.realpath(os.getcwd()))
        with open('condorfile.submit','w') as filout:
            logging.debug('[DEBUG] Writing condor submit file at {}'.format(os.path.realpath('./condorfile.submit')))
            filout.write(condor_submit_txt)

        # This must be done otuside the "with" statement otherwise the file is still buffered and not really available
        logging.debug("[DEBUG] Executing condor_sbmit file : \n {}".format(condor_submit_txt))
        call_process('condor_submit {}'.format(os.path.realpath("./condorfile.submit")))


        os.chdir(HERE)
        chunkid +=1

    logging.info('[INFO] SMSN ended.')

