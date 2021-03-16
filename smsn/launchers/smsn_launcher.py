#!/usr/bin/env python

__author__ = "DELEVOYE Guillaume"
__credits__ = ["BAHIN Mathieu", "MEYER Eric"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "DELEVOYE Guillaume"
__email__ = "delevoye@ens.fr"
__status__ = "Developpment"

import smsn
from datetime import datetime
import time
import argparse
import os
import sys
import logging
import psutil

def handle_chemistry_and_compatibility(modeluser,bamfile):
    """Returns the chemistry to use. Logs errors and warnings if potential compatibility problems are detected, or if we
    can give some advices."""
    logging.debug("[DEBUG] Handling the ipdmodel of user compared to the type of its data")

    picked = ""

    logging.warning("[WARNING] SMSN is not available to parse all the subtle differences between the Sequel Versions. Please ensure that you've read https://github.com/GDelevoye/ipdtools#-which-model-should-i-use-")

    parsed_header = smsn.bam_toolbox.parse_header(bamfile)
    if modeluser != "auto":
        logging.warning("[WARNING] The model is not selected automatically, which is deprecated outside of very precise studies.")
        picked = modeluser
    elif modeluser == "auto":
        if parsed_header["PM"] == "SEQUEL":
            logging.info("[INFO] Sequel data detected. Model SP2-C2 will be used.")
            picked = "SP2-C2"
        elif parsed_header["PM"] == "RS":
            logging.info("[INFO] RS data detected. Model P5-C3 will be used.")
            picked = "P5-C3"
        else:
            logging.error("[ERROR] Platform unknown : {} (never tested). SP2-C2 (Sequel I/II.2) will be used as default. Results could be totally wrong.".format(parsed_header["PM"]))
            picked = "SP2-C2"
            logging.info("[INFO] Trying to use SP2-C2")

    return picked

def main():

    parser = argparse.ArgumentParser()

    def check_positive(value):
        ivalue = int(value)
        if ivalue <= 0:
            raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
        return ivalue

    parser.add_argument("--bam","-b",
                            help="Path to a .bam file with all the subreads (adapters sequences must be already removed)",
                            required=True
                        )

    parser.add_argument("--reference","-r",
                            help='Path to a genome reference (fasta file).',
                            required=True,
                            default=None)

    parser.add_argument("--model","-m",
                            help='Choose the model for IPD prediction. See https://github.com/GDelevoye/ipdtools#-which-model-should-i-use- for more info. DEFAULT: Automatic guess according to the input .bam',
                            required=False,
                            choices=["SP2-C2","C2","P4-C2","P5-C3","P6-C4","XL-C2","XL-XL","auto"],
                            default="auto")

    parser.add_argument("--frequency","-f",
                            help="[NOT MANDATORY] Frequency of the sequencer (Hz). Default = AUTO. In the tested datasets, Sequel I had 80Hz and RSII had 75Hz.",
                            required=False,
                            default = "auto")

    parser.add_argument('--output_csv',"-o",
                            help='Ouput file (csv) of the methylation analysis. See the README for further details on the output\'s format.',
                            required=True)

    parser.add_argument('--tmpdir','-t',
                        help="Tmp directory",
                        default="smsn_tmpdir_[DATE_HOUR]",
                        required=False)


    parser.add_argument('--verbosity',"-v",
                            help='Choose your verbosity. Default: INFO',
                            required=False,
                            default="INFO",
                            choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])

    parser.add_argument('--progress_bar',"-p",
                            help='Displays a progress bar',
                            default=False,
                            required=False,
                            action='store_true')

    parser.add_argument('--n_proc',"-n",
                            help="""Multiprocessing on n CPU. Default: 1.
                            If set to -1, will use all available cores.""",
                            required=False,
                            default=1,
                            type=check_positive)

    args = parser.parse_args()

    fastafile = os.path.realpath(args.fastafile)
    bamfile = os.path.realpath(args.bam)
    output_csv = os.path.realpath(args.output_csv)

    ######## Handling framerate

    bam_framerate = smsn.bam_toolbox.parse_header(bamfile)["FRAMERATEHZ"]
    if args.frequency != "auto":
        if int(args.frequency) != int(bam_framerate):
            logging.error('[ERROR] User framerate(HZ) and bam framerate(HZ) are not identical ({} VS {})'.format(bam_framerate,int(args.frequency)))
            logging.warning("[WARNING] Framerate {} is used instead of {}".format(int(bam_framerate),int(args.frequency)))
            args.frequency = float(args.frequency)
        else:
            args.frequency = float(args.frequency)

    else:
        args.frequency = float(bam_framerate)

    assert args.frequency > 0

    ########## HANDLING CHEMISTRY

    handle_chemistry_and_compatibility(args.model, bamfile)

    ########## HANDLING tmp_dir



    ##### Handling outputdir name

    if args.tmpdir == "smsn_tmpdir_[DATE_HOUR]":
        now = datetime.now()
        # dd/mm/YY H:M:S
        # time.sleep(1)
        dt_string = now.strftime("%d_%m_%Y_%H_%M_%S")
        outdir = "smsn_tmpdir_"+dt_string
        tmpdirname = os.path.join(os.path.realpath(os.getcwd()),args.tmpdir)
        args.tmpdir = tmpdirname

    else:
        tmpdirname = os.path.realpath(args.tmpdir)
        args.tmpdir = tmpdirname

    logging.debug("[DEBUG] Creating tmpdir = {}".format(args.tmpdir))

    try:
        os.mkdir(args.tmpdir)
    except OSError as error:
        logging.error("[ERROR] {}".format(error)) # This can occur, especially if many instances of pt_sort are launched simultaneously
        # (Because they would have the same outputdir_name)
        # When this happen we can try another time later:
        try:
            time.sleep(2)
            os.mkdir(args.tmpdir)
        except OSError as error:
            logging.critical("[CRITICAL] {}".format(error)) # Then it's bad
            raise OSError
            # exit(-1) # And the user should be alerted


    if args.tmpdir:
        tmpdir = os.path.realpath(args.tmpdir)
    else:
        tmpdir = os.getcwd()

    ####### Handling verbosity


    if args.verbosity != "NONE": # Configure the loglevel according to user's preference
        verboselevel = "logging."+str(args.verbosity)
        logging.basicConfig(stream=sys.stdout, level=eval(verboselevel),
                            format='%(asctime)s %(message)s')

    if args.progress_bar:
        show_progress_bar = True


    ###### Handling processor informations

    try:
        nb_available = psutil.cpu_count(logical=True)
        logging.info("[INFO] Available cores: {}".format(nb_available))
    except:
        nb_available = 1
        logging.error("[ERROR] Failed to seek the available number of CPU cores with psutils. Max is set to 1")

    if int(args.n_proc) == -1:
        args.n_proc = nb_available
        logging.info("[INFO] Using all known cores: {} CPU cores.".format(args.n_proc))

    elif int(args.n_proc) < 1:
        logging.warning("[WARNING] Number of CPU cores cannot be set to null or negative values. Will be set to 1")
        args.n_proc = 1

    elif args.n_proc > nb_available:
        logging.warning("[WARNING] You asked for {} CPU cores but only {} are available. n_proc is set to {}".format(args.n_proc,nb_available,nb_available))
        args.n_proc = nb_available

    else: # If there is no problem
        logging.info("[INFO] Use CPU cores: {}".format(args.n_proc))

    logging.critical("[CRITICAL] Not implemented")

if __name__ == "__main__":
    main()
