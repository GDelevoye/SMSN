#!/usr/bin/env python

__author__ = "DELEVOYE Guillaume"
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

def main():
    """Parses the user's arguments, makes all the required checkings and warnings (if needed), and then launches the
    pipeline """

    parser = argparse.ArgumentParser()

    def check_positive(value):
        ivalue = int(value)
        if ivalue <= 0:
            raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
        return ivalue

    def check_positive_or_null(value):
        ivalue = int(value)
        if ivalue < 0:
            raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
        return ivalue


    parser.add_argument("--bam","-b",
                            help="Path to a .bam file with all the subreads (adapters sequences must be already "
                                 "removed)",
                            required=True
                        )

    parser.add_argument("--reference","-r",
                            help='Path to a genome reference (fasta file).',
                            required=True,
                            default=None)

    parser.add_argument("--model","-m",
                            help='Choose the model for IPD prediction.',
                            required=False,
                            choices=["SP2-C2","SP3-C3","P6-C4"],
                            default="SP2-C2")
    #
    # parser.add_argument("--strategies","-s",
    #                         help="Strategies to investivate DNA methylation.'
    #                              'PacBio: Default PacBio scores'
    #                              'Beulaurier: Beaulaurier 2015 SMsn scores'
    #                              'MannWhitney_sequential: Non-parametric and Sequential analysis based on median IPD',
    #                         required=False,
    #                         choices=["PacBio","Beaulaurier","SequentialMannWhitney"],
    #                         nargs="+",
    #                         default="auto")
    #
    # parser.add_argument('--MannWhitneySequential',
    #                         help='Minimum threshold on the sequential Mann-Whitney test',
    #                         required=False,
    #                         default=0.01,
    #                         type=float)
    #
    # parser.add_argument("--frequency","-f",
    #                         help="[NOT MANDATORY] Frequency of the sequencer (Hz). "
    #                              "This will only be used if strategy 'beaulaurier' or 'SequentialMannWhitney' are used"
    #                              "Default = AUTO. In the tested datasets, Sequel I had 80Hz and RSII had 75Hz.",
    #                         required=False,
    #                         default = "PacBio")

    parser.add_argument('--output_csv',"-o",
                            help='Ouput file (csv) of the methylation analysis. See the README for further details on '
                                 'the output\'s format.',
                            required=True)

    parser.add_argument("--CCS","-c",
                        help="[FACULTATIVE] Path to the circular consensus corresponding to the .bam subreads. "
                             "Default = CCS will be recreated from scratch.",
                        required=False,
                        default=None)

    parser.add_argument('--min_identity',"-i",
                            help='minimum identity (percentage) of the CCS required to launch analysis on a hole.',
                            required=False,
                            default=0.99,
                            type=float)


    parser.add_argument('--min_subreads',
                            help='Minimum number of subreads required to launch analysis on the hole. DEFAULT = 50 ('
                                 'so that its possible to have >=25X per strand on at least one position) ',
                            required=False,
                            default=50,
                            type=check_positive_or_null)


    parser.add_argument('--tmpdir','-t',
                        help="Tmp directory",
                        default="smsn_tmpdir_[DATE_HOUR]",
                        required=False)


    parser.add_argument('--verbosity',"-v",
                            help='Choose your verbosity. Default: DEBUG',
                            required=False,
                            default="INFO",
                            choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])

    parser.add_argument('--progress_bar',"-p",
                            help='Displays a progress bar. Disabled automatically if verbosity is set to debug'\
                                 '[DEFAULT]: Set to True if DEBUG as verbosity, otherwise FALSE',
                            default=False,
                            required=False,
                            action='store_true')

    parser.add_argument('--nb_proc',"-n",
                            help="""Multiprocessing on n CPU. Default: 1.""",
                            required=False,
                            default=1,
                            type=check_positive)

    parser.add_argument('--sizechunks',"-k",
                            help="""Because the subreads .bam will often not fit entirely in RAM and because the 
                            methylation analysis itself generates lots of data, SMSN will pause and perform intensive 
                            I/O operation every S holes it has analyzed. Lower values are better for machines that 
                            are limited in RAM. The optimal nb_proc/sizehunks ratio will vary from one computer to 
                            another. In case sizechunks < nb_proc, SMSN will use sizechunks = 20x nb_proc instead.""",
                            required=False,
                            default=5000,
                            type=check_positive)

    parser.add_argument('--add_context',
                            help='In the output .csv file, displays the +12/-12 context around the nucleotide.'
                                 '(Files generated can be sensitively heavier) [DEFAUTL: TRUE]',
                            default=True,
                            required=False,
                            type=bool)

    parser.add_argument('--preserve_tmpdir',
                            help="""Forbids deletion of tmp dir (experimental / deprecated / debug only)""",
                            required=False,
                            default=False,
                            action="store_true")

    parser.add_argument('--idQvs',
                            help="""Outputs PacBio's identificationQV [DEFAULT: TRUE]
                            A value of -1 means that the identificationQv was not computed by the PacBio softwares""",
                            required=False,
                            default=True,
                            type=bool)


    args = parser.parse_args()

    assert args.min_identity > 0 and args.min_identity <= 1

    ####### Handling verbosity


    if args.verbosity != "NONE": # Configure the loglevel according to user's preference
        verboselevel = "logging."+str(args.verbosity)
        logging.basicConfig(stream=sys.stdout, level=eval(verboselevel),
                            format='%(asctime)s %(message)s')
    if args.verbosity == "DEBUG":
        args.progress_bar = True
    elif args.progress_bar:
        show_progress_bar = True

    ######## Passing all the abspath arguments into real paths

    args.reference = os.path.realpath(args.reference)
    args.bamfile = os.path.realpath(args.bam)
    args.output_csv = os.path.realpath(args.output_csv)

    ######## Handling framerate of the sequencer
    # This is basically just warning the user if it's frameratehz are a bit weird.

    # IF "Auto" is not selected

    # if args.frequency != "auto":
    #     try: # Even if not in auto mode, we'll still look into the .bam just to be sure
    #         bam_framerate = smsn.bam_toolbox.parse_header(bamfile)[
    #             "FRAMERATEHZ"]  # This is the framerate given in the PacBio input file
    #
    #         if int(args.frequency) != int(bam_framerate): # Make sure that the user really wants to use another framerate than what's indicated because it's weird
    #             logging.warning('[WARNING] User framerate(HZ) and bam framerate(HZ) are not identical ({} VS {})'.format(bam_framerate,int(args.frequency)))
    #             logging.warning("[WARNING] Framerate {} is used instead of {}".format(int(bam_framerate),int(args.frequency)))
    #             args.frequency = float(args.frequency) # Respect the user's will anyway
    #         else: # If the user gave the same framerate as what's in the .bam then just stay quiet and go on it's ok
    #             args.frequency = float(args.frequency)
    #
    #     except Exception as e:
    #         # Since the user specified the framerate explicitely, there is probably a reason for that
    #         # The program should keep going on, no matter what exception has happened before
    #         # In case it fails reading the framerate of the .bam file:
    #         logging.error("[ERROR] Exception : \n{}".format(e)) # Users should still be warned that something went wrong because that's weird
    #         logging.warning('[WARNING] Failed to inder framerate from input bamfile')
    #         logging.warning('[WARNING] Custom framerate will be used: {}.'.format(args.frequency)) # Indicate that we respect their will after all, even if that's a strange request
    #         args.frequency = float(args.frequency)
    #
    # # When in auto mode, things are easier
    # else:
    #     try: # Just try to read the Framerate from the file
    #         bam_framerate = smsn.bam_toolbox.parse_header(bamfile)[
    #             "FRAMERATEHZ"]  # This is the framerate given in the PacBio input file
    #         args.frequency = float(bam_framerate)
    #
    #     except Exception as e: # If it doesn't work, indicate the problem clearly to the user
    #         logging.CRITICAL("[CRITICAL] Exception : \n{}".format(e))
    #         # Propose a solution
    #         logging.critical("[CRITICAL] Can't read the framerate from .bamfile. The .bam file might be corrupted or "
    #                          "ill-formated. You can try to specify your framerate manually using --frequency. Usual "
    #                          "framerate are around 75-80HZ for respectively RSII and Sequel I sequencers")
    #         # Then quit, it's not worth going any further if methylation analysis is impossible anyway
    #         exit(-1)
    #
    # assert args.frequency > 0 # Seems legitimate

    ########## HANDLING CHEMISTRY

    # args.model = handle_chemistry_and_compatibility(args.model, bamfile)

    ########## HANDLING tmp_dir

    logging.debug('[DEBUG] tmpdir argument recieved = {}'.format(args.tmpdir))

    if args.tmpdir == "smsn_tmpdir_[DATE_HOUR]":
        now = datetime.now()
        # dd/mm/YY H:M:S
        # time.sleep(1)
        dt_string = now.strftime("%d_%m_%Y_%H_%M_%S")
        tmpdirname = os.path.join(os.path.realpath(os.getcwd()),"smsn_tmpdir_"+dt_string)
        args.tmpdir = os.path.realpath(tmpdirname)

    else:
        tmpdirname = os.path.realpath(args.tmpdir)
        args.tmpdir = os.path.realpath(tmpdirname)

    logging.debug("[DEBUG] Creating tmpdir = {}".format(args.tmpdir))

    try:
        logging.debug("[DEBUG] Creating {}".format(args.tmpdir))
        os.mkdir(args.tmpdir)
    except OSError as error:
        logging.error("[ERROR] {}".format(error)) # This can occur, especially if many several of smsn with default tmpdir are launched simultaneously
        # (Because they would have the same outputdir_name)
        # When this happen we can try another time later:
        logging.warning(
            "[WARNING] It seems that another tmpdir has already been produced at the same time. Waiting 2s to try again")

        try:
            time.sleep(2)
            os.mkdir(args.tmpdir)
            logging.debug("[DEBUG] Seems alright finally. Tmpdir created: {}".format(args.tmpdir))
        except OSError as error:
            logging.critical("[CRITICAL] {}".format(error)) # Then it's bad
            raise OSError

    ###### Handling processor informations

    try:
        nb_available = psutil.cpu_count(logical=True)
        logging.info("[INFO] Available cores: {}".format(nb_available))
    except:
        nb_available = 1
        logging.error("[ERROR] Failed to seek the available number of CPU cores with psutils. Max is set to 1")

    if int(args.nb_proc) == -1:
        args.nb_proc = nb_available
        logging.info("[INFO] Using all known cores: {} CPU cores.".format(args.nb_proc))

    elif int(args.nb_proc) < 1:
        logging.warning("[WARNING] Number of CPU cores cannot be set to null or negative values. Will be set to 1")
        args.nb_proc = 1

    elif args.nb_proc > nb_available:
        logging.warning("[WARNING] You asked for {} CPU cores but only {} are available. n_proc is set to {}".format(args.nb_proc,nb_available,nb_available))
        args.nb_proc = nb_available

    else: # If there is no problem
        logging.info("[INFO] Use CPU cores: {}".format(args.nb_proc))

    logging.debug("[DEBUG] Launching the analysis with the following arguments: {}".format(args))

    if args.sizechunks < args.nb_proc:
        logging.warning("[WARNING] sizechunk being < to nb_proc, we'll switch to sizechunk = np_proc x 20 by default")
        args.sizechunks = (args.nb_proc) * 20

    # And finally when everything is set up correctly, we just have to do, very simply:
    smsn.pipeline.launch_smsn(vars(args))

if __name__ == "__main__":
    main()
