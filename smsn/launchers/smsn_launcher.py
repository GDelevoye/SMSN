#!/usr/bin/env python

import smsn
import argparse
import os
import sys
import logging

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
                            help='Choose the model for IPD prediction. See the package\'s README for more info. DEFAULT: SP2-C2 (PacBio Sequel I)',
                            required=False,
                            choices=["SP2-C2","C2","P4-C2","P5-C3","P6-C4","XL-C2","XL-XL"],
                            default="SP2-C2")


    # parser.add_argument("--alignment_CSV","-a",
    #                     help="smrt_alignment .csv \n"
    #                          "[IN THE FINAL PIPELINE, CCS will be recreated and aligned]",
    #                     required=True
    #                     )

    parser.add_argument("--frequency","-f",
                            help="Frequency of the sequencer (Hz). Default [AUTO = Will parse the .bam to find it]. Suggested DEFAULT value : 80 Hz.",
                            required=True,
                            default = 80,
                            type = check_positive
                        )

    parser.add_argument('--output_csv',"-o",
                            help='Ouput file (csv) of the methylation analysis. See the README for further details on the output\'s format.',
                            required=True)

    parser.add_argument('--tmpdir','-t',
                        help="Tmp directory [DEFAULT = Current]",
                        required=False)

    parser.add_argument('--verbosity',"-v",
                            help='Choose your verbosity. Default: INFO',
                            required=False,
                            default="INFO",
                            choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])

    parser.add_argument('--progress_bar',"-p",
                            help='Displays a progress bar',
                            required=False,
                            action='store_true')

    parser.add_argument('--nproc','-n',
                        help="Max number of processors for parallelism.",
                        required=False,
                        default=1,
                        type=check_positive)

    args = parser.parse_args()

    fastafile = os.path.realpath(args.fastafile)
    bamfile = os.path.realpath(args.bam)
    output_csv = os.path.realpath(args.output_csv)
    # alignment_csv = os.path.realpath(args.alignment_csv)

    if args.tmpdir:
        tmpdir = os.path.realpath(tmpdir)
    else:
        tmpdir = os.getcwd()


    verboselevel = "logging."+str(args.verbosity)
    logging.basicConfig(level=eval(verboselevel),
                        format='%(asctime)s %(message)s',
                        stream=sys.stdout)


    show_progress_bar = False
    if args.progress_bar:
        show_progress_bar = True

    logging.critical("[CRITICAL] Not implemented")
    # ipdtools.ipdModel.compute_fasta_to_csv(modelname=args.model, fastafile=fastafile, csv_out=output_csv,
    #                                        show_progress_bar=show_progress_bar,nproc=args.nproc,indexing=args.indexing)

if __name__ == "__main__":
    main()
