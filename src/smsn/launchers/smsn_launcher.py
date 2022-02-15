#!/usr/bin/env python

__author__ = "DELEVOYE Guillaume"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "DELEVOYE Guillaume"
__email__ = "delevoye@ens.fr"
__status__ = "Developpment"

import logging
import smsn
from datetime import datetime
import argparse
from argparse import HelpFormatter
import os
import sys
import textwrap
import psutil
import copy

# https://stackoverflow.com/questions/3853722/how-to-insert-newlines-on-argparse-help-text
class RawFormatter(HelpFormatter):
    """Ensures that the description text of the parser fits to standards in terms of width and indents"""

    def _fill_text(self, text, width, indent):
        # Strip the indent from the original python definition that plagues
        # most of us.
        text = textwrap.dedent(text)
        text = textwrap.indent(text, indent)  # Apply any requested indent.
        text = text.splitlines()  # Make a list of lines
        text = [textwrap.fill(line, width) for line in text]  # Wrap each line
        text = "\n".join(text)  # Join the lines again
        return text


def handle_cpu(args):
    """Handles the CPU number and the chunksize to ensure that paralellism works correctly

    >>> from dataclasses import dataclass
    >>> import logging
    >>> @dataclass
    ... class MockerArgparse:
    ...     nb_proc: int
    ...     sizechunks: int
    ...
    >>> mocked_args = MockerArgparse(12,5000)
    >>> logging.basicConfig(stream=sys.stdout, level="ERROR",format='%(asctime)s %(message)s')
    >>> # To avoid logging messages to pollute the docstring test
    >>> import psutil
    >>> psutil.cpu_count = lambda logical: 12 # Simulate platform with 12 CPU
    >>> psutil.cpu_count(logical=True)
    12
    >>> handle_cpu(mocked_args).nb_proc
    12
    >>> psutil.cpu_count = lambda logical: 6 # Simulate platform with 6 CPU
    >>> psutil.cpu_count(logical=True)
    6
    >>> mocked_args = MockerArgparse(12,5000) # When asking too much, a capping should occur
    >>> handle_cpu(mocked_args).nb_proc
    6
    >>> mocked_args = MockerArgparse(5,5000) # Otherwise, same value
    >>> handle_cpu(mocked_args).nb_proc
    5
    >>> mocked_args = MockerArgparse(-1,5000) # Otherwise, same value
    >>> handle_cpu(mocked_args).nb_proc == psutil.cpu_count(logical=True)
    True
    >>> mocked_args = MockerArgparse(-2,5000)
    >>> handle_cpu(mocked_args).nb_proc
    1
    >>> logging.basicConfig(stream=sys.stdout, level="WARNING",format='%(message)s')
    >>> mocked_args = MockerArgparse(6,5)
    >>> results = handle_cpu(mocked_args)
    >>> psutil.cpu_count = lambda logical: 6 # Simulate platform with 6 CPU
    >>> results.sizechunks == 20 * results.nb_proc
    True
    """
    # Handling processor informations

    try:
        nb_available = psutil.cpu_count(logical=True)
        logging.info("[INFO] Available cores: {}".format(nb_available))
    except BaseException:
        nb_available = 1
        logging.error(
            "[ERROR] Failed to seek the available number of CPU cores with psutils. Max is set to 1")

    if int(args.nb_proc) == -1:
        args.nb_proc = nb_available
        logging.info(
            "[INFO] Using all known cores: {} CPU cores.".format(
                args.nb_proc))

    elif int(args.nb_proc) < 1:
        logging.warning(
            "[WARNING] Number of CPU cores cannot be set to null or negative values. Will be set to 1")
        args.nb_proc = 1

    elif args.nb_proc > nb_available:
        logging.warning(
            "[WARNING] You asked for {} CPU cores but only {} are available. n_proc is set to {}".format(
                args.nb_proc, nb_available, nb_available))
        args.nb_proc = nb_available

    else:  # If there is no problem
        logging.info("[INFO] Use CPU cores: {}".format(args.nb_proc))

    logging.debug(
        "[DEBUG] Launching the analysis with the following arguments: {}".format(args))

    if args.sizechunks < args.nb_proc:
        logging.warning(
            "[WARNING] sizechunk being < to nb_proc, we'll switch to sizechunk = np_proc x 20 by default")
        args.sizechunks = args.nb_proc * 20

    args_copy = copy.deepcopy(args)
    return args_copy


def check_positive(value):
    """Checks that a value is strictly positive when converted to int

    >>> try:
    ...     check_positive(float(0.1))
    ... except Exception as e:
    ...     print(True)
    ...
    True
    >>> try:
    ...     check_positive(int(0))
    ... except Exception as e:
    ...     print(True)
    ...
    True
    >>> check_positive(3)
    3"""

    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(
            "%s is an invalid (strictly) positive int value" % value)
    return ivalue

def check_pos_or_zero(value):
    """ Checks that an integer value is positive or null

    >>> try:
    ...     check_pos_or_zero(int(-1))
    ... except Exception as e:
    ...     print(True)
    ...
    True
    >>> try:
    ...     check_pos_or_zero(int(0))
    ... except Exception as e:
    ...     print(True)
    ...
    0
    >>> check_pos_or_zero(int(10))
    10
    >>> try:
    ...     check_pos_or_zero(int(0.1))
    ... except Exception as e:
    ...     print(True)
    ...
    0
    """
    ivalue = int(value)
    if ivalue < 0:
        raise argparse.ArgumentTypeError(
            "%s is an invalid positive (or null) int value" % value)
    return ivalue

def between_0_and_1(value) -> float:
    """Checks that a value is in [0,1]

    >>> try:
    ...     foo = between_0_and_1(float(0.1))
    ... except Exception as e:
    ...     print("error")
    ...
    >>> try:
    ...     foo = between_0_and_1(float(1))
    ... except Exception as e:
    ...     print("error")
    ...
    >>> try:
    ...     foo = between_0_and_1(float(0))
    ... except Exception as e:
    ...     print("error")
    ...
    >>> try:
    ...     foo = between_0_and_1(float(10))
    ... except Exception as e:
    ...     print("error")
    ...
    error
    >>> try:
    ...     foo = between_0_and_1(float(-4))
    ... except Exception as e:
    ...     print("error")
    ...
    error
    >>> try:
    ...     foo = between_0_and_1(float(0.5))
    ... except Exception as e:
    ...     print(error)
    ...
    """
    ivalue = float(value)
    if ivalue < 0 or ivalue > 1:
        raise argparse.ArgumentTypeError(
            "%s is an invalid value between ]0:1]" % value)

    return ivalue

def str2bool(v):
    """
    >>> str2bool("YES")
    True
    >>> str2bool("YeS")
    True
    >>> str2bool("1")
    True
    >>> str2bool("y")
    True
    >>> str2bool("No")
    False
    >>> str2bool("None")
    False
    >>> str2bool("nope")
    Traceback (most recent call last):
     ...
    argparse.ArgumentTypeError: Boolean value expected.
    >>> str2bool("N")
    False
    >>> str2bool("nO")
    False
    >>> str2bool(True)
    True
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0', "none"):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def parseargs(args):

    description = """This software implements the SMSN approach of PacBio sequencing in the same way as Beaulaurier
    et al 2015, but contrary to Beaulaurier's software, this one can be used even in the absence of a PCR-amplified
    control, and on more recent data (tested on Sequel I).

    DELEVOYE Guillaume 2020.
    thttps://github.com/EMeyerLab/SMSN
    """

    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawFormatter)

    # Otherwise optional arguments can preceed required ones, which is absurd
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')


    # Required arguments >

    required.add_argument("--bam", "-b",
                          help="Path to a .bam file with all the subreads (adapters sequences must be already "
                               "removed)",
                          required=True
                          )

    required.add_argument("--reference", "-r",
                          help='Path to a genome reference (fasta file).',
                          required=True,
                          default=None)

    required.add_argument("--model", "-m",
                          help="Choose the model for IPD prediction. [DEFAULT: auto (PacBio's kineticsTools software"
                               "choses after it has parsed the input file)]",
                          required=True,
                          choices=["SP2-C2", "SP3-C3", "P6-C4"])

    required.add_argument('--output_csv', "-o",
                          help='Ouput file (csv) of the methylation analysis. See the README.md for further details on '
                               'the output\'s format.',
                          required=True)

    # Optional arguments >

    optional.add_argument("--CCS", "-c",
                          help="[FACULTATIVE] Path to the circular consensus corresponding to the .bam subreads. "
                               "Default = CCS will be recreated from the subreads provided.",
                          required=False,
                          default=None)

    optional.add_argument('--min_identity', "-i",
                          help='minimum identity (percentage) of the CCS required to launch analysis on a hole.'
                               '[DEFAULT : 0.99]. Must be in [0;1]',
                          required=False,
                          default=0.99,
                          type=between_0_and_1)

    optional.add_argument('--min_subreads',
                          help='Minimum number of subreads required to launch analysis on the hole. DEFAULT = 50 ('
                               'so that its possible to have >=25X per strand on at least one position).',
                          required=False,
                          default=50,
                          type=check_pos_or_zero)

    optional.add_argument('--tmpdir', '-t',
                          help="Tmp directory (DEFAULT : smsn_tmpdir_[DATE_HOUR])",
                          default="smsn_tmpdir_[DATE_HOUR]",
                          required=False)

    optional.add_argument('--verbosity', "-v",
                          help='Choose your verbosity. Default: INFO',
                          required=False,
                          default="INFO",
                          choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])

    optional.add_argument('--progress_bar', "-p",
                          help='Displays a progress bar. Disabled automatically if verbosity is set to debug'
                               '[DEFAULT]: False',
                          default=False,
                          type=str2bool)

    optional.add_argument('--nb_proc', "-n",
                          help="""Multiprocessing on n CPU. Default: 1.""",
                          required=False,
                          default=1,
                          type=check_positive)

    optional.add_argument('--sizechunks', "-k",
                          help="""The subreads will often not fit entirely in RAM and the methylation analysis itself
                            generates lots of I/O usage. Because of it, smsn will pause and perform intensive
                            I/O operation every S holes it has analyzed. Lower values are better for machines that
                            are limited in RAM. The optimal nb_proc/sizehunks ratio will vary from one computer to
                            another. In case sizechunks < nb_proc, SMSN will use sizechunks = 20x nb_proc instead.

                            DEFAULT : 5000""",
                          required=False,
                          default=5000,
                          type=check_positive)

    optional.add_argument('--add_context',
                          help='In the output .csv file, displays the +12/-12 context around the nucleotide.'
                               '(Files generated can be sensitively heavier) [DEFAULT: True, choices = True or False]',
                          default=True,
                          required=False,
                          type=str2bool)

    optional.add_argument('--idQvs',
                          help="""Outputs PacBio's identificationQV [DEFAULT: True, choices = True or False]""",
                          required=False,
                          default=True,
                          type=str2bool)

    return parser.parse_args(args)

def main(testargs=None):
    """Launches the parser, Handles verbosity, paths, tmpdir, and finally launches the pipeline
    testargs : Serves to simulate fake sys.argv[1:] to test the CLI

    Optionnaly, this can also serve to call the software directly from python if the module is imported"""


    # Parsing the arguments

    if not testargs:
        args = parseargs(sys.argv[1:])
    else: # To test the CLI programatically
        args = parseargs(testargs)

    # Handling verbosity

    if args.verbosity != "NONE":  # Configure the loglevel according to user's preference
        verboselevel = "logging." + str(args.verbosity)
        logging.basicConfig(stream=sys.stdout, level=eval(verboselevel),
                            format='%(asctime)s %(message)s')


    logging.debug('[DEBUG] Verbosity set to {}'.format(args.verbosity))

    # Handling CPU

    args = handle_cpu(args)

    # All the abspath arguments must be transformed into real paths

    args.reference = os.path.realpath(args.reference)
    args.bam = os.path.realpath(args.bam)
    args.output_csv = os.path.realpath(args.output_csv)

    # HANDLING tmp_dir

    logging.debug('[DEBUG] tmpdir argument recieved = {}'.format(args.tmpdir))

    # Create an automated name if the user did not specify any
    if args.tmpdir == "smsn_tmpdir_[DATE_HOUR]":
        now = datetime.now()
        dt_string = now.strftime("%d_%m_%Y_%H_%M_%S")
        tmpdirname = os.path.join(
            os.path.realpath(
                os.getcwd()),
            "smsn_tmpdir_" +
            dt_string)
        args.tmpdir = os.path.realpath(tmpdirname)

    # Otherwise, create the tmp_dir as the user asks
    else:
        tmpdirname = os.path.realpath(args.tmpdir)
        args.tmpdir = os.path.realpath(tmpdirname)

    logging.debug("[DEBUG] Creating tmpdir = {}".format(args.tmpdir))

    try:
        logging.debug("[DEBUG] Creating {}".format(args.tmpdir))
        os.mkdir(args.tmpdir)
    except OSError as error:
        # This can occur, especially if many several of smsn with default
        # tmpdir are launched simultaneously
        logging.error("[ERROR] {}".format(error))

    # And finally...
    smsn.pipeline.launch_smsn(vars(args))


if __name__ == "__main__":
    main()
