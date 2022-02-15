"""Implementation of the pipeline itself (preprocessing,
create CCS, align them, launch the methylation analysis)"""

__author__ = "DELEVOYE Guillaume"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "DELEVOYE Guillaume"
__email__ = "delevoye@ens.fr"
__status__ = "Production"

import subprocess
import os
import smsn
import logging
import pandas as pd
from pandarallel import pandarallel
from pkg_resources import Requirement, resource_filename
import gc
import typing


def _getAbsPath(fname: str) -> os.path.realpath:
    """Returns the path of the package

    >>> import smsn
    >>> smsn.__file__.split('/')[-2] == "smsn"
    True
    >>> smsn.pipeline._getAbsPath("resources").split('/')[-2::]
    ['smsn', 'resources']
    >>> import os
    >>> os.path.isdir(_getAbsPath("resources"))
    True
    """
    return resource_filename(Requirement.parse('smsn'), 'smsn/%s' % fname)


def transform_model_name(modelname: str) -> os.path.realpath:
    """Returns the path to a PacBio model (among SP2-C2, P6-C4, SP3-C3)

    >>> import smsn
    >>> import os
    >>> listmodels = sorted(os.listdir(smsn.pipeline._getAbsPath("resources/models")))
    >>> listmodels
    ['P6-C4.npz.gz', 'SP2-C2.npz.gz', 'SP3-C3.npz.gz']
    >>> listmodels = [os.path.join(smsn.pipeline._getAbsPath("resources/models"),x) for x in listmodels]
    >>> listmodels[0] == smsn.pipeline.transform_model_name("P6-C4")
    True
    >>> listmodels[1] == smsn.pipeline.transform_model_name("SP2-C2")
    True
    >>> listmodels[2] == smsn.pipeline.transform_model_name("SP3-C3")
    True
    """
    resources_dir = _getAbsPath("/resources/models/")
    modifiedmodelname = modelname + ".npz.gz"
    modelname = os.path.join(resources_dir, modifiedmodelname)
    return modelname


def call_process(cmd: str) -> typing.Tuple[str, str]:
    """Simple wrapper to call the PacBio tools. The tests do not only check the function itself,
    but also the presence of these tools.

    Returns a tuple (stdout, stderr)

    Since PacBio tools abuse of stderr, a re-attribution is made and some errors are redirected
    to stdout

    >>> import smsn
    >>> import sys
    >>> liste_cmd = ["ccs --version","pbindex --version","ipdSummary --version","blasr --version"]
    >>> for mycmd in liste_cmd:
    ...     print(smsn.pipeline.call_process(mycmd)[0].__repr__())
    ...
    'ccs 3.1.0 (commit v3.1.0)'
    'pbindex 1.3.0'
    '3.0'
    'blasr\\t5.3.5-5.3.5'
    """

    return_stdout = []
    return_stderr = []

    processname = cmd.split()[0]
    logging.debug(
        '[DEBUG] (call_process) Processname = {}'.format(processname))
    logging.debug('[DEBUG] (call_process) cmd = {}'.format(cmd))

    process = subprocess.run(cmd, shell=True, capture_output=True)
    std_output = process.stdout.decode('utf-8')
    error_output = process.stderr.decode('utf-8')

    if std_output:
        if "ERROR" in std_output.upper():
            logging.error(
                '[ERROR] (stdout of process {} yielded an error: {}'.format(
                    processname, std_output.strip()))
            return_stderr.append(std_output.strip())
        else:
            logging.debug(
                "[DEBUG] (stdout of process {} : {}\n".format(processname, std_output.strip()))
            return_stdout.append(std_output.strip())

    if error_output:
        # In order to not be spammed especially by CCS which has a buggy output
        if not ("ERROR" in error_output.upper()) and not ("CRITICAL" in error_output.upper()):
            logging.debug(
                '[DEBUG] (stderr of process {} : {}\n'.format(processname, error_output.strip()))
            return_stdout.append(error_output.strip())
        else:
            logging.error(
                "[ERROR] (stderr output of process {}) : {}\n".format(
                    processname, error_output.strip()))
            return_stderr.append(error_output.strip())
            # This will show you eventual errors + will force the kernel to
            # wait the end of the process

    return "\n".join(return_stdout), "\n".join(return_stderr)


def handle_spaces_reference(referencepath: os.path.realpath,
                            tmpdir: os.path.realpath) -> os.path.realpath:
    """A preprocessing is needed because none of the PacBio tools (BLASR, CCS, kineticsTools) handles
    whitespaces the identifiers of the scaffolds in the .fasta reference file.

    In case the user's reference contains whitespaces in its identifiers, we'll automatically switch them by "_"
    in another .fasta file located in the tmpdir. Since this can break future analysis, user will be properly
    warned, and asked to provide whitespace-free references next time."""

    loaded_fasta = smsn.fasta_tools.load_fasta(referencepath)

    contains_some = False
    list_specials = [" ", ".", ",", ",", ",", "-", "&", "|", "[", "]",
                     "{", "}", "/", "\t"]

    for key in list(loaded_fasta):
        if any([character in list_specials for character in key]):
            contains_some = True

    if contains_some:
        logging.warning(
            '[WARNING] At least one ID in the .fasta file contains whitespaces or unsupported characters.')
        referencename = os.path.basename(referencepath)
        newfasta = os.path.join(tmpdir, "whitespacefree_" + referencename)
        loaded_fasta = smsn.fasta_tools.load_fasta_special(referencepath)
        smsn.fasta_tools.dict_to_fasta(loaded_fasta, newfasta)
        logging.warning(
            '[WARNING] Since the identifiers in the .fasta reference contain special characters or whitespaces, '
            'we built another reference at {}, where identifiers are whitespace-free without special characters. This '
            'might force you to make additionnal preprocessing before analyzing your results, since the scaffold '
            'names are now changed. Please rename your scaffolds to avoid future pitfalls. You can replace the '
            'problematic characters by the underscore \'_\' if needed.'.format(
                newfasta))
        return newfasta
    else:
        return referencepath


def generator_to_chunk(generator, chunksize):
    """From a generator, yields chunked lists

    >>> mygenerator = (n for n in [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])
    >>> [x for x in generator_to_chunk(mygenerator,3)]
    [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11], [12, 13, 14], [15, 16, 17]]
    >>> mygenerator = (n for n in [0,1,2,3,4])
    >>> [x for x in generator_to_chunk(mygenerator,2)]
    [[0, 1], [2, 3], [4]]
    >>> mygenerator = (n for n in [0,1,2,3,4])
    >>> [x for x in generator_to_chunk(mygenerator,5)]
    [[0, 1, 2, 3, 4]]
    """

    assert chunksize > 0
    i = 0
    chunknb = 0
    total = 0
    items = []
    while True:
        try:
            nextitem = next(generator)
            items.append(nextitem)
            i += 1

            if i >= chunksize:
                yield items
                items = []
                total += chunksize
                chunknb += 1
                i = 0
        except StopIteration:
            if len(items) > 0:
                yield items
            total += i
            break

    logging.debug('[DEBUG] Generator {} was divided in {} chunks for a total of {} items'.format(
        generator, chunknb, total
    ))
    return


def launch_smsn(args):
    """Launchs the whole pipeline. That's where every functions are called in order.
    Args are parsed from user input through the smsn_launcher.py script"""
    moviename = smsn.bam_toolbox.parse_header(args["bam"])["PU"]
    bam_header = smsn.bam_toolbox.read_header(args["bam"])

    args["moviename"] = moviename
    args["bam_header"] = bam_header

    args["pathmodel"] = transform_model_name(args["model"])
    logging.info('[INFO] Using model = {}'.format(args["pathmodel"]))

    logging.debug("[DEBUG] Ensuring that the .fasta doesn't contained empty spaces in the references")
    args["reference"] = handle_spaces_reference(args["reference"], args["tmpdir"])
    # If there is nothing to change, the reference will stay the same
    # Otherwise, it will be changed to another ref that will be in the tmpdir

    ##########################################################################
    logging.info("[INFO] Step 1 - 1 : Creating the consensus")

    if args["CCS"]:
        # This step is really expensive computationnally. We'll let the user give its own consensus if they're already
        # built
        args["CCS"] = os.path.realpath(args["CCS"])
        logging.info("[INFO] CCS {} will be used - CCS creation is skipped.".format(args["CCS"]))
    else:
        recreate_CCS(args["bam"], args["nb_proc"], tmpdir=args["tmpdir"])
        args["CCS"] = os.path.join(args["tmpdir"], moviename + ".CCS.bam")
    ##########################################################################
    logging.info("[INFO] Step 1 - 1 : Preprocessing -- Ensuring the subreads are sorted")

    sorted_subreads_path = os.path.join(args["tmpdir"], "sorted_subreads_{}.bam".format(moviename))

    smsn.bam_toolbox.sort_by_name(bamfilepath=args["bam"],
                                  newpath=sorted_subreads_path,
                                  nbcore=args["nb_proc"])

    logging.debug('[DEBUG] Sorted subreads of {} to put them in {}'.format(args["bam"],sorted_subreads_path))

    args["bam"] = sorted_subreads_path
    ##########################################################################
    logging.info(
        "[INFO] Step 2 - 1: Aligning all the CCS against the reference using BLASR")

    args["aligned_CCS"] = os.path.join(args["tmpdir"], "alignedCCS_" + os.path.basename(args["CCS"]))
    args["notaligned_CCS"] = os.path.join(
        args["tmpdir"], "NOTalignedCCS_" + moviename + ".fasta")

    align_CCS(args["CCS"],
              args["reference"],
              args["aligned_CCS"],
              args["notaligned_CCS"],
              args["tmpdir"],
              args["nb_proc"])

    logging.debug("[DEBUG] We make the list of all holes to be analyzed")

    ##########################################################################
    logging.info('[INFO] Step 2 - 2: Loading the alignment produced')
    df_aligned_CCS = smsn.bam_toolbox.detailed_get_pd_dataframe_alignment(
        args["aligned_CCS"])

    assert len(df_aligned_CCS) > 0

    # Unaligned CCS should figure unaligned in the output with our parameters
    try:  # Safety check
        assert len(df_aligned_CCS) == len(
            df_aligned_CCS['HoleID'].unique())  # Just checking that BLASR did everything right
    except AssertionError:
        logging.critical('[CRITICAL] There is not as much alignment as CCS. The rest of the analysis might be '
                         'corrupted')

    output_alignmentcsv = os.path.join(args["tmpdir"], "aligned_CCS_" + args["moviename"] + ".csv")
    logging.debug("[DEBUG] Saving alignment (.csv) at {}".format(output_alignmentcsv))
    df_aligned_CCS.to_csv(output_alignmentcsv, index=False)

    ##########################################################################
    logging.info('[INFO] Step 2 - 3: Keeping only alignments with >= {} identity'.format(args["min_identity"]))
    df_aligned_CCS = df_aligned_CCS[df_aligned_CCS["identity_score"] >= args["min_identity"]]
    df_aligned_CCS.index = [int(x) for x in df_aligned_CCS["HoleID"]]

    set_holes_to_analyse = set(df_aligned_CCS["HoleID"].unique())

    logging.debug('[DEBUG] {} holes are kept after filtering upon identity (Before filtering we had {})'.format(
        len(set_holes_to_analyse), len(df_aligned_CCS)))

    ##########################################################################
    logging.debug("[DEBUG] Initializing pandarallel")
    if args["progress_bar"]:  # Initialize pandarallel
        logging.debug('[DEBUG] Pandarallel will use a progress_bar')
        pandarallel.initialize(progress_bar=True, nb_workers=args["nb_proc"])
    else:
        logging.debug('[DEBUG] Pandarallel will not use a progress_bar')
        pandarallel.initialize(progress_bar=False, nb_workers=args["nb_proc"])

    ##########################################################################
    logging.debug("[DEBUG] Entering chunked data + Parallelism")
    logging.info("[INFO] Starting analysis")
    reader = smsn.bam_toolbox.yield_hole_by_hole(args["bam"],
                                                 restrictionlist=set_holes_to_analyse,
                                                 min_subreads=args["min_subreads"])

    chunks = generator_to_chunk(reader, chunksize=args["sizechunks"])
    nb_analyzed = compute_chunks(chunks,df_aligned_CCS,args)
    compile_results(nb_analyzed,args)

def compute_chunks(chunks,df_aligned_CCS,args):
    proto_df = []
    chunknb = 0
    nb_analyzed = 0

    for chunk in chunks:
        for samtxt in chunk:  # Preparing for pandarallel
            holeID = smsn.bam_toolbox.get_hole_id(samtxt)

            alignment_this_hole = df_aligned_CCS.loc[holeID]

            # Preparing a list of dictionnaries that can be used to feed pandarallel
            # with the arguments
            newdict = {"HoleID": holeID, "sam": args["bam_header"] + samtxt,
                       "scaffold": alignment_this_hole["scaffold"], "start": alignment_this_hole["start"],
                       "end": alignment_this_hole["end"]}

            proto_df.append(newdict.copy())

        logging.debug('[DEBUG] Analyzing chunk {} of size {}'.format(chunknb, len(chunk)))
        df = pd.DataFrame(proto_df)

        outputs = df.parallel_apply(lambda x: smsn.single_hole.analyze_singleHole(x["HoleID"],
                                                                                  x["sam"],
                                                                                  x["scaffold"],
                                                                                  x["start"],
                                                                                  x["end"],
                                                                                  args.copy()),
                                    axis=1)  # Simple trick to make it parallel
        nb_analyzed += len(outputs)

        chunk_csvpath = os.path.join(
            args["tmpdir"], "tmp_analysis_chunk_" + str(chunknb) + ".csv")
        logging.debug(
            '[DEBUG] Compiling all results for chunk n° {} into {}'.format(
                chunknb, chunk_csvpath))

        # #TODO Note somewhere safe the following line because it saved my life
        # outputs is a pd.Series
        # Each item is a DataFrame though
        # You cannot access these DataFrames using a list(outputs) or anything
        # Just with outputs.values

        if isinstance(outputs,pd.Series):
            outputs = pd.concat([x for x in outputs.values], ignore_index=True)
        else:
            outputs = pd.concat([x for x in outputs], ignore_index=True)

        outputs.to_csv(chunk_csvpath,index=False)

        logging.info("[INFO] Chunk nb {}, size = {}".format(chunknb,len(chunk)))

        chunknb += 1

        del df
        del outputs
        gc.collect()

        proto_df = []

    return nb_analyzed

def compile_results(nb_analyzed,args):
    logging.info("[INFO] Compiling all results of all chunks into a single file")

    if nb_analyzed < 0:
        logging.error('[ERROR] Not a single molecule could be analyzed. Exiting with error')
        logging.info('[INFO] SMSN ended.')
        exit(-1)

    outputs = [elt for elt in os.listdir(args["tmpdir"]) if 'tmp_analysis_chunk' in elt]

    if len(outputs) < 1:
        logging.error("[ERROR] No output .csv file")


    compilation = []
    for elt in os.listdir(args["tmpdir"]):
        if "tmp_analysis_chunk" in elt:
            tmp_path = os.path.join(args["tmpdir"], elt)
            compilation.append(pd.read_csv(tmp_path))

    try:
        # If compilation fails for some reason, be certain that we don't destroy the temporary files
        # Computation time is expensive and boring for everyone, and also it kills polar bears
        # So, give a chance to the user to get the best out of what has already
        # been computed

        output_df = pd.concat(compilation,ignore_index=True).drop_duplicates().reset_index(drop=True)
        output_df.sort_values(["HoleID", "tpl", "strand"], inplace=True)

        list_interest = ["HoleID", "scaffold", "tpl", "strand", "base", "score", "tMean", "tErr", "modelPrediction",
                         "ipdRatio", "coverage", "isboundary"]  # identificationQv, context
        if args["idQvs"]:
            list_interest.append("identificationQv")

        if args["add_context"]:
            list_interest.append("context")

        output_df[list_interest].to_csv(args["output_csv"], index=False)

        logging.info('[INFO] Result saved at {}'.format(args["output_csv"]))

        for filetmp in os.listdir(args["tmpdir"]):
            if "tmp_analysis_chunk" in filetmp:
                os.remove(os.path.join(args["tmpdir"], filetmp))

    # Warn when operation failed (ex: I/O error, server down, whatever) :
    except Exception as e:
        logging.error(
            '[ERROR] Compilation of all the output csv files failed with the following error: \n{}'.format(e))
        logging.error("[ERROR] Problem with compiling into {}".format(args["output_csv"]))
        logging.info('[INFO] For some reason, the compilation was not possible. However, it should still possible to '
                     'retrieve the tmp .csv files at {}'.format(args["tmpdir"]))

    logging.info('[INFO] SMSN ended.')


# STEP 1

def recreate_CCS(subread_file, nbproc, tmpdir) -> None:
    """Recreates the circular consensus (CCS) from a .subreads.bam (will be named [moviename].CCS.bam, in the same
    directory as the subreads. Very laxist parameters will be used in order to yield as many CCS as possible. Because we
    are not interested here in doing an assembly, it is the mapping will then determine whether a given CCS is a good
    one or not. Some files produced by PacBio will be interminently found outside of tmpdir"""

    logging.info(
        '[INFO] Recreating the circular consensus from subreads : {}, using {} CPU.'.format(subread_file, nbproc))

    subread_file = os.path.realpath(subread_file)
    subread_dir = os.path.dirname(subread_file)
    moviename = smsn.bam_toolbox.parse_header(subread_file)["PU"]
    outputbam = os.path.join(os.path.realpath(tmpdir), moviename + ".CCS.bam")

    # Here I use very laxist parameters so that I'm sure that almost every
    # hole is going to give a consensus
    cmd = "ccs --minLength 0 --minIdentity 0 --minZScore NaN --minSnr 1 --polish --minReadScore 0 --richQVs " \
          "--maxLength 5000000 --minPasses 0 --minPredictedAccuracy 0 --numThreads " + str(nbproc) + \
          " " + subread_file + " " + outputbam
    logging.debug("[DEBUG] Launching cmd = {}".format(cmd))
    std_out, std_err = call_process(cmd)

    logging.debug('[DEBUG] {}'.format(std_out))
    if std_err:
        logging.debug('[DEBUG] {}'.format(std_err))


# STEP 2

def align_CCS(CCS, reference, aligned_CCS, notaligned_CCS, tmpdir, nb_proc) -> None:
    """Aligns the PacBio consensus against a reference genome"""

    os.chdir(tmpdir)

    # # When debugging : verbosity
    # cmd = "blasr -v " + CCS + " " + reference + " --bestn 1 --hitPolicy leftmost --bam --out " + \
    #       aligned_CCS + " --clipping none --nproc " + str(nb_proc) + " --unaligned " + \
    #       os.path.join(tmpdir, notaligned_CCS)

    cmd = "blasr " + CCS + " " + reference + " --bestn 1 --hitPolicy leftmost --bam --out " + \
          aligned_CCS + " --clipping none --nproc " + str(nb_proc) + " --unaligned " + \
          os.path.join(tmpdir, notaligned_CCS)

    std_out, std_err = call_process(cmd)

    print(std_out)

    print(std_err)

    logging.debug(
        '[DEBUG] Mapping DONE for {} on {}. Output can be found at {}'.format(
            CCS, reference, aligned_CCS))
