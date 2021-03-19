import os, logging
import multiprocessing as mp
import pandas as pd
from Joblib import delayed , Parallel


def true_smrt(experiment, bamfilePath, alignmentCSV, fastafile, outputdir, tmpdir, genome, csv_out, nbproc=1,
              sizechunks=None):
    bamfilePath = os.path.realpath(bamfilePath)
    subreads_bamfile = bamfilePath
    alignmentCSV = os.path.realpath(alignmentCSV)
    outputdir = os.path.realpath(outputdir)
    tmpdir = os.path.realpath(tmpdir)
    csv_out = os.path.realpath(csv_out)

    # Makes sure that all the required directories exist
    os.system('mkdir -p ' + str(outputdir))
    os.system('mkdir -p ' + str(tmpdir))
    csv_outdir = os.path.dirname(os.path.realpath(csv_out))
    os.system('mkdir -p ' + str(csv_outdir))

    args = [experiment, bamfilePath, alignmentCSV, fastafile, outputdir, tmpdir, genome, csv_out, nbproc]

    logging.info('[INFO] (true_smrt) Launching true_smrt() with arguments: {}'.format(args))

    ##########################################
    # Launching the writer process ...       #
    ##########################################


    logging.debug('[DEBUG] (true_smrt) Creating the manager queue')
    manager = mp.Manager()
    queue = manager.Queue()

    logging.debug('[DEBUG] (true_smrt) Launching the writer_process')
    writer_process = mp.Process(target=result_writer, args=(queue, (os.path.realpath(os.path.join(outputdir, csv_out)))))
    writer_process.start()  # This will write every worker's result that has been pushed in the queue

    logging.debug('[DEBUG] Waiting for the workers to be done')

    ##########################################
    # Launching analysis of holes in parallel#
    ##########################################
    Parallel(n_jobs=nbproc - 1)(delayed(worker_perform_one_hole_analysis)(**argdict) for argdict in
                                YAD_perform_analysis_one_hole(experiment, alignmentCSV, bamfilePath, fastafile, outputdir,
                                                              queue, tmpdir, genome))

    ################################################
    # END OF THE ANALYSIS / WAITING FOR THE WRITER #
    ################################################

    logging.info('[INFO] (true_smrt) Analysis done. Waiting for the writer_process to be done')
    queue.put('poison')  # Poisoning the end of the queue
    writer_process.join()  # Waiting writer_process to reach the poison
    logging.info('[INFO] (true_smrt) Analysis DONE and saved')  # You're welcome


def YAD_perform_analysis_one_hole(experiment, alignmentCSV, bamfilePath, fastafile, outputdir, queueresults, tmpdir,
                                  genome):
    """
    Yields arguments for the function 'worker_perform_one_hole (multiprocessing implementation)
    Assumes that the bamfile is correctly sorted
    """

    logging.debug('[DEBUG] (YAD_perform_analysis_one_hole) Generating arguments for the analysing function')
    pdAlignment = pd.read_csv(os.path.realpath(alignmentCSV))

    # Getting the header only once for all
    # Inmemory_asbam translates .sam to .bam, reshape_header removes all the header's line that will cause errors, and readheader justs reads the header as a .sam file
    # ipdSummary cannot work if the header contains reference lines that do not exist in the reference (which is +/- logical actually) so this is why I use reshape_header, that will remove all un-necessary lines

    header = st.read_header(bamfilePath)

    pdAlignment["Best_alignment_this_hole"] = pdAlignment.groupby(['HoleID'])["identity_score"].transform(max)

    # Getting the name of the bamfile
    basename = os.path.basename(bamfilePath)
    if len(basename.split('.')) > 1:  # if it has at least one extension
        name_without_extension = '.'.join(x for x in basename.split('.')[:-1])  # We remove the last one

    for sam in yield_hole_by_hole(bamfilePath):
        dico = {}
        dico[
            'samseq'] = sam  # It is already cleaned from any previously existing alignment or any anomaly in the header linked to a previous alignment
        dico['workdir'] = tmpdir
        dico["genome"] = genome
        logging.debug(
            '[DEBUG] (YAD_perform_analysis_one_hole) yielding for holeID {} and genome {}'.format(int(get_hole_id(sam)),
                                                                                                  genome))

        logging.debug('[DEBUG] (YAD_perform_analysis_one_hole) Aligment file = {}'.format(alignmentCSV))
        # logging.debug("[DEBUG] (YAD_perform_analysis_one_hole) pdAlignment = {}".format(pdAlignment))
        try:
            (real_start, real_end, scaffold) = get_best_alignment(pdAlignment, int(get_hole_id(sam)), genome)
        except:
            logging.error("[ERROR] (YAD) ERROR for hole {} --> Could not retrieve the alignment within {}".format(
                int(get_hole_id(sam))), alignmentCSV)

        dico['real_start'] = real_start
        dico['real_end'] = real_end
        dico['scaffold'] = scaffold
        dico['fastafile'] = os.path.realpath(fastafile)
        dico['QueueResults'] = queueresults
        dico["experiment"] = experiment

        logging.debug(
            '[DEBUG] (YAD_perform_analysis_one_hole) Caring about : Hole =  {} Name =  {}'.format(int(get_hole_id(sam)),
                                                                                                  name_without_extension))

        yield dico


def result_writer(queue, fileOut):
    """ Will be launched as an independant process, writing every that will work until it reads 'poison' on the queue. After it's done, he will close the files"""

    os.system('mkdir -p ' + str(os.path.dirname(fileOut)))
    list_df = []

    logging.debug('[DEBUG] (result_writer) result_writer is ready to write (destination = {})'.format(fileOut))

    while True:  # Continuously listens to the queue
        if not queue.empty():  # most of the time it will be empty
            to_write = queue.get()  # This is a .csv string

            if to_write == 'poison':  # If the queue is poisoned
                logging.debug('[DEBUG] (result_writer) Last line received, writer will end')
                logging.debug(
                    '[DEBUG] (result_writer) WAIT FOR THE WRITER (don\'t interrupt otherwise all computed value will be lost)')
                final_result = pd.concat(list_df, ignore_index=True)
                with open(os.path.realpath(fileOut), "w") as myfile:
                    final_result.to_csv(myfile, sep=';', index=False)
                logging.debug('[DEBUG] (result_writer) Written in {}'.format(fileOut))

                break

            else:

                df = to_write  # It's a pd.DataFrame
                list_df.append(df.copy())

                logging.debug('[DEBUG] (result_writer) Cleaning hole {}'.format(str(os.path.dirname(to_write))))
                # os.system('rm -rf '+str(os.path.dirname(to_write)))
                logging.debug('[DEBUG] (result_writer) Hole cleaned --> {}'.format(os.path.dirname(to_write)))

                logging.debug('[DEBUG] (result_writer) Results collected at : {}'.format(to_write))

                del df
                del to_write

    logging.debug('[DEBUG] (result_writer) Result saved here = {}'.format(os.path.realpath(fileOut)))

    return


def yield_hole_by_hole(subreads_bamfile):
    """ Loads the whole .bam in memory, sorts it inplace by HoleID and yields HoleID by HoleID"""
    samseq = []
    header = st.read_header(os.path.realpath(subreads_bamfile))

    for line in st.read_bam(os.path.realpath(subreads_bamfile)):
        samseq.append(line)
    samseq.sort(key=lambda x: int(x.split()[0].split('/')[1]))  # Sorting by holeID

    # Yielding the samseq hole by hole
    current_id = 0
    current_sam = []
    first = True

    for line in samseq:
        if first:
            first = False
            id = int(line.split()[0].split('/')[1])
            current_id = int(line.split()[0].split('/')[1])
            current_sam.append(line)
        else:
            current_id = int(line.split()[0].split('/')[1])
            if current_id == id:
                current_sam.append(line)
            else:
                # print (st.reshape_header(str(header)+"".join(current_sam)))
                yield (st.reshape_header(str(header) + "".join(current_sam)))
                del current_sam
                current_sam = []
                current_sam.append(line)
                id = int(line.split()[0].split('/')[1])
                current_id = int(line.split()[0].split('/')[1])


def get_hole_id(samseq):
    try:
        holeID = "".join([x for x in samseq.split('\n') if '@' not in x]).split()[0].split('/')[1]
    except:
        holeID = -1
    return holeID


def compute_chunk_infos(real_start, real_end, fasta_this_scaffold):
    """Creates a chunk to build a false ref at +100 / -100 of the begin/end where the CCS mapped"""
    if real_start - 100 < 0:
        chunk_start = 0
    else:
        chunk_start = real_start - 100
    if real_end + 100 > len(fasta_this_scaffold):
        chunk_end = len(fasta_this_scaffold)
    else:
        chunk_end = real_end + 100

    chunk_size = chunk_end - chunk_start
    offset = real_start - 100
    sequence = fasta_this_scaffold[chunk_start:chunk_end]

    return (chunk_start, chunk_end, chunk_size, offset, sequence)


def worker_perform_one_hole_analysis(experiment, samseq, workdir, QueueResults, real_start, real_end, scaffold, genome,
                                     fastafile):
    fasta = load_fasta_special(os.path.realpath(fastafile))
    HERE = os.getcwd()
    os.chdir(workdir)

    holeID = get_hole_id(samseq)
    holeNumber = holeID  # just because I tend to use one alias or another after copying/cutting/pasting codes...

    os.system('mkdir -p ' + str(holeID))
    os.chdir(str(holeID))

    filout = open(str(holeNumber) + '.bam', 'wb')  # Let's drop the .bam file here
    filout.write(st.inmemory_asbam(samseq))  # We need a clean unaligned .bam
    filout.close()

    (chunk_start, chunk_end, chunk_size, offset, sequence) = compute_chunk_infos(int(real_start), int(real_end),
                                                                                 str(fasta[scaffold]))

    # Taking only a sub-reference (+ 100 nt / -100 nt) to re-align all the subreads precisely where the CCS mapped

    dict_to_fasta({str(scaffold): str(sequence)}, './chunked_ref.fasta', specify_HoleID=True)

    # Map all the subreads restrictively on the scaffold
    cmd = 'blasr ' + str(holeNumber) + '.bam ./chunked_ref.fasta ' + \
          ' --useccs --bestn 1 --clipping none --bam --out aligned_on_restrictedscaffold_' + str(holeNumber) + \
          '.bam --unaligned ' + str(holeNumber) + '.unaligned.fasta'

    aligned_on_restrictedscaffold = os.getcwd() + '/' + 'aligned_on_restrictedscaffold_' + str(holeNumber) + '.bam'
    logging.debug('Executing {}'.format(cmd))
    os.system(cmd)

    # Indexing the mapped .bam
    logging.debug('[DEBUG] (worker_perform_analysis_one_hole) Generating index for the aligned file')
    os.system('pbindex aligned_on_restrictedscaffold_' + str(holeNumber) + '.bam')
    os.system('samtools faidx ./chunked_ref.fasta')

    # Perform the analysis itself
    # We don't switch the mode of ipdSummary with the hack for it has already been made before in the 'true_smrt' function
    results = get_ipdSummary_details('./' + str(holeNumber) + '.bam', './chunked_ref.fasta', workdir=None,
                                     threshold_coverage=25, single_hole=str(holeNumber), nbcore=1,
                                     chunk_size=chunk_size)

    results = pd.DataFrame(results)
    try:
        results['tpl'] = results['tpl'] + 1 + offset  # Returning the results into the proper coordinates
    except:
        results['tpl'] = 'NaN'

    results["HoleID"] = int(holeID)
    results["scaffold"] = str(scaffold)
    results["genome"] = genome
    results["experiment"] = experiment

    # Write result
    logging.debug(
        '[DEBUG] (worker_perform_analysis_one_hole) Queuing the result with the writer on hole {}'.format(holeNumber))
    # See https://stackoverflow.com/questions/10880813/typeerror-sequence-item-0-expected-string-int-found for more details

    with open('./methylation_analysis.csv', 'w') as methylationfile:
        results.to_csv(methylationfile, sep=';', index=False)
    # filout = open('./methylation_analysis.csv', 'r')

    QueueResults.put(os.path.realpath(
        os.getcwd() + "/methylation_analysis.csv"))  # The writer process will automatically write anything pushed in the queue

    # os.chdir('./../')
    os.system('rm -rf ' + str(holeID))

    os.chdir(HERE)
