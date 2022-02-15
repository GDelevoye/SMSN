"""Tests of the argument parsing"""

__author__ = "DELEVOYE Guillaume"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "DELEVOYE Guillaume"
__email__ = "delevoye@ens.fr"
__status__ = "Production"

import os
import smsn
import psutil

import smsn.launchers.smsn_launcher
from ..testtools.t_tools import zip_equal
from ..fixtures.files import get_resource

def test_main(get_resource,tmpdir):

    d = tmpdir
    os.chdir(tmpdir)

    nb_proc = psutil.cpu_count(logical=True)

    # this test is in two series :
        # Serie 1 : Regular usage
        # Serie 2 :  Weird usage (corrupted files, no output expected, absurd parameters etc)

    # The first serie contains three test cases :
        #  2 without the CCS
        #  1 with the CCS already provided (HTVEG)

    testnames = ["HT2","HT6","HTVEG"]

    # subreads

    htveg_subreads = os.path.join(*[get_resource, "resources_tests", "subreads_coli",
                                    "HTVEG_Coli.subreads.bam"])
    ht2_subreads = os.path.join(*[get_resource, "resources_tests", "subreads_coli",
                                    "HT2_Coli.subreads.bam"])
    ht6_subreads = os.path.join(*[get_resource, "resources_tests", "subreads_coli",
                                    "HT6_Coli.subreads.bam"])

    list_subreads = [ht2_subreads,ht6_subreads,htveg_subreads]


    # CCS

    ht2_ccs = None
    ht6_ccs = None

    htveg_ccs = os.path.join(*[get_resource, "resources_tests", "subreads_coli",
                               "CCS_HTVEG.bam"])


    list_ccs = [ht2_ccs,ht6_ccs,htveg_ccs]

    # reference

    ecoli_ref = os.path.join(*[get_resource, "resources_tests","e_coli_O157.fasta"])
    refs = [ecoli_ref]*len(testnames)

    # outputs

    outputdirs = [os.path.join(d,x) for x in testnames]
    outputs_csv = [os.path.join(x,"output.csv") for x in outputdirs]

    # min_identities

    min_identities = [0.95,0.95,0.99]

    # min_subreads

    min_subreads = [0,0,0]

    # chunksizes

    chunksizes = [1000,1000,1000]

    # nbprocs

    nb_proc = [nb_proc, nb_proc, nb_proc]

    test_cases = zip_equal(list_subreads,
                           list_ccs,
                           refs,
                           outputs_csv,
                           min_identities,
                           min_subreads,
                           chunksizes,
                           nb_proc,
                           testnames)

    for (subreadfile,
         ccsfile,
         reffile,
         outputfile,
         min_identity,
         min_subread,
         chuksize,
         nbproc,
         name) in test_cases:


        if name in ["HT2","HT6"]:
            args = ["--bam", str(subreadfile),
                    "--CCS", ccsfile,
                    "--reference", str(reffile),
                    "--output_csv", str(outputfile),
                    "--model", "SP2-C2",
                    "--min_identity", str(min_identity),
                    "--verbosity","INFO",
                    "--tmpdir", str(os.path.join(d,name)),
                    "--min_subreads", str(min_subread),
                    "--verbosity", "INFO",
                    "--nb_proc",str(nbproc),
                    "--sizechunks", str(chuksize)]


            smsn.launchers.smsn_launcher.main(testargs=args)


