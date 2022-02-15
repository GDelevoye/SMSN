#!/usr/bin/env python

__author__ = "DELEVOYE Guillaume"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "DELEVOYE Guillaume"
__email__ = "delevoye@ens.fr"
__status__ = "Production"

from .testtools.t_tools import zip_equal
from .fixtures.files import get_resource

from testfixtures import TempDirectory
import shutil
import smsn
import os


# def test_getAbsPath():  # Already a docstring
#     pass
# def test_transform_model_name():  # Already a docstring
#     pass
# def test_call_process():  # Already a docstring
#     pass

def clean_from_PG(sam):
    """The @PG header line put here by CCS contains the directory
    This directory changes from one test to another (a .tmp dir)
    This is why I have to remove it. Hence this function
    """
    toreturn = []
    for line in sam.split('\n'):
        if line[0:3] == "@PG":
            continue
        else:
            toreturn.append(line.strip())
    return "\n".join(toreturn)


def test_handle_spaces_references(get_resource):

    with TempDirectory() as d:
        inputtest = os.path.join(*[get_resource,"resources_tests","mocked","misformated.fa"])
        truthpath = os.path.join(*[get_resource,"resources_tests","mocked","reformated.fa"])

        smsn.pipeline.handle_spaces_reference(inputtest,d.path)

        newfasta = os.path.join(d.path, "whitespacefree_" + "misformated.fa")

        assert open(newfasta,"r").read() == open(truthpath,"r").read()
        # Check when nothing is wrong with the reference
        renew_fasta = smsn.pipeline.handle_spaces_reference(newfasta,d)
        assert open(newfasta,"r").read() == open(renew_fasta,"r").read()

        inputtest = os.path.join(*[get_resource,"resources_tests","mocked","reformated.fa"])

def test_recreate_CCS(get_resource):

    with TempDirectory() as d:
        inputtest = os.path.join(*[get_resource,"resources_tests","mocked","1hole_yieldholebyhole_3subreadsmin.bam"])
        truthpath = os.path.join(*[get_resource,
                                   "resources_tests",
                                   "mocked",
                                   "holeID_4260642_m54063_170928_100556.CCS.sam"])

        os.chdir(d.path) # The function is sensitive to the directory where it is called from

        smsn.pipeline.recreate_CCS(inputtest,
                                   nbproc=1,
                                   tmpdir=d.path)

        # Check that all files are generated in the tmpdir, CCS report included
        assert sorted(os.listdir(d.path)) == sorted(['m54063_170928_100556.CCS.bam.pbi',
                                      'ccs_report.txt',
                                      'm54063_170928_100556.CCS.bam'])

        moviename = smsn.bam_toolbox.parse_header(inputtest)["PU"]
        outputbam = os.path.join(os.path.realpath(d.path), moviename + ".CCS.bam")

        header_out = smsn.bam_toolbox.read_header(outputbam)
        sam_out = "\n".join([x for x in smsn.bam_toolbox.read_bam(outputbam)])
        sam_out = header_out.strip() + "\n" + sam_out


        truth = clean_from_PG(open(truthpath,"r").read())
        test = clean_from_PG(sam_out)

        assert (truth == test)



def test_align_CCS(get_resource):
    with TempDirectory() as d:
        samtest = os.path.join(*[get_resource,"resources_tests","mocked",
                                   "test_ccs_toalign.sam"])

        fastaref = os.path.join(*[get_resource,"resources_tests","mocked",
                                  "test_falseref_toalign.fa"])

        truthpath = os.path.join(*[get_resource,"resources_tests","mocked",
                                  "example_aligned_CCS.sam"])

        inputtest_bamstring = smsn.bam_toolbox.inmemory_asbam(samtest)

        # Write it here as a .bam file (input is a sam to facilitate debugging)
        samtxt = open(samtest,"r").read()
        bamtxt = smsn.bam_toolbox.inmemory_asbam(samtxt)
        bamfile = os.path.join(d.path,"input.bam")

        with open(bamfile,"wb") as bf:
            bf.write(bamtxt)

        notaligned_path = os.path.join(d.path,"notaligned.ccs.fa")
        output_alignmentfile = os.path.join(d.path,"aligned_CCS.bam")

        smsn.pipeline.align_CCS(bamfile,
                                fastaref,
                                output_alignmentfile,
                                notaligned_path,
                                d.path,
                                1)

        assert sorted(os.listdir(d.path)) == sorted(
            ['aligned_CCS.bam', 'input.bam', 'notaligned.ccs.fa'])

        truth = open(truthpath,"r").read()
        test = "\n".join([x for x in smsn.bam_toolbox.read_bam(output_alignmentfile)])
        # Careful : The header is not read by the function, so we should add it !
        test = smsn.bam_toolbox.read_header(output_alignmentfile).strip() + "\n" + test

        test = clean_from_PG(test).strip()
        truth = clean_from_PG(truth).strip()

        for line1,line2 in zip_equal(test.split('\n'),truth.split('\n')):
            assert line1 == line2

