import smsn
import os
from testfixtures import TempDirectory
import shutil

from .testtools.t_tools import zip_equal
from .fixtures.files import get_resource

def test_check_sorted(get_resource):
    testpath = os.path.join(*[get_resource,"resources_tests","mocked","unsorted_holes.bam"])
    assert not smsn.bam_toolbox.check_sorted(testpath)
    testpath = os.path.join(*[get_resource,"resources_tests","mocked","sorted_holes.bam"])
    assert smsn.bam_toolbox.check_sorted(testpath)


def test_sort_by_name(get_resource):
    with TempDirectory() as d:
        src = os.path.join(*[get_resource,"resources_tests","mocked","unsorted_holes.bam"])
        dst = os.path.join(d.path,"unsorted_holes.bam")

        shutil.copyfile(src,dst)

        output = os.path.join(d.path,"sorted_holes.bam")
        smsn.bam_toolbox.sort_by_name(dst,output)

        truthpath = os.path.join(*[get_resource,"resources_tests","mocked","sorted_holes.bam"])
        truth = "\n".join([x for x in smsn.bam_toolbox.read_bam(truthpath)])
        test = "\n".join([x for x in smsn.bam_toolbox.read_bam(output)])

        assert test == truth

def test_read_bam(get_resource):
    testpath = os.path.join(*[get_resource,"resources_tests","mocked","false_100_subreads_noheader.bam"])
    truthpath = os.path.join(*[get_resource,"resources_tests","mocked","false_100_subreads_noheader.txt"])

    linesbam = []
    linestxt = []
    i = 0

    for (linebam,linetxt) in zip_equal(smsn.bam_toolbox.read_bam(testpath),
                                 open(truthpath,"r")):
        linetxt = linetxt.strip()

        assert linebam.strip().__repr__() == linebam.__repr__()
        assert linebam == linetxt

        i+=1

        linesbam.append(linebam)
        linestxt.append(linetxt)

    assert "\n".join(linesbam) == "\n".join(linestxt)
    assert i == 100


def test_read_header(get_resource):
    testpath = os.path.join(*[get_resource,"resources_tests","subreads_coli","HTVEG_Coli.subreads.bam"])
    truthpath = os.path.join(*[get_resource,"resources_tests","mocked","false_header.txt"])
    truth = open(truthpath,"r").read().strip()
    got = smsn.bam_toolbox.read_header(testpath).strip()

    for truthline, gotline in zip_equal(truth.split('\n'),got.split('\n')):
        assert gotline.strip() == truthline.strip()

def test_get_hole_id(get_resource):
    testpath = os.path.join(*[get_resource,"resources_tests","mocked","false_subread.txt"])
    test = open(testpath,"r").read()
    truth = int(4260642)

    assert smsn.bam_toolbox.get_hole_id(test) == truth

def test_inmemory_asbam(get_resource):
    testpath =  os.path.join(*[get_resource,"resources_tests","mocked","false_100_lines_with_5linesheader.txt"])
    truthpath = os.path.join(*[get_resource,"resources_tests","mocked","false_100_lines_with_5linesheader.bam"])

    truth = open(truthpath,"br").read()
    samtxt = open(testpath,"r").read()
    test = smsn.bam_toolbox.inmemory_asbam(samtxt)

    assert truth == test


def test_parse_header(get_resource):
    testpath = os.path.join(*[get_resource,"resources_tests","subreads_coli","HTVEG_Coli.subreads.bam"])
    truthpath = os.path.join(*[get_resource,"resources_tests","mocked","false_header.txt"])

    # #TODO Someday The Ipd and pulsewidth codec semm wrong and should mabe be fixed
    # #FIXME
    # I did not care until now because I don't use this information myself

    # On the other hand I also wrote the following code someday :

    # logging.warning("[WARNING] Your data seems to contain raw frames of inter-pulse durations (IPDs), rather than \
    # "the usual lossy-encoding (CodecV1). See : \
    # "https://pacbiofileformats.readthedocs.io/en/3.0/BAM.html. The program will continue anyway, \
    # "and might work properly. Be carefull, however, when analyzing your data because the pipeline \
    # "was never tested without CodecV1. If your data was generated from old h5 files, \
    # "it is possible to re-use bax2bam to encode the IPDs with the CodecV1.")

    # So... Perhaps I was right at the time, but I'm not sure
    # I leave it alone for now

    truth = {'ID': '3c2d9223',
     'PL': 'PACBIO',
     'PU': 'm54063_170928_100556',
     'PM': 'SEQUEL',
     'READTYPE': 'SUBREAD',
     'Ipd:CodecV1': 'ip',
     'PulseWidth:CodecV1': 'pw',
     'BINDINGKIT': '100-862-200',
     'SEQUENCINGKIT': '101-093-700',
     'BASECALLERVERSION': '5.0.0.6236',
     'FRAMERATEHZ': 80.0}

    test = smsn.bam_toolbox.parse_header(testpath)
    test = smsn.bam_toolbox.parse_header(testpath)

    assert truth == test


def test_listholes(get_resource):
    truth = set()
    for x in [4260763,4260642,4325503]:
        truth.add(int(x))

    testpath = os.path.join(*[get_resource,"resources_tests","mocked","mocked_holes.bam"])

    assert truth == smsn.bam_toolbox.listholes(testpath)


def test_yield_hole_by_hole(get_resource):
    """Tests that each hole is yielded correctly when we iterate the .bam file"""

    # There might be some errors due to string comparisons and differences in escape encoding when
    # python reads from a file and when it gets the samtools output
    # https://stackoverflow.com/questions/23788038/python-compare-n-escape-chars-with-n

    testpath = os.path.join(*[get_resource,"resources_tests","mocked","mocked_holes.bam"])

    truthpaths = [
        "hole1.sam",
        "hole2.sam",
        "hole3.sam"
    ]

    truthpaths = [os.path.join(*[get_resource,"resources_tests","mocked",x]) for x in truthpaths]

    for testhole,truthpath in zip_equal([x for x in smsn.bam_toolbox.yield_hole_by_hole(testpath)],truthpaths):
        assert(testhole == open(truthpath,"r").read())

    # Now with a restrictionlist
    # without the headers

    restrictionlist = {4260763,4325503}
    testpath = os.path.join(*[get_resource,"resources_tests","mocked","mocked_holes.bam"])
    truthpath = os.path.join(*[get_resource,"resources_tests","mocked","restrictionlist_2holes_noheader.sam"])
    truth = open(truthpath,"r").read()
    test = "\n".join(list(smsn.bam_toolbox.yield_hole_by_hole(subreads_bamfile=testpath,
                                               restrictionlist=restrictionlist,
                                               include_header=False)))


    assert test == truth

    # With a restriction list, the header, and a minimal coverage
    restrictionlist = {4260763}
    testpath = os.path.join(*[get_resource,"resources_tests","mocked","mocked_holes.bam"])
    truth = ""
    test = "\n".join(list(smsn.bam_toolbox.yield_hole_by_hole(subreads_bamfile=testpath,
                                               restrictionlist=restrictionlist,
                                               include_header=False,
                                               min_subreads=3)))
    assert test == truth

    # When only one is returned

    # With a restriction list, the header, and a minimal coverage
    testpath = os.path.join(*[get_resource,"resources_tests","mocked","mocked_holes.bam"])
    truthpath = os.path.join(*[get_resource,"resources_tests","mocked","1hole_yieldholebyhole_3subreadsmin.sam"])
    test = "\n".join(list(smsn.bam_toolbox.yield_hole_by_hole(subreads_bamfile=testpath,
                                               restrictionlist=None,
                                               include_header=True,
                                               min_subreads=3)))

    truth = open(truthpath,"r").read()
    assert test == truth

def test_detailed_get_pd_dataframe_alignment(get_resource):
    """Ensures particularly that the integers are conserved"""
    bamfile = os.path.join(*[get_resource,"resources_tests","mocked","mocked_alignment.bam"])

    df = smsn.bam_toolbox.detailed_get_pd_dataframe_alignment(bamfile)
    assert len(df) == 3

    movienames = ["m54063_170928_100556",
                  "m54063_170928_100556",
                  "m54063_170928_100556"]

    holeIDs = [4260642,
               4260763,
               4325503]

    scaffolds = ["NZ_AVCD01000005.1"]*3
    

    starts = [400918,
              3263306,
              1998406]

    qvs = [254] * 3

    second_seq = 'GTGGGCGCCCTCGATAGCGAAGTAAAAAGTCTGCACGACATAACGTGCGTATGCGATTATTGGCGATACCAGTCGCTTTAACCGGTTGCAAGGAACGGAT'\
    'TCGTAAATCTGAAGCCAACAGCCGGGAATACAGGGTCTGACGCTGAATATTGCGGCGAACTACGGTGGACGTTGCGGATATAGTCCAGGGAGTCAGGCAACTGGCTGAAAAGG'\
    'TGCAGAGGAAACCGCAAACCAGATCAGATAGCATGAAGAGATGCTAAACCAGCCATGTCTGTATGCATGAAC'
    'TGCAGAGGAAACCGCAAACCAGATCAGATAGCATGAAGAGATGCTAAACCAGCCATGTCTGTATGCATGAAC'

    clipped = [3,9,2]

    flags = [0,16,16]

    matching = [328,267,342]


    for i in range(0,3):
        for valuelist,column in zip([holeIDs, starts, qvs, clipped, flags, matching],
                                    ["HoleID","start","mapQV","clipped_bases","alignment_flag","matching_bases"]):
            assert df.iloc[i][column] == valuelist[i]

    assert (not (any([x!=y for (x,y) in zip(movienames,df["moviename"])])))
    assert (not (any([x!=y for (x,y) in zip(scaffolds,df["scaffold"])])))
    assert df.iloc[1]["sequence"] == second_seq

    bamfile = os.path.join(*[get_resource,"resources_tests","mocked","mocked_holes.bam"])
    df = smsn.bam_toolbox.detailed_get_pd_dataframe_alignment(bamfile)
    assert len(df) == 11
