import smsn
import os
from testfixtures import TempDirectory

from .testtools.t_tools import zip_equal
from .fixtures.files import get_resource


def test_parse_gff(get_resource):  # Not really a .fasta tool but there's no better place than here
    """According to https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md"""
    testnames = ["emptygff.gff3",  # Necessary because pacbio can output empty files
                 "standardgff.gff3",
                 "pacbiogff.gff3"]

    testcases = [os.path.join(*[get_resource, "resources_tests", "mocked", x]) for x in testnames]
    minimum_list_interest = ["scaffold", "source", "type", "start",
                             "end", "score", "strand", "phase", "attributes"]

    for testname, testfile in zip_equal(testnames, testcases):
        if testname == "emptygff.gff3":
            output = smsn.fasta_tools.parse_gff(testfile)
            assert len(output) == 0
            assert list(output) == minimum_list_interest

        elif testname == "standardgff.gff3":
            output = smsn.fasta_tools.parse_gff(testfile, parse_attributes=False)
            assert len(output) == 23
            assert list(output) == minimum_list_interest

        elif testname == "pacbiogff.gff3":
            output = smsn.fasta_tools.parse_gff(testfile, parse_attributes=True)
            assert [x for x in output["identificationQv"]] == ["20", "22"]
            assert len(output) == 2
            assert list(output) == (minimum_list_interest + sorted(["identificationQv", "context", "modQV"]))

            true_sequences = ["AAACCAACTCCTGCGCTGTCGGAGA",
                              "TCTCCGACAGCGCAGGAGTTGGTTT"]

            assert [x == y for (x, y) in zip_equal(output["context"], true_sequences)]


def test_load_fasta(get_resource):
    truthdict = {"sequence1": "ATCGCGGATTCGATTCTGAGCTGATCGATGCCTGATAGCTGATATTCGATCTTCTGAGTCGA",
                 "sequence 2": "ATGC",
                 "sequence 4": "GTAGCTAggatcga"}

    testpath = os.path.join(*[get_resource, "resources_tests", "mocked", "testfasta.fa"])
    test = smsn.fasta_tools.load_fasta(testpath)

    assert test == truthdict


def test_load_fasta_special(get_resource):
    truthdict = {"sequence1": "ATCGCGGATTCGATTCTGAGCTGATCGATGCCTGATAGCTGATATTCGATCTTCTGAGTCGA",
                 "sequence_2": "ATGC",
                 "sequence_4": "GTAGCTAGGATCGA"}

    testpath = os.path.join(*[get_resource, "resources_tests", "mocked", "testfasta.fa"])
    test = smsn.fasta_tools.load_fasta_special(testpath)

    assert test == truthdict


# def test_reverse_cpl(): # Docstring already exists
#   pass

# def test_get_snipet(): # Docstring already exists
#     pass


def test_dict_to_fasta(get_resource):
    truthdict = {"sequence1": 'ATCGCGGATTCGATTCTGAGCTGATCGATGCCTGATAGCTGATATTCGATCTTCTGAGTCGA',
                 "sequence 2": "ATGC",
                 "sequence 4": "GTAGCTAggatcga"}

    with TempDirectory() as d:
        testpath = os.path.join(d.path, "reformed_fasta_from_dict.fa")
        smsn.fasta_tools.dict_to_fasta(truthdict, testpath)
        assert smsn.fasta_tools.load_fasta(testpath) == truthdict


def test_compute_chunk_infos(get_resource):
    fastapath = os.path.join(*[get_resource, "resources_tests", "mocked", "backtofasta.fa"])
    scaffold_sequence = smsn.fasta_tools.load_fasta(fastapath)["sequence1"]
    # same as :
    # scaffold_sequence = "ATCGCGGATTCGATTCTGAGCTGATCGATGCCTGATAGCTGATATTCGATCTTCTGAGTCGA"
    assert len(scaffold_sequence) == 62

    assert scaffold_sequence[0] == "A"
    assert scaffold_sequence[-1] == "A"

    # I'm gonna pretend that I want the chunk at +2 / -2 nt where the CCS mapped
    distance = 2
    CCS_start = 6  # First position in the reference where a nucleotide is aligned
    CCS_end = 10  # Last position in the reference where a nucleotide is aligned
    assert scaffold_sequence[CCS_start:CCS_end + 1] == "GATTC"

    chk_start, chk_end, chk_size, offset, chk_squence = \
        smsn.fasta_tools.compute_chunk_infos(real_start=CCS_start,
                                             real_end=CCS_end + 1,
                                             fasta_this_scaffold=scaffold_sequence,
                                             distance=distance)

    # /!\ #WARNING Do not forget that the function was called with CCS_end+1, not CCS_end

    # by definition :
    shouldreturn = scaffold_sequence[CCS_start-distance:CCS_end + 1 + distance]
    assert shouldreturn == "CGGATTCGA"  # What sequence the function should return

    assert chk_squence == shouldreturn
    assert chk_squence == scaffold_sequence[chk_start:chk_end] # This is how the user must use it

    # Test overflow on left

    distance = 10
    CCS_start = 6  # included
    CCS_end = 10  # included

    assert scaffold_sequence[CCS_start:CCS_end + 1] == "GATTC"

    chk_start, chk_end, chk_size, offset, sequence = \
        smsn.fasta_tools.compute_chunk_infos(real_start=CCS_start,
                                             real_end=CCS_end + 1,
                                             fasta_this_scaffold=scaffold_sequence,
                                             distance=distance)

    truth_sequence = "ATCGCGGATTCGATTCTGAGC"
    assert sequence == truth_sequence
    assert sequence == scaffold_sequence[chk_start:chk_end]

    # Test overflow on right

    distance = 10
    CCS_start = 55
    CCS_end = 62

    assert scaffold_sequence[CCS_start:CCS_end + 1] == "GAGTCGA"

    chk_start, chk_end, chk_size, offset, sequence = \
        smsn.fasta_tools.compute_chunk_infos(real_start=CCS_start,
                                             real_end=CCS_end + 1,
                                             fasta_this_scaffold=scaffold_sequence,
                                             distance=distance)

    truth_sequence = "TCGATCTTCTGAGTCGA"
    assert sequence == truth_sequence
    assert sequence == scaffold_sequence[chk_start:chk_end]

    # Test overflow on both sides

    distance = 100
    CCS_start = 55
    CCS_end = 62

    assert scaffold_sequence[CCS_start:CCS_end + 1] == "GAGTCGA"

    chk_start, chk_end, chk_size, offset, sequence = \
        smsn.fasta_tools.compute_chunk_infos(real_start=CCS_start,
                                             real_end=CCS_end + 1,
                                             fasta_this_scaffold=scaffold_sequence,
                                             distance=distance)

    truth_sequence = "ATCGCGGATTCGATTCTGAGCTGATCGATGCCTGATAGCTGATATTCGATCTTCTGAGTCGA"
    assert sequence == truth_sequence
    assert sequence == scaffold_sequence[chk_start:chk_end]
