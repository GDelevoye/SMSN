#!/usr/bin/env python

__author__ = "DELEVOYE Guillaume"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "DELEVOYE Guillaume"
__email__ = "delevoye@ens.fr"
__status__ = "Developpment"

"""Get the in-sillico control values from kineticsTools"""

from kineticsTools.ipdModel import IpdModel
from pkg_resources import Requirement, resource_filename
import os

def _getAbsPath(fname):
    return resource_filename(Requirement.parse('smsn'), 'smsn/%s' % fname)

def transform_model_name(modelname):
    resources_dir = _getAbsPath("/resources/")
    modifiedmodelname = modelname + ".npz.gz"
    modelname = os.path.join(resources_dir, modifiedmodelname)
    return modelname

class Contig():
    def __init__(self, id, seq):
        self.id = id
        self.ID = id
        self.sequence = seq
        self.alignmentID = id

class Str2IPD():
    """This is the object I'd recommend to use to predict ipds according to an in-sillico model.
    The in-sillico model predicts IPDs on any nucleotide using the N-5/N+4 nucleotide context

    When it is not possible to get the snipet of a given sequence (example: Trying to predict ipds outside of the
    sequence, or near the end of a sequence), note that a padding with Adenines will be used by PacBio's program.

    Usage:
    >>> hacked_model = Str2IPD("ATCGATGCGGATTGCGTTGT",model="SP2-C2",indexing=1)
    >>> hacked_model.predict(position=3,strand=0)
    1.1981472
    >>>
    """

    def __init__(self, sequence, name="NO_ID", model="SP2-C2",indexing=1,secure_check=False):
        if secure_check:
            assert not any([character not in ["A","T","C","G","N","a","t","c","g","n"] for character in sequence])
            assert indexing in [0,1]

        self.indexing = indexing
        assert sequence != None and sequence != "" and len(sequence) > 0
        self.character_seq = sequence

        self.sequence = [Contig(name, sequence)]

        self.model = IpdModel(fastaRecords=self.sequence, modelFile=transform_model_name(model))
        self.predictfunc = self.model.predictIpdFuncModel(refId=name)

    def predict(self, position, strand=0):
        assert isinstance(position,int)
        assert strand in [0,1]

        if self.indexing == 1:
            assert position >= 1
            assert position+1 <= len(self.character_seq)
            return self.predictfunc(position+1, strand)
        elif self.indexing == 0:
            assert position >= 0
            assert position <= len(self.character_seq)
            return self.predictfunc(position, strand)



