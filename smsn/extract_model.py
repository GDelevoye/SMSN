#!/usr/bin/env python

__author__ = "DELEVOYE Guillaume"
__credits__ = ["BAHIN Mathieu", "MEYER Eric"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "DELEVOYE Guillaume"
__email__ = "delevoye@ens.fr"
__status__ = "Developpment"

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



