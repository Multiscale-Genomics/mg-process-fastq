#!/bin/sh
#!usr/bin/python

"""
.. Copyright 2017 EMBL-European Bioinformatics Institute
 
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at 

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

import pytest
import random
import os.path

from tool.kallisto_quant import kallistoQuantificationTool

def test_kallisto_quant():
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    
    kqft = kallistoQuantificationTool()
    fastq1 = resource_path+"kallisto.Human.ERR030872_1.fastq"
    fastq2 = resource_path+"kallisto.Human.ERR030872_2.fastq"
    kqft.run ([resource_path+"kallisto.Human.GRCh38.fasta.idx", fastq1, fastq2],{}, )
    
    print(__file__)
    
