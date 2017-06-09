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

from tool.kallisto_indexer import kallistoIndexerTool

def test_kallisto_indexer():
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    ki= kallistoIndexerTool()
    
    fasta_file = resource_path+"chr7_gene_Gene.fasta"
    ki.run ([fasta_file],{}, [resource_path+"chr7_geneForRNAseq.fasta.idx"] )
    
    print(__file__)
    
