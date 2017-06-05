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

import os, shutil, shlex, subprocess

try:
    from pycompss.api.parameter import FILE_IN, FILE_OUT, FILE_INOUT
    from pycompss.api.task import task
except ImportError :
    print ("[Warning] Cannot import \"pycompss\" API packages.")
    print ("          Using mock decorators.")
    
    from dummy_pycompss import *

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

from . import common

# ------------------------------------------------------------------------------

class kallistoIndexerTool(Tool):
    """
    Tool for running indexers over a genome FASTA file
    """
    
    def __init__(self):
        """
        Init function
        """
        print ("Kallisto Indexer")
    
    @task(cdna_file_loc=FILE_IN, cdna_idx_file=FILE_INOUT)
    def kallisto_indexer(self, cdna_file_loc, cdna_idx_file):
        """
        Kallisto Indexer
        
        Parameters
        ----------
        file_loc : str
            Location of the cDNA FASTA file for a genome
        idx_loc : str
            Location of the output index file
        """
        
        command_line = 'kallisto index -i ' + cdna_idx_file + ' ' + cdna_file_loc
        
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()
        
        return True
    
    
    def run(self, input_files, metadata):
        """
        Tool for generating assembly aligner index files for use with Kallisto
        
        Parameters
        ----------
        input_files : list
            FASTA file location will all the cDNA sequences for a given genome
        metadata : list
        
        Returns
        -------
        array : list
            First element is a list of the index files. Second element is a
            list of the matching metadata
        """
        
        file_name = input_files[0]
        output_file = input_files[1]
        
        # input and output share most metadata
        output_metadata = {}
        
        # handle error
        if not self.kallisto_indexer(input_files[0], output_file):
            output_metadata.set_exception(
                Exception(
                    "kallisto_indexer: Could not process files {}, {}.".format(*input_files)))
        output_file = None
        return ([output_file], [output_metadata])

# ------------------------------------------------------------------------------
