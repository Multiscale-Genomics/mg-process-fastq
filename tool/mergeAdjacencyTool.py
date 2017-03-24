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

import os

try:
    from pycompss.api.parameter import FILE_IN, FILE_OUT
    from pycompss.api.task import task
except ImportError :
    print "[Warning] Cannot import \"pycompss\" API packages."
    print "          Using mock decorators."
    
    from dummy_pycompss import *

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

from common import common

# ------------------------------------------------------------------------------

class mergeAdjacencyTool(Tool):
    """
    Tool for running indexers over a genome FASTA file
    """
    
    def __init__(self):
        """
        Init function
        """
        print "GEM Indexer"
    
    def merge_adjacency_files(self, adjlist_files = []):
        """
        Merge adjacency files

        Used to merge multiple datasets into a single adjacency file

        Parameters
        ----------
        """

        from fastq2adjacency import fastq2adjacency
        
        f2a = fastq2adjacency()
        
        f2a.load_hic_read_data(adj_list[0])
        new_hic_data = f2a.hic_data
        
        for i in range(1,len(adj_list)):
            f2a = fastq2adjacency()
            f2a.load_hic_read_data(adj_list[i])
        
            new_hic_data += f2a.hic_data
        
        f2a = fastq2adjacency()
        f2a.hic_data = new_hic_data

        save_hic_data()


    def run(self, input_files, metadata):
        """
        Tool for generating assembly aligner index files for use with the GEM
        indexer
        
        Parameters
        ----------
        input_files : list
            List with a single str element with the location of the genome
            assembly FASTA file
        metadata : list
        
        Returns
        -------
        array : list
            First element is a list of the index files. Second element is a
            list of the matching metadata
        """
        
        file_name = input_files[0].split('/')
        output_file = '/'.join(file_name[-1].replace('.fa', ''))
        
        # input and output share most metadata
        output_metadata = dict(
            data_type=metadata[0]["data_type"],
            file_type=metadata[0]["file_type"],
            meta_data=metadata[0]["meta_data"])
        
        # handle error
        if not self.gem_indexer(input_files[0], output_file):
            output_metadata.set_exception(
                Exception(
                    "gem_indexer: Could not process files {}, {}.".format(*input_files)))
        output_file = None
        return ([output_file], [output_metadata])

# ------------------------------------------------------------------------------
