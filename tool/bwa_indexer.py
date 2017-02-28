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

from ..common import common

# ------------------------------------------------------------------------------

class bwaIndexerTool(Tool):
    """
    Tool for running indexers over a genome FASTA file
    """
    
    def __init__(self):
        """
        Init function
        """
        print "BWA Indexer"
    
    @task(file_loc=FILE_IN, amb_loc=FILE_OUT, ann_loc=FILE_OUT, bwt_loc=FILE_OUT, pac_loc=FILE_OUT, sa_loc=FILE_OUT)
    def bwa_indexer(self, file_loc, amb_loc, ann_loc, bwt_loc, pac_loc, sa_loc):
        """
        BWA Indexer
        
        Parameters
        ----------
        file_loc : str
            Location of the genome assebly FASTA file
        amb_loc : str
            Location of the output file
        ann_loc : str
            Location of the output file
        bwt_loc : str
            Location of the output file
        pac_loc : str
            Location of the output file
        sa_loc : str
            Location of the output file
        """
        cf = common()
        amb_loc, ann_loc, bwt_loc, pac_loc, sa_loc = cf.bwa_index_genome(file_loc)
        return True
    
    def run(self, input_files, metadata):
        """
        Function to run the BWA over a genome assembly FASTA file to generate
        the matching index for use with the aligner
        
        Parameters
        ----------
        input_files : list
            List containing the location of the genome assembly FASTA file
        meta_data : list
        
        Returns
        -------
        list
            amb_loc : str
                Location of the output file
            ann_loc : str
                Location of the output file
            bwt_loc : str
                Location of the output file
            pac_loc : str
                Location of the output file
            sa_loc : str
                Location of the output file
        """
        genome_file = input_files[0]
        
        amb_name = genome_file.split("/")
        amb_name[-1].replace('.fa', '.amb')
        amb_loc = '/'.join(amb_name)
        
        ann_name = genome_file.split("/")
        ann_name[-1].replace('.fa', '.ann')
        ann_loc = '/'.join(ann_name)
        
        bwt_name = genome_file.split("/")
        bwt_name[-1].replace('.fa', '.bwt')
        bwt_loc = '/'.join(bwt_name)
        
        pac_name = genome_file.split("/")
        pac_name[-1].replace('.fa', '.pac')
        pac_loc = '/'.join(pac_name)
        
        sa_name = genome_file.split("/")
        sa_name[-1].replace('.fa', '.pac')
        sa_loc = '/'.join(sa_name)
        
        # handle error
        if not self.bwa_indexer(genome_file, amb_loc, ann_loc, bwt_loc, pac_loc, sa_loc):
            output_metadata.set_exception(
                Exception(
                    "bwa_indexer: Could not process files {}, {}.".format(*input_files)))
        output_file = None
        return ([amb_loc, ann_loc, bwt_loc, pac_loc, sa_loc], [output_metadata])

# ------------------------------------------------------------------------------
