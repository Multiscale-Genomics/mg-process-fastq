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
    from pycompss.api.parameter import FILE_IN, FILE_OUT, FILE_INOUT
    from pycompss.api.task import task
except ImportError :
    print "[Warning] Cannot import \"pycompss\" API packages."
    print "          Using mock decorators."
    
    from dummy_pycompss import *

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

from common import common

# ------------------------------------------------------------------------------

class bwaAlignerTool(Tool):
    """
    Tool for aligning sequence reads to a genome using BWA
    """

    def __init__(self):
        """
        Init function
        """
        print "BWA Aligner"
    
    @task(genome_file_loc=FILE_IN, reads_file_loc=FILE_IN, bam_loc=FILE_OUT, amb_loc=FILE_IN, ann_loc=FILE_IN, bwt_loc=FILE_IN, pac_loc=FILE_IN, sa_loc=FILE_IN)
    def bwa_aligner(self, genome_file_loc, read_file_loc, bam_loc, amb_loc, ann_loc, bwt_loc, pac_loc, sa_loc):
        """
        BWA Aligner
        
        Parameters
        ----------
        genome_file_loc : str
            Location of the genomic fasta
        read_file_loc : str
            Location of the FASTQ file
        
        Returns
        -------
        bam_loc : str
            Location of the output file
        """
        
        cf = common()
        bam_loc = cf.bwa_align_reads(genome_file_loc, read_file_loc, bam_loc)
        
        return True
    
    def run(self, input_files, metadata):
        """
        The main function to align bam files to a genome using BWA
        
        Parameters
        ----------
        input_files : list
            File 0 is the genome file location, file 1 is the FASTQ file
        
        Returns
        -------
        output : list
            First element is a list of output_bam_files, second element is the
            matching meta data
        """
        
        output_bam_file = input_files[1].replace('.fastq', '.bam')

        output_metadata = {}
        
        results = self.bwa_aligner(input_files[0], input_files[1], output_bam_file, input_files[2], input_files[3], input_files[4], input_files[5], input_files[6])
        
        # handle error
        #if not self.bwa_aligner(input_files[0], input_files[1], output_bam_file, input_files[2], input_files[3], input_files[4], input_files[5], input_files[6]):
        #    output_metadata['error'] = Exception(
        #            "bwa_indexer: Could not process files {}, {}.".format(*input_files)
        #    )
        #    output_bam_file = None
        return ([output_bam_file], [output_metadata])

# ------------------------------------------------------------------------------
