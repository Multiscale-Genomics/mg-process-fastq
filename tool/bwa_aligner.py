"""
.. Copyright 2016 EMBL-European Bioinformatics Institute
 
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

from pycompss.api.parameter import FILE_IN, FILE_OUT
from pycompss.api.task import task

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

from .. import common

# ------------------------------------------------------------------------------

class bwaAlignerTool(Tool):
    """
    Tool for aligning sequence reads to a genome using BWA
    """
    
    @task(genome_file_loc=FILE_IN, reads_file_loc=FILE_IN, bam_loc=FILE_OUT)
    def bwa_aligner(self, genome_file_loc, read_file_loc):
        """
        BWA Aligner
        """
        cf = common()
        bam_loc = cf.bwa_align_reads(genome_file_loc, reads_file_loc)
        return True
    
    def run(self, input_files, metadata):
         """
        Standard function to call a task
        """
        
        # handle error
        if not self.bwa_indexer(input_files[0]):
            output_metadata.set_exception(
                Exception(
                    "bwa_indexer: Could not process files {}, {}.".format(*input_files)))
output_file = None
        return ([output_file], [output_metadata])

# ------------------------------------------------------------------------------
