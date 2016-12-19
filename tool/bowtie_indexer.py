import os

from pycompss.api.parameter import FILE_IN, FILE_OUT
from pycompss.api.task import task

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

from .. import common

class bowtieIndexerTool(Tool):
    """
    Tool for running indexers over a genome FASTA file
    """
    
    @task(file_loc=FILE_IN, idx_loc=FILE_OUT)
    def bowtie2_indexer(self, file_loc):
        """
        Bowtie2 Indexer
        """
        cf = common()
        idx_loc = cf.bowtie_index_genome(file_loc)
        return True
    
    def run(self, input_files, metadata):
         """
        Standard function to call a task
        """
        
        # handle error
        if not self.bowtie_indexer(input_files[0]):
            output_metadata.set_exception(
                Exception(
                    "bowtie2_indexer: Could not process files {}, {}.".format(*input_files)))
output_file = None
        return ([output_file], [output_metadata])
