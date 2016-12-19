import os

from pycompss.api.parameter import FILE_IN, FILE_OUT
from pycompss.api.task import task

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

from .. import common

class bwaIndexerTool(Tool):
    """
    Tool for running indexers over a genome FASTA file
    """
    
    @task(file_loc=FILE_IN, amb_loc=FILE_OUT, ann_loc=FILE_OUT, pac_loc=FILE_OUT)
    def bwa_indexer(self, file_loc):
        """
        BWA Indexer
        """
        cf = common()
        amb_loc, ann_loc, pac_loc = cf.bwa_index_genome(file_loc)
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
