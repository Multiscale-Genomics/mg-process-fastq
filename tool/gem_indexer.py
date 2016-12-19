import os

from pycompss.api.parameter import FILE_IN, FILE_OUT
from pycompss.api.task import task

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

from .. import common

class gemIndexerTool(Tool):
    """
    Tool for running indexers over a genome FASTA file
    """
    
    @task(file_loc=FILE_IN, idx_loc=FILE_OUT)
    def gem_indexer(self, file_loc):
        """
        GEM Indexer
        """
        cf = common()
        idx_loc = cf.gem_index_genome(file_loc)
        return True
    
    
    def run(self, input_files, metadata):
         """
        Standard function to call a task
        """
        
        # handle error
        if not self.gem_indexer(input_files[0]):
            output_metadata.set_exception(
                Exception(
                    "gem_indexer: Could not process files {}, {}.".format(*input_files)))
        output_file = None
        return ([output_file], [output_metadata])
