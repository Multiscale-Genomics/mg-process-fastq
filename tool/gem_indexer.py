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

from __future__ import print_function

import sys

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT
    from pycompss.api.task import task
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT
    from dummy_pycompss import task

#from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

from tool.common import common

# ------------------------------------------------------------------------------

class gemIndexerTool(Tool):
    """
    Tool for running indexers over a genome FASTA file
    """

    def __init__(self):
        """
        Init function
        """
        print ("GEM Indexer")
        Tool.__init__(self)

    @task(file_loc=FILE_IN, idx_loc=FILE_OUT)
    def gem_indexer(self, file_loc, idx_loc): # pylint: disable=unused-argument
        """
        GEM Indexer

        Parameters
        ----------
        file_loc : str
            Location of the genome assembly FASTA file
        idx_loc : str
            Location of the output index file
        """
        common_handle = common()
        idx_loc = common_handle.gem_index_genome(file_loc, file_loc)
        return True


    def run(self, input_files, output_files, metadata=None):
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

        file_out = input_files[0] + ".gem"

        # input and output share most metadata
        output_metadata = {}

        result = self.gem_indexer(input_files[0], file_out)

        return ([file_out], [output_metadata])

# ------------------------------------------------------------------------------
