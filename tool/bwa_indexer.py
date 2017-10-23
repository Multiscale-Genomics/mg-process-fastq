"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.

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
    from pycompss.api.api import compss_wait_on
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT
    from dummy_pycompss import task
    from dummy_pycompss import compss_wait_on

from basic_modules.tool import Tool

from tool.common import common

# ------------------------------------------------------------------------------

class bwaIndexerTool(Tool):
    """
    Tool for running indexers over a genome FASTA file
    """

    def __init__(self):
        """
        Init function
        """
        print("BWA Indexer")
        Tool.__init__(self)

    @task(file_loc=FILE_IN, amb_loc=FILE_OUT, ann_loc=FILE_OUT,
          bwt_loc=FILE_OUT, pac_loc=FILE_OUT, sa_loc=FILE_OUT)
    def bwa_indexer(self, file_loc, amb_loc, ann_loc, bwt_loc, pac_loc, sa_loc): # pylint: disable=unused-argument
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
        common_handler = common()
        amb_loc, ann_loc, bwt_loc, pac_loc, sa_loc = common_handler.bwa_index_genome(file_loc)
        return True

    def run(self, input_files, output_files, metadata=None):
        """
        Function to run the BWA over a genome assembly FASTA file to generate
        the matching index for use with the aligner

        Parameters
        ----------
        input_files : list
            List containing the location of the genome assembly FASTA file
        meta_data : list
        output_files : list
            List of outpout files generated

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
        output_metadata = {}

        output_files = [
            input_files[0] + ".amb",
            input_files[0] + ".ann",
            input_files[0] + ".bwt",
            input_files[0] + ".pac",
            input_files[0] + ".sa"
        ]

        results = self.bwa_indexer(
            input_files[0],
            input_files[0] + ".amb",
            input_files[0] + ".ann",
            input_files[0] + ".bwt",
            input_files[0] + ".pac",
            input_files[0] + ".sa"
        )

        results = compss_wait_on(results)

        return (output_files, [output_metadata])

# ------------------------------------------------------------------------------
