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

import shlex
import subprocess
import sys

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task
    from utils.dummy_pycompss import compss_wait_on

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

# ------------------------------------------------------------------------------

class kallistoIndexerTool(Tool):
    """
    Tool for running indexers over a genome FASTA file
    """

    def __init__(self):
        """
        Init function
        """
        logger.info("Kallisto Indexer")
        Tool.__init__(self)

    @task(cdna_file_loc=FILE_IN, cdna_idx_file=FILE_OUT)
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
        logger.info("command : " + command_line)

        try:
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()
        except Exception:
            return False

        return True

    def run(self, input_files, input_metadata, output_files):
        """
        Tool for generating assembly aligner index files for use with Kallisto

        Parameters
        ----------
        input_files : list
            FASTA file location will all the cDNA sequences for a given genome
        input_metadata : list

        Returns
        --------
        array : list
            First element is a list of the index files. Second element is a
            list of the matching metadata
        """

        # file_name = input_files[0]
        # genome_idx_loc = file_name.replace('.fasta', '.idx')
        # genome_idx_loc = genome_idx_loc.replace('.fa', '.idx')

        # input and output share most metadata
        output_metadata = {}

        results = self.kallisto_indexer(
            input_files["cdna"],
            output_files["index"]
        )
        results = compss_wait_on(results)

        if results is False:
            logger.fatal("Kallisto Indexer: run failed")
            return {}, {}

        output_metadata = {
            "index": Metadata(
                data_type="index_kallisto",
                file_type="",
                file_path=output_files["index"],
                sources=[input_metadata["cdna"].file_path],
                taxon_id=input_metadata["cdna"].taxon_id,
                meta_data={
                    "assembly": input_metadata["cdna"].meta_data["assembly"],
                    "tool": "kallisto_indexer"
                }
            )
        }

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
