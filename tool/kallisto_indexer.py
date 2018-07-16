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
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task  # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

# ------------------------------------------------------------------------------


class kallistoIndexerTool(Tool):  # pylint: disable=invalid-name
    """
    Tool for running indexers over a genome FASTA file
    """

    def __init__(self, configuration=None):
        """
        Initialise the tool with its configuration.


        Parameters
        ----------
        configuration : dict
            a dictionary containing parameters that define how the operation
            should be carried out, which are specific to each Tool.
        """
        logger.info("Kallisto Indexer")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(cdna_file_loc=FILE_IN, cdna_idx_file=FILE_OUT)
    def kallisto_indexer(self, cdna_file_loc, cdna_idx_file):  # pylint: disable=no-self-use
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
        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}) - KALLISTO INDEX CMD: {1}\n{2}".format(
                msg.errno, msg.strerror, command_line))
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

        # input and output share most metadata
        output_metadata = {}

        self.kallisto_indexer(
            input_files["cdna"],
            output_files["index"]
        )

        output_metadata = {
            "index": Metadata(
                data_type="sequence_mapping_index_kallisto",
                file_type="IDX",
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
