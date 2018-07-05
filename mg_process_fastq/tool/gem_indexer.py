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

    from utils.dummy_pycompss import FILE_IN, FILE_OUT  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import compss_wait_on  # pylint: disable=ungrouped-imports

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

from mg_process_fastq.tool.aligner_utils import alignerUtils
from mg_process_fastq.tool.common import common


# ------------------------------------------------------------------------------

class gemIndexerTool(Tool):  # pylint: disable=invalid-name
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
        logger.info("GEM Indexer")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(genome_file=FILE_IN, index_loc=FILE_OUT)
    def gem_indexer(self, genome_file, index_loc):  # pylint: disable=unused-argument, no-self-use
        """
        GEM Indexer

        Parameters
        ----------
        genome_file : str
            Location of the genome assembly FASTA file
        idx_loc : str
            Location of the output index file
        """
        try:
            au_handle = alignerUtils()
            au_handle.gem_index_genome(genome_file)
        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}): {1}".format(
                msg.errno, msg.strerror))
            return False

        common.zip_file(genome_file + ".gem")

        if genome_file + ".gem.gz" != index_loc:
            with open(index_loc, "wb") as f_out:
                with open(genome_file + ".gem.gz") as f_in:
                    f_out.write(f_in.read())

        return True

    def run(self, input_files, input_metadata, output_files):
        """
        Tool for generating assembly aligner index files for use with the GEM
        indexer

        Parameters
        ----------
        input_files : list
            List with a single str element with the location of the genome
            assembly FASTA file
        input_metadata : list

        Returns
        -------
        array : list
            First element is a list of the index files. Second element is a
            list of the matching metadata
        """

        # input and output share most metadata
        results = self.gem_indexer(
            input_files['genome'],
            input_files['genome'] + ".gem.gz"
        )
        results = compss_wait_on(results)

        if results is False:
            logger.fatal("GEM Indexer: run failed")
            return {}, {}

        output_files["index"] = input_files['genome'] + ".gem.gz"

        output_metadata = {
            "index": Metadata(
                data_type="sequence_mapping_index_gem",
                file_type="GEM",
                file_path=output_files['index'],
                sources=[input_files['genome']],
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "gem_indexer"
                }
            )
        }

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
