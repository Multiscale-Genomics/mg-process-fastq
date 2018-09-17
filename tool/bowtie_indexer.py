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

import os
import shutil
import sys
import tarfile

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT
    from pycompss.api.task import task
    # from pycompss.api.api import compss_wait_on
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task  # pylint: disable=ungrouped-imports
    # from utils.dummy_pycompss import compss_wait_on  # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

from tool.aligner_utils import alignerUtils
from tool.common import common

# ------------------------------------------------------------------------------


class bowtieIndexerTool(Tool):  # pylint: disable=invalid-name
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
        logger.info("Bowtie2 Indexer")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(file_loc=FILE_IN, index_loc=FILE_OUT)
    def bowtie2_indexer(self, file_loc, index_loc):  # pylint: disable=unused-argument, no-self-use
        """
        Bowtie2 Indexer

        Parameters
        ----------
        file_loc : str
            Location of the genome assembly FASTA file
        idx_loc : str
            Location of the output index file
        """

        au_handle = alignerUtils()
        bt2_1, bt2_2, bt2_3, bt2_4, bt2_rev1, bt2_rev2 = au_handle.bowtie_index_genome(file_loc)

        try:
            # tar.gz the index
            logger.info("BOWTIE2 - index_loc", index_loc, index_loc.replace('.tar.gz', ''))
            idx_out_pregz = index_loc.replace('.tar.gz', '.tar')

            index_dir = index_loc.replace('.tar.gz', '')
            os.mkdir(index_dir)

            shutil.move(bt2_1, index_dir)
            shutil.move(bt2_2, index_dir)
            shutil.move(bt2_3, index_dir)
            shutil.move(bt2_4, index_dir)
            shutil.move(bt2_rev1, index_dir)
            shutil.move(bt2_rev2, index_dir)

            tar = tarfile.open(idx_out_pregz, "w")
            tar.add(index_dir, arcname=os.path.split(index_dir)[1])
            tar.close()

        except (OSError, IOError) as error:
            logger.fatal("I/O error({0}): {1}".format(error.errno, error.strerror))
            return False

        common.zip_file(idx_out_pregz)
        shutil.rmtree(index_dir)

        return True

    def run(self, input_files, input_metadata, output_files):
        """
        Tool for generating assembly aligner index files for use with the
        Bowtie 2 aligner

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

        self.bowtie2_indexer(
            input_files["genome"],
            output_files["index"]
        )

        output_metadata = {
            "index": Metadata(
                data_type="sequence_mapping_index_bowtie",
                file_type="TAR",
                file_path=output_files["index"],
                sources=[input_metadata["genome"].file_path],
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "bowtie_indexer"
                }
            )
        }

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
