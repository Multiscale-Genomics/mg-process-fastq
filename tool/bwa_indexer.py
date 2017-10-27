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
import shlex
import shutil
import subprocess
import sys
import tarfile

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
from basic_modules.metadata import Metadata
from utils import logger

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
        logger.info("BWA Indexer")
        Tool.__init__(self)

    @task(file_loc=FILE_IN, idx_out=FILE_OUT)
    def bwa_indexer(self, file_loc, idx_out): # pylint: disable=unused-argument
        """
        BWA Indexer

        Parameters
        ----------
        file_loc : str
            Location of the genome assebly FASTA file
        idx_out : str
            Location of the output index file

        Returns
        -------
        bool
        """
        try:
            common_handler = common()
            amb_loc, ann_loc, bwt_loc, pac_loc, sa_loc = common_handler.bwa_index_genome(file_loc)

            # tar.gz the index
            print("BS - idx_out", idx_out, idx_out.replace('.tar.gz', ''))
            idx_out_pregz = idx_out.replace('.tar.gz', '.tar')

            index_dir = idx_out.replace('.tar.gz', '')
            os.mkdir(index_dir)

            idx_split = index_dir.split("/")

            shutil.move(amb_loc, index_dir)
            shutil.move(ann_loc, index_dir)
            shutil.move(bwt_loc, index_dir)
            shutil.move(pac_loc, index_dir)
            shutil.move(sa_loc, index_dir)

            index_folder = idx_split[-1]

            tar = tarfile.open(idx_out_pregz, "w")
            tar.add(index_dir, arcname=index_folder)
            tar.close()

            command_line = 'pigz ' + idx_out_pregz
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()

            return True
        except Exception:
            return False

    def run(self, input_files, input_metadata, output_files):
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
        output_files : dict
            index : str
                Location of the index file defined in the input parameters
        output_metadata : dict
            index : Metadata
                Metadata relating to the index file
        """
        results = self.bwa_indexer(
            input_files["genome"],
            output_files["index"]
        )
        results = compss_wait_on(results)

        if results is False:
            logger.fatal("BWA Indexer: run failed")
            return {}, {}

        output_metadata = {
            "index": Metadata(
                data_type="sequence_mapping_index_bwa",
                file_type="TAR",
                file_path=output_files["index"],
                sources=[input_metadata["genome"].file_path],
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "bwa_indexer"
                }
            )
        }

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
