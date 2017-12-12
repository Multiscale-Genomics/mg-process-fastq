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

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT
    from utils.dummy_pycompss import task
    from utils.dummy_pycompss import compss_wait_on

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool
from utils import logger

# ------------------------------------------------------------------------------

class fastqcTool(Tool):
    """
    Tool for running indexers over a genome FASTA file
    """

    def __init__(self):
        """
        Init function
        """
        logger.info("FastQC")
        Tool.__init__(self)

    @task(fastq_file=FILE_IN, report_file=FILE_OUT)
    def validate(self, fastq_file, report_file): # pylint: disable=unused-argument
        """
        GEM Indexer

        Parameters
        ----------
        FastQC_file : str
            Location of the FastQ file
        report_loc : str
            Location of the output report file
        """
        try:
            command_line = 'fastqc ' + fastq_file
            args = shlex.split(command_line)
            with open(report_file, "w") as f_out:
                process = subprocess.Popen(args, stdout=f_out)
                process.wait()

            return True
        except Exception:
            return False

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
        results = self.validate(
            output_files['fastq_file'],
            output_files['report']
        )
        results = compss_wait_on(results)

        if results is False:
            logger.fatal("FASTQC: run failed")
            return {}, {}

        output_metadata = {
            "report": Metadata(
                data_type="xml",
                file_type="HTML",
                file_path=output_files['report'],
                sources=[input_metadata["fastq_file"].file_path],
                taxon_id=input_metadata["fastq_file"].taxon_id,
                meta_data={
                    "tool": "fastqc_validator"
                }
            )
        }

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
