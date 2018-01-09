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
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task
    from utils.dummy_pycompss import compss_wait_on

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

# ------------------------------------------------------------------------------

class filterReadsTool(Tool):
    """
    Script from BS-Seeker2 for filtering FASTQ files to remove repeats
    """

    def __init__(self, configuration=None):
        """
        Init function
        """
        logger.info("BS-Seeker FilterReads wrapper")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(infile=FILE_IN, outfile=FILE_OUT, bss_path=IN)
    def bss_seeker_filter(self, infile, outfile, bss_path):
        """
        This is optional, but removes reads that can be problematic for the
        alignment of whole genome datasets.

        If performing RRBS then this step can be skipped

        This is a function that is installed as part of the BS-Seeker
        installation process.

        Parameters
        ----------
        infile : str
            Location of the FASTQ file

        Returns
        -------
        outfile : str
            Location of the filtered FASTQ file
        """
        command_line = (
            "python " + bss_path + "/FilterReads.py"
            " -i " + infile + ""
            " -o " + outfile + ".tmp"
        ).format()

        logger.info(command_line)

        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        try:
            with open(outfile, 'wb') as f_out:
                with open(outfile + '.tmp', 'rb') as f_in:
                    f_out.write(f_in.read())
        except ImportError:
            return False

        return True

    def run(self, input_files, input_metadata, output_files):
        """
        Tool for filtering duplicate entries from FASTQ files using BS-Seeker2

        Parameters
        ----------
        input_files : list
            FASTQ file
        input_metadata : list

        Returns
        -------
        array : list
            Location of the filtered FASTQ file
        """

        try:
            if "bss_path" not in input_metadata:
                raise KeyError
        except KeyError:
            logger.fatal("WGBS - BS SEEKER2: Unassigned configuration variables")
            return {}, {}

        output_metadata = {}

        print(
            "BS FILTER PARAMS:",
            input_files["fastq"],
            output_files["fastq_filtered"],
            input_metadata["bss_path"])

        results = self.bss_seeker_filter(
            input_files["fastq"],
            output_files["fastq_filtered"],
            input_metadata["bss_path"]
        )
        results = compss_wait_on(results)

        if results is False:
            logger.fatal("BS SEEKER2 Filter: run failed")
            return {}, {}

        output_metadata = {
            "fastq_filtered": Metadata(
                data_type="data_wgbs",
                file_type="fastq",
                file_path=output_files["fastq_filtered"],
                sources=[input_metadata["fastq"].file_path],
                taxon_id=input_metadata["fastq"].taxon_id,
                meta_data={
                    "tool": "bs_seeker_filter"
                }
            )
        }

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
