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

class biobambam(Tool):
    """
    Tool to sort and filter bam files
    """

    def __init__(self, configuration=None):
        """
        Init function
        """
        print("BioBamBam2 Filter")
        Tool.__init__(self)

    @task(returns=bool, bam_file_in=FILE_IN, bam_file_out=FILE_OUT,
          isModifier=False)
    def biobambam_filter_alignments(self, bam_file_in, bam_file_out):
        """
        Sorts and filters the bam file.

        It is important that all duplicate alignments have been removed. This
        can be run as an intermediate step, but should always be run as a check
        to ensure that the files are sorted and duplicates have been removed.

        Parameters
        ----------
        bam_file_in : str
            Location of the input bam file
        bam_file_out : str
            Location of the output bam file
        tmp_dir : str
            Tmp location for intermediate files during the sorting

        Returns
        -------
        bam_file_out : str
            Location of the output bam file
        """

        td_list = bam_file_in.split("/")
        logger.info("BIOBAMBAM: bam_file_in: " + bam_file_in)
        logger.info("BIOBAMBAM: bam_file_out: " + bam_file_out)
        tmp_dir = "/".join(td_list[0:-1])

        command_line = 'bamsormadup --tmpfile=' + tmp_dir
        args = shlex.split(command_line)

        bam_tmp_out = tmp_dir + '/' + td_list[-1] + '.filtered.tmp.bam'

        logger.info("BIOBAMBAM: command_line: " + command_line)

        try:
            with open(bam_file_in, "r") as f_in:
                with open(bam_tmp_out, "w") as f_out:
                    process = subprocess.Popen(args, stdin=f_in, stdout=f_out)
                    process.wait()
        except IOError as error:
            logger.fatal("I/O error({0}): {1}".format(error.errno, error.strerror))
            return False

        try:
            with open(bam_file_out, "wb") as f_out:
                with open(bam_tmp_out, "rb") as f_in:
                    f_out.write(f_in.read())
        except IOError as error:
            logger.fatal("I/O error({0}): {1}".format(error.errno, error.strerror))
            return False

        return True

    def run(self, input_files, input_metadata, output_files):
        """
        The main function to run BioBAMBAMfilter to remove duplicates and
        spurious reads from the FASTQ files before analysis.

        Parameters
        ----------
        input_files : dict
            List of input bam file locations where 0 is the bam data file
        metadata : dict
        output_files : dict

        Returns
        -------
        output_files : dict
            Filtered bam fie.
        output_metadata : dict
            List of matching metadata dict objects
        """
        logger.info("BIOBAMBAM FILTER: Ready to run")

        results = self.biobambam_filter_alignments(input_files['input'], output_files['output'])
        results = compss_wait_on(results)

        if results is False:
            logger.fatal("BIOBAMBAM: run failed")
            return {}, {}

        logger.info("BIOBAMBAM FILTER: completed")

        output_metadata = {
            "bam": Metadata(
                data_type="data_chip_seq",
                file_type="BAM",
                file_path=output_files["output"],
                sources=[input_metadata["input"].file_path],
                taxon_id=input_metadata["input"].taxon_id,
                meta_data={
                    "assembly": input_metadata["input"].meta_data["assembly"],
                    "tool": "biobambam_filter"
                }
            )
        }

        return (
            {"bam": output_files['output']},
            output_metadata
        )

# ------------------------------------------------------------------------------
