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
import subprocess
import sys

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


# ------------------------------------------------------------------------------

class chromSizeTool(Tool):  # pylint: disable=invalid-name
    """
    Tool for peak calling for iDamID-seq data
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
        logger.info("Forge BSgenome")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @staticmethod
    def genome_to_2bit(genome, genome_2bit):
        """
        Generate the 2bit genome file from a FASTA file

        Parameters
        ----------
        genome : str
            Location of the FASRA genome file
        genome_2bit : str
            Location of the 2bit genome file

        Returns
        -------
        bool
            True if successful, False if not.
        """
        command_line_2bit = "faToTwoBit " + genome + " " + genome_2bit

        try:
            logger.info("faToTwoBit ...")
            # args = shlex.split(command_line_2bit)
            process = subprocess.Popen(
                command_line_2bit, shell=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            process.wait()
            out, err = process.communicate()
            logger.info(out)
            if process.returncode > 0:
                logger.warn(err)
                return False
        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}) - faToTwoBit: {1}\n{2}".format(
                msg.errno, msg.strerror, command_line_2bit))
            return False

        return True

    @staticmethod
    def get_chrom_size(genome_2bit, chrom_size):
        """
        Generate the chrom.size file and identify the available chromosomes in
        the 2Bit file.

        Parameters
        ----------
        genome_2bit : str
            Location of the 2bit genome file
        chrom_size : str
            Location to save the chrom.size file to
        circ_chrom : list
            List of chromosomes that are known to be circular

        Returns
        -------
        If successful 2 lists:
            [0] : List of the linear chromosomes in the 2bit file
            [1] : List of circular chromosomes in the 2bit file

        Returns (False, False) if there is an IOError
        """
        command_line_chrom_size = "twoBitInfo " + genome_2bit + " stdout | sort -k2rn"

        try:
            logger.info("twoBitInfo ...")
            # args = shlex.split(command_line_chrom_size)
            with open(chrom_size, "w") as f_out:
                sub_proc_1 = subprocess.Popen(
                    command_line_chrom_size, shell=True,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                sub_proc_1.wait()
                out, err = sub_proc_1.communicate()
                f_out.write(out)
                logger.info(out)
                logger.warn(err)
        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0} - twoBitInfo): {1}\n{2}".format(
                msg.errno, msg.strerror, command_line_chrom_size))
            out, err = sub_proc_1.communicate()
            logger.info(out)
            logger.warn(err)
            return False

        return True

    @task(returns=int,
          genome=FILE_IN, genome_2bit=FILE_OUT, chrom_size=FILE_OUT,
          isModifier=False)
    def create_2bit_index(self, genome, genome_2bit, chrom_size):
        """
        Make BSgenome index files.Uses an R script that wraps the required code.

        Parameters
        ----------
        genome : str
        circo_chrom : str
            Comma separated list of chromosome ids that are circular in the genome
        seed_file_param : dict
            Parameters required for the function to build the seed file
        genome_2bit : str
        chrom_size : str
        seed_file : str
        bsgenome : str

        Returns
        -------

        """

        if self.genome_to_2bit(genome, genome_2bit) is False:
            return False

        if self.get_chrom_size(genome_2bit, chrom_size) is False:
            return False

        return True

    def run(self, input_files, input_metadata, output_files):
        """
        Run function to generate a 2bit genome index file
        """
        self.create_2bit_index(
            input_files["genome"],
            output_files["genome_2bit"],
            output_files["chrom_size"])

        output_metadata = {
            "genome_2bit": Metadata(
                data_type=input_metadata['genome'].data_type,
                file_type="2BIT",
                file_path=output_files["genome_2bit"],
                sources=[input_metadata["genome"].file_path],
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "create_2bit_genome"
                }
            ),
            "chrom_size": Metadata(
                data_type=input_metadata['genome'].data_type,
                file_type="TXT",
                file_path=output_files["chrom_size"],
                sources=[input_metadata["genome"].file_path],
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "create_2bit_genome"
                }
            )
        }

        return (output_files, output_metadata)
