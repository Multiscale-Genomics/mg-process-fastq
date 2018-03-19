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
import tarfile

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

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import compss_wait_on  # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata


# ------------------------------------------------------------------------------

class bssIndexerTool(Tool):
    """
    Script from BS-Seeker2 for building the index for alignment. In this case
    it uses Bowtie2.
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
        logger.info("BS-Seeker Indexer wrapper")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(
        fasta_file=FILE_IN, aligner=IN, aligner_path=IN, bss_path=IN, params=IN,
        idx_out=FILE_OUT)
    def bss_build_index(
            self, fasta_file, aligner, aligner_path, bss_path, params, idx_out
    ):  # pylint disable=no-self-use
        """
        Function to submit the FASTA file for the reference sequence and build
        the required index file used by the aligner.

        Parameters
        ----------
        fasta_file : str
            Location of the genome FASTA file
        aligner : str
            Aligner to use by BS-Seeker2. Currently only bowtie2 is available in
            this build
        aligner_path : str
            Location of the aligners binary file
        bss_path
            Location of the BS-Seeker2 libraries
        idx_out : str
            Location of the output compressed index file

        Returns
        -------
        bam_out : str
            Location of the output bam alignment file
        """

        ff_split = fasta_file.split("/")

        command_line = (
            "python " + bss_path + "/bs_seeker2-build.py " + ""
            " ".join(params) + ""
            " -f " + fasta_file + ""
            " --aligner " + aligner + " --path " + aligner_path + ""
            " --db " + "/".join(ff_split[:-1])
        ).format()

        try:
            logger.info("BS - INDEX CMD:", command_line)
            # args = shlex.split(command_line)
            process = subprocess.Popen(command_line, shell=True)
            process.wait()
        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}) - BS - INDEX CMD: {1}\n{2}".format(
                msg.errno, msg.strerror, command_line))
            return False

        try:
            # tar.gz the index
            logger.info("BS - idx_out:", idx_out, idx_out.replace('.tar.gz', ''))
            idx_out_pregz = idx_out.replace('.tar.gz', '.tar')

            logger.info("BS - idx archive:", idx_out_pregz)
            logger.info("BS - idx folder to add:", fasta_file + "_" + aligner)
            logger.info("BS - idx folder arcname:", ff_split[-1] + "_" + aligner)
            tar = tarfile.open(idx_out_pregz, "w")
            tar.add(fasta_file + "_" + aligner, arcname=ff_split[-1] + "_" + aligner)
            tar.close()

            command_line = 'pigz ' + idx_out_pregz
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()
        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, command_line))
            return False

        return True

    @staticmethod
    def get_bss_index_params(params):
        """
        Function to handle to extraction of commandline parameters and formatting
        them for use in the aligner for BWA ALN

        Parameters
        ----------
        params : dict

        Returns
        -------
        list
        """
        command_params = []

        command_parameters = {
            "bss_rrbs_param": ["-r", False],
            "bss_lower_bound_param": ["-l", True],
            "bss_upper_bound_param": ["-u", True],
            "bss_cut_format_param": ["-c", True],
        }

        for param in params:
            if param in command_parameters:
                if command_parameters[param][1]:
                    command_params = command_params + [command_parameters[param][0], params[param]]
                else:
                    if command_parameters[param][0]:
                        command_params.append(command_parameters[param][0])

        return command_params

    def run(self, input_files, input_metadata, output_files):
        """
        Tool for indexing the genome assembly using BS-Seeker2. In this case it
        is using Bowtie2

        Parameters
        ----------
        input_files : list
            FASTQ file
        metadata : list

        Returns
        -------
        array : list
            Location of the filtered FASTQ file
        """

        logger.info("WGBS - Index output files:", output_files)

        try:
            if "bss_path" in self.configuration:
                bss_path = self.configuration["bss_path"]
            else:
                raise KeyError
            if "aligner_path" in self.configuration:
                aligner_path = self.configuration["aligner_path"]
            else:
                raise KeyError
            if "aligner" in self.configuration:
                aligner = self.configuration["aligner"]
            else:
                raise KeyError
        except KeyError:
            logger.fatal("WGBS - BS SEEKER2: Unassigned configuration variables")

        command_params = self.get_bss_index_params(self.configuration)

        # handle error
        logger.info("FASTA: " + str(input_files["genome"]))
        logger.info("ALIGNER: " + str(aligner))
        logger.info("ALIGNER PATH: " + str(aligner))
        logger.info("BSS PATH: " + str(aligner))
        results = self.bss_build_index(
            input_files["genome"],
            aligner, aligner_path, bss_path, command_params,
            output_files["index"])
        results = compss_wait_on(results)

        if results is False:
            logger.fatal("BS SEEKER2 Indexer: run failed")
            return {}, {}

        output_metadata = {
            "index": Metadata(
                data_type="sequence_mapping_index_bowtie",
                file_type="TAR",
                file_path=output_files["index"],
                sources=[input_metadata["genome"].file_path],
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "bs_seeker_indexer"
                }
            )
        }

        return output_files, output_metadata

# ------------------------------------------------------------------------------
