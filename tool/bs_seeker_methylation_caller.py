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
    # from pycompss.api.api import compss_wait_on
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task  # pylint: disable=ungrouped-imports
    # from utils.dummy_pycompss import compss_wait_on  # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

from tool.bam_utils import bamUtilsTask


# ------------------------------------------------------------------------------


class bssMethylationCallerTool(Tool):
    """
    Script from BS-Seeker2 for methylation calling
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
        logger.info("BS-Seeker Methylation Caller")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(
        bss_path=IN, bam_file=FILE_IN, genome_idx=FILE_IN, params=IN,
        wig_file=FILE_OUT, cgmap_file=FILE_OUT, atcgmap_file=FILE_OUT)
    def bss_methylation_caller(  # pylint disable=no-self-use
            self, bss_path, bam_file, genome_idx, params, wig_file, cgmap_file, atcgmap_file):
        """
        Takes the merged and sorted bam file and calls the methylation sites.
        Generates a wig file of the potential sites.
        This is performed by running the external program rather than
        reimplementing the code from the main function to make it easier when
        it comes to updating the changes in BS-Seeker2

        Parameters
        ----------
        bss_path : str
            Location of the Methylation caller script
        bam_file : str
            Location of the sorted bam alignment file
        genome_idx : str
            Location of the FASTA file
        wig_file : str
            Location of the wig results file
        cgmap_file : str
            Location of the CGmap results file
        atcgmap_file : str
            Location of the ATCGmap results file

        A full description of the

        Returns
        -------
        The wig, CTmap and ATCGmap files are returned to the matching locations.
        This is managed by pyCOMPS
        """

        g_dir = genome_idx.split("/")
        g_dir = "/".join(g_dir[:-1])
        gi_dir = "/".join(g_dir[:-1])

        untar_idx = True
        if "no-untar" in self.configuration and self.configuration["no-untar"] is True:
            untar_idx = False

        try:
            tar = tarfile.open(genome_idx)
            for member in tar.getmembers():
                if member.isdir():
                    gi_dir = g_dir + "/" + member.name
                    break
            logger.info("EXTRACTING " + genome_idx + " to " + g_dir)
            if untar_idx is True:
                tar.extractall(path=g_dir)
            tar.close()
        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, g_dir))

            return False

        command_line = (
            "python " + bss_path + "/bs_seeker2-call_methylation.py "
            "--sorted --input " + str(bam_file) + " --wig " + str(wig_file) + " "
            "--CGmap " + str(cgmap_file) + " --ATCGmap " + str(atcgmap_file) + " "
            "--db " + gi_dir + " ".join(params)).format()
        logger.info("command for methyl caller: " + command_line)
        try:
            args = shlex.split(command_line)
            process = subprocess.Popen(
                args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            process.wait()

            proc_out, proc_err = process.communicate()
            logger.info("PEAK CALLER STDOUT:" + proc_out)
            logger.info("PEAK CALLER STDERR:" + proc_err)
        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, command_line))

            proc_out, proc_err = process.communicate()
            logger.fatal("IO ERROR - PEAK CALLER STDOUT:" + proc_out)
            logger.fatal("IO ERROR - PEAK CALLER STDERR:" + proc_err)
            return False

        return True

    @staticmethod
    def get_params(params):
        """
        Function to handle to extraction of commandline parameters and formatting
        them for use in the aligner for Bowtie2

        Parameters
        ----------
        params : dict
        paired : bool
            Indicate if the parameters are paired-end specific. [DEFAULT=False]

        Returns
        -------
        list
        """
        command_params = []

        bss_pc_command_parameters = {
            "bss_pc_rm_SX_param": ["--rm-SX", False],
            "bss_pc_rm_CCGG_param": ["--rm-CCGG", False],
            "bss_pc_rm_overlap_param": ["--rn-overlap", False],
            "bss_pc_read_no_param": ["--read_no", True],
        }

        for param in params:
            if param in bss_pc_command_parameters:
                if bss_pc_command_parameters[param][1]:
                    command_params = command_params + [
                        bss_pc_command_parameters[param][0], params[param]]
                else:
                    if bss_pc_command_parameters[param][0]:
                        command_params.append(bss_pc_command_parameters[param][0])

        return command_params

    def run(self, input_files, input_metadata, output_files):
        """
        Tool for methylation calling using BS-Seeker2.

        Parameters
        ----------
        input_files : list
            Sorted BAM file with the sequence alignments
        metadata : list

        Returns
        -------
        array : list
            Location of the output wig file
        """

        try:
            if "bss_path" in self.configuration:
                bss_path = self.configuration["bss_path"]
            else:
                raise KeyError
        except KeyError:
            logger.fatal("WGBS - BS SEEKER2: Unassigned configuration variables")

        bam_handler = bamUtilsTask()
        bam_handler.check_header(input_files["bam"])

        self.bss_methylation_caller(
            bss_path,
            input_files["bam"],
            input_files["index"],
            self.get_params(self.configuration),
            output_files["wig_file"],
            output_files["cgmap_file"],
            output_files["atcgmap_file"]
        )

        output_metadata = {
            "wig_file": Metadata(
                data_type="data_wgbs",
                file_type="wig",
                file_path=output_files["wig_file"],
                sources=input_metadata["bam"].sources,
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "bs_seeker_methylation_caller"
                }
            ),
            "cgmap_file": Metadata(
                data_type="data_wgbs",
                file_type="tsv",
                file_path=output_files["cgmap_file"],
                sources=input_metadata["bam"].sources,
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "bs_seeker_methylation_caller"
                }
            ),
            "atcgmap_file": Metadata(
                data_type="data_wgbs",
                file_type="tsv",
                file_path=output_files["atcgmap_file"],
                sources=input_metadata["bam"].sources,
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "bs_seeker_methylation_caller"
                }
            )
        }

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
