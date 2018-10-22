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
import tarfile

from utils import logger
from utils import remap

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on, compss_delete_file
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import compss_wait_on, compss_delete_file  # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

from mg_process_files.tool.wig_indexer import wigIndexerTool
from tool.bam_utils import bamUtilsTask
from tool.create_chrom_size import chromSizeTool

# ------------------------------------------------------------------------------


class bssMethylationCallerTool(Tool):  # pylint: disable=invalid-name
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
        returns=bool, bss_path=IN, bam_file=FILE_IN, genome_idx=FILE_IN, params=IN,
        wig_file=FILE_OUT, cgmap_file=FILE_OUT, atcgmap_file=FILE_OUT)
    def bss_methylation_caller(  # pylint: disable=no-self-use, too-many-locals, too-many-arguments
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

        g_dir = os.path.split(genome_idx)
        gi_dir = g_dir[0]

        untar_idx = True
        if "no-untar" in self.configuration and self.configuration["no-untar"] is True:
            untar_idx = False

        try:
            tar = tarfile.open(genome_idx)
            for member in tar.getmembers():
                if member.isdir():
                    gi_dir = os.path.join(g_dir[0], member.name)
                    break
            logger.info("EXTRACTING " + genome_idx + " to " + g_dir[0])
            if untar_idx is True:
                tar.extractall(path=g_dir[0])
            tar.close()
        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, g_dir[0]))

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
                if bss_pc_command_parameters[param][1] and params[param] != "":
                    command_params = command_params + [
                        bss_pc_command_parameters[param][0], params[param]]
                else:
                    if bss_pc_command_parameters[param][0] and params[param] is not False:
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
        command_params = self.get_params(self.configuration)
        self.bss_methylation_caller(
            bss_path,
            input_files["bam"],
            input_files["index"],
            command_params,
            output_files["wig_file"] + "_tmp.wig",
            output_files["cgmap_file"],
            output_files["atcgmap_file"]
        )

        chrom_size_handle = chromSizeTool(self.configuration)
        chrom_size_files, chrom_size_meta = chrom_size_handle.run(  # pylint: disable=unused-variable
            remap(input_files, "genome"),
            remap(input_metadata, "genome"),
            {
                "genome_2bit": os.path.join(
                    self.configuration["execution"], "genome.2bit"),
                "chrom_size": os.path.join(
                    self.configuration["execution"], "chrom.size")
            }
        )

        wig2bw_handle = wigIndexerTool()

        print(output_files["wig_file"] + "_tmp.wig")
        wig2bw_handle.wig2bigwig(
            output_files["wig_file"] + "_tmp.wig",
            chrom_size_files["chrom_size"],
            output_files["wig_file"]
        )
        wig2bw_handle = compss_wait_on(wig2bw_handle)

        compss_delete_file(output_files["wig_file"] + "_tmp.wig")
        compss_delete_file(chrom_size_files["chrom_size"])
        compss_delete_file(chrom_size_files["genome_2bit"])

        output_metadata = {
            "wig_file": Metadata(
                data_type="data_wgbs",
                file_type="bw",
                file_path=output_files["wig_file"],
                sources=input_metadata["bam"].sources,
                taxon_id=input_metadata["bam"].taxon_id,
                meta_data={
                    "assembly": input_metadata["bam"].meta_data["assembly"],
                    "tool": "bs_seeker_methylation_caller",
                    "parameters": command_params
                }
            ),
            "cgmap_file": Metadata(
                data_type="data_wgbs",
                file_type="tsv",
                file_path=output_files["cgmap_file"],
                sources=input_metadata["bam"].sources,
                taxon_id=input_metadata["bam"].taxon_id,
                meta_data={
                    "assembly": input_metadata["bam"].meta_data["assembly"],
                    "tool": "bs_seeker_methylation_caller",
                    "parameters": command_params
                }
            ),
            "atcgmap_file": Metadata(
                data_type="data_wgbs",
                file_type="tsv",
                file_path=output_files["atcgmap_file"],
                sources=input_metadata["bam"].sources,
                taxon_id=input_metadata["bam"].taxon_id,
                meta_data={
                    "assembly": input_metadata["bam"].meta_data["assembly"],
                    "tool": "bs_seeker_methylation_caller",
                    "parameters": command_params
                }
            )
        }

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
