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

from tool.common import cd

# ------------------------------------------------------------------------------

class bsgenomeTool(Tool):
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
    def get_chrom_size(genome_2bit, chrom_size, circ_chrom):
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
            return (False, False)

        chrom_seq_list = []
        chrom_circ_list = []
        with open(chrom_size, "r") as f_in:
            for line in f_in:
                line = line.split("\t")
                if any(cc in line[0] for cc in circ_chrom):
                    chrom_circ_list.append(line[0])
                else:
                    chrom_seq_list.append(line[0])

        return (chrom_seq_list, chrom_circ_list)

    @task(
        returns=int,
        genome=FILE_IN, circ_chrom=IN, seed_file_param=IN,
        genome_2bit=FILE_OUT, chrom_size=FILE_OUT, seed_file=FILE_OUT, bsgenome=FILE_OUT,
        isModifier=False)
    def bsgenome_creater(  # pylint disable=no-self-use
            self, genome, circ_chrom, seed_file_param,
            genome_2bit, chrom_size, seed_file, bsgenome):
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

        chrom_seq_list, chrom_circ_list = self.get_chrom_size(genome_2bit, chrom_size, circ_chrom)
        if chrom_seq_list is False:
            return False

        genome_split = genome_2bit.split("/")
        seed_file_param["seqnames"] = 'c("' + '","'.join(chrom_seq_list) + '")'
        if len(chrom_circ_list) > 0:
            seed_file_param["circ_seqs"] = 'c("' + '","'.join(chrom_circ_list) + '")'
        seed_file_param["circ_seqs"] = ""
        seed_file_param["seqs_srcdir"] = "/".join(genome_split[0:-1])
        seed_file_param["seqfile_name"] = genome_split[-1]

        # Create the seed file
        seed_order = [
            "Package", "Title", "Description", "Author", "Maintainer", "License", "Version",
            "organism", "common_name", "provider", "provider_version", "release_date",
            "release_name", "organism_biocview", "BSgenomeObjname", "seqnames", "circ_seqs",
            "seqs_srcdir", "seqfile_name"
        ]
        with open(seed_file, "wb") as f_out:
            for seed_key in seed_order:
                logger.info(seed_key + ": " + seed_file_param[seed_key])
                f_out.write(seed_key + ": " + seed_file_param[seed_key] + "\n")

        # Forge the BSgenomedirectory
        rscript = os.path.join(os.path.dirname(__file__), "../scripts/forge_bsgenome.R")
        command_line = "Rscript " + rscript + " --file " + seed_file
        logger.info("BSGENOME CMD: Rscript scripts/forge_bsgenome.R --file " + seed_file)

        with cd(seed_file_param["seqs_srcdir"]):
            # args = shlex.split(command_line)
            try:
                logger.info("command_line")
                process = subprocess.Popen(
                    command_line, shell=True,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                process.wait()
                out, err = process.communicate()
                logger.info(out)
            except (IOError, OSError) as msg:
                logger.fatal("I/O error({0} - forge_bsgenome.R): {1}\n{2}".format(
                    msg.errno, msg.strerror, command_line))
                out, err = process.communicate()
                logger.info(out)
                logger.warn(err)
                return False

            package_build = seed_file_param["seqs_srcdir"] + "/" + seed_file_param["Package"]
            command_line_build = "R CMD build " + package_build
            command_line_check = "R CMD check " + package_build

            try:
                logger.info(command_line_build)
                # args = shlex.split(command_line_build)
                process = subprocess.Popen(
                    command_line_build, shell=True,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                process.wait()
                out, err = process.communicate()
                logger.info(out)
            except (IOError, OSError) as msg:
                logger.fatal("I/O error({0} - BUILD): {1}\n{2}".format(
                    msg.errno, msg.strerror, command_line_build))
                out, err = process.communicate()
                logger.info(out)
                logger.warn(err)
                return False

            try:
                logger.info(command_line_check)
                # args = shlex.split(command_line_check)
                process = subprocess.Popen(
                    command_line_check, shell=True,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                process.wait()
                out, err = process.communicate()
                logger.info(out)

                with open(package_build + ".Rcheck/00install.out", "r") as f_in:
                    logger.warn(f_in.read())
                with open(package_build + ".Rcheck/00check.log", "r") as f_in:
                    logger.warn(f_in.read())
            except (IOError, OSError) as msg:
                logger.fatal("I/O error({0} - CHECK): {1}\n{2}".format(
                    msg.errno, msg.strerror, command_line_check))
                out, err = process.communicate()
                logger.info(out)
                logger.warn(err)

                with open(package_build + ".Rcheck/00install.out", "r") as f_in:
                    logger.warn(f_in.read())
                with open(package_build + ".Rcheck/00check.log", "r") as f_in:
                    logger.warn(f_in.read())
                return False

            try:
                with open(
                    package_build + "_" + seed_file_param["Version"] + ".tar.gz", "rb"
                ) as f_in:
                    with open(bsgenome, "wb") as f_out:
                        f_out.write(f_in.read())
            except (IOError, OSError) as msg:
                logger.fatal("BSgenome failed to generate the index file")
                logger.fatal("I/O error({0}) - Change Package name: {1}\n{2}\n{3}".format(
                    msg.errno, msg.strerror,
                    package_build + "_" + seed_file_param["Version"] + ".tar.gz",
                    bsgenome))
                return False

        return True

    def run(self, input_files, input_metadata, output_files):
        """
        The main function to run iNPS for peak calling over a given BAM file
        and matching background BAM file.

        Parameters
        ----------
        input_files : list
            List of input bam file locations where 0 is the bam data file and 1
            is the matching background bam file
        metadata : dict

        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects

        """

        seed_param = {}

        seed_param["Title"] = str(self.configuration["idear_title"])
        seed_param["Description"] = str(self.configuration["idear_description"])
        seed_param["Author"] = str(self.configuration["idear_provider"])

        maintainer = str(self.configuration["idear_provider"]) + " <datasubs@ebi.ac.uk>"
        seed_param["Maintainer"] = maintainer

        seed_param["License"] = "Apache 2.0"
        seed_param["common_name"] = str(self.configuration["idear_common_name"])
        seed_param["BSgenomeObjname"] = str(self.configuration["idear_common_name"])
        seed_param["assembly"] = input_metadata["genome"].meta_data["assembly"]
        seed_param["release_name"] = input_metadata["genome"].meta_data["assembly"]
        seed_param["organism"] = str(self.configuration["idear_organism"])

        org_split = str(self.configuration["idear_organism"]).split(" ")
        seed_param["organism_biocview"] = "_".join(org_split)

        seed_param["release_date"] = str(self.configuration["idear_release_date"])
        seed_param["provider"] = str(self.configuration["idear_provider"])

        seed_param["Package"] = (
            "BSgenome." + seed_param["common_name"] + "." + seed_param["assembly"]
        )
        seed_param["Version"] = "1.4.2"
        seed_param["provider_version"] = seed_param["assembly"]

        circ_chroms = []
        if "idear_circ_chrom" in self.configuration:
            circ_chroms = str(self.configuration["idear_circ_chrom"]).split(",")

        results = self.bsgenome_creater(
            input_files["genome"],
            circ_chroms,
            seed_param,
            output_files["genome_2bit"],
            output_files["chrom_size"],
            output_files["seed_file"],
            output_files["bsgenome"]
        )
        results = compss_wait_on(results)

        output_metadata = {
            "bsgenome": Metadata(
                data_type=input_metadata['genome'].data_type,
                file_type="TAR",
                file_path=output_files["bsgenome"],
                sources=[input_metadata["genome"].file_path],
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "forge_bsgenome"
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
                    "tool": "forge_bsgenome"
                }
            ),
            "genome_2bit": Metadata(
                data_type=input_metadata['genome'].data_type,
                file_type="2BIT",
                file_path=output_files["genome_2bit"],
                sources=[input_metadata["genome"].file_path],
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "forge_bsgenome"
                }
            ),
            "seed_file": Metadata(
                data_type=input_metadata['genome'].data_type,
                file_type="TXT",
                file_path=output_files["seed_file"],
                sources=[input_metadata["genome"].file_path],
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "forge_bsgenome"
                }
            )
        }

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
