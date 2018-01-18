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
    from utils.dummy_pycompss import task
    from utils.dummy_pycompss import compss_wait_on

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
        Init function
        """
        print("Forge BSgenome")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(
        returns=int,
        genome=FILE_IN, circ_chrom=IN, seed_file_param=IN,
        genome_2bit=FILE_OUT, chrom_size=FILE_OUT, seed_file=FILE_OUT, bsgenome=FILE_OUT,
        isModifier=False)
    def bsgenome_creater(
            self, genome, circ_chrom, seed_file_param,
            genome_2bit, chrom_size, seed_file, bsgenome):
        """
        Make iDamID-seq peak calls. These are saved as bed files That can then
        get displayed on genome browsers. Uses an R script that wraps teh iDEAR
        protocol.

        Parameters
        ----------
        sample_name : str
        bg_name : str
        sample_bam_file_1 : str
            Location of the aligned sequences in bam format
        sample_bam_file_2 : str
            Location of the aligned sequences in bam format
        bg_bam_file_1 : str
            Location of the aligned background sequences in bam format
        bg_bam_file_2 : str
            Location of the aligned background sequences in bam format
        species : str
            Species name for the alignments
        assembly : str
            Assembly used for the aligned sequences
        peak_bed : str
            Location of the peak bed file

        Returns
        -------
        peak_bed : str
            Location of the collated bed file
        """

        command_line_2bit = "faToTwoBit " + genome + " " + genome_2bit
        command_line_chrom_size = "twoBitInfo " + genome_2bit + " stdout | sort -k2rn"

        try:
            args = shlex.split(command_line_2bit)
            process = subprocess.Popen(args)
            process.wait()

            args = shlex.split(command_line_chrom_size)
        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}) - faToTwoBit: {1}\n{2}".format(
                msg.errno, msg.strerror, command_line_2bit))
            return False

        try:
            with open(chrom_size, "w") as f_out:
                sub_proc_1 = subprocess.Popen(command_line_chrom_size, shell=True, stdout=f_out)
                sub_proc_1.wait()
        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0} - twoBitInfo): {1}\n{2}".format(
                msg.errno, msg.strerror, command_line_chrom_size))
            return False

        chrom_seq_list = []
        chrom_circ_list = []
        with open(chrom_size, "r") as f_in:
            for line in f_in:
                line = line.split("\t")
                if any(cc in line[0] for cc in circ_chrom):
                    chrom_circ_list.append(line[0])
                else:
                    chrom_seq_list.append(line[0])

        genome_split = genome_2bit.split("/")
        seed_file_param["seqnames"] = 'c("' + '","'.join(chrom_seq_list) + '")'
        if len(chrom_circ_list) > 0:
            seed_file_param["circ_seqs"] = 'c("' + '","'.join(chrom_circ_list) + '")'
        seed_file_param["circ_seqs"] = ""
        seed_file_param["seqs_srcdir"] = "/".join(genome_split[0:-1])
        seed_file_param["seqfile_name"] = genome_split[-1]

        # Create the seed file
        seed_order = [
            "Package", "Title", "Description", "Version", "organism", "common_name",
            "provider", "provider_version", "release_date", "release_name", "organism_biocview",
            "BSgenomeObjname", "seqnames", "circ_seqs", "seqs_srcdir", "seqfile_name"
        ]
        with open(seed_file, "wb") as f_out:
            for seed_key in seed_order:
                f_out.write(seed_key + ": " + seed_file_param[seed_key] + "\n")

        # Forge the BSgenomedirectory
        rscript = os.path.join(os.path.dirname(__file__), "../scripts/forge_bsgenome.R")
        command_line = "Rscript " + rscript + " --file " + seed_file
        logger.info("BSGENOME CMD: Rscript scripts/forge_bsgenome.R --file " + seed_file)

        with cd(seed_file_param["seqs_srcdir"]):
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()

            package_build = seed_file_param["seqs_srcdir"] + "/" + seed_file_param["Package"]
            command_line_build = "R CMD build " + package_build
            command_line_check = "R CMD check " + package_build

            print(command_line_build)
            args = shlex.split(command_line_build)
            process = subprocess.Popen(args)
            process.wait()

            print(command_line_check)
            args = shlex.split(command_line_check)
            process = subprocess.Popen(args)
            process.wait()

            try:
                with open(
                    package_build + "_" + seed_file_param["Version"] + ".tar.gz", "rb"
                ) as f_in:
                    with open(bsgenome, "wb") as f_out:
                        f_out.write(f_in.read())
            except IOError:
                logger.fatal("BSgenome failed to generate the index file")
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
