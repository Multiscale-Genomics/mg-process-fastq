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

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task
    from utils.dummy_pycompss import compss_wait_on

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

# ------------------------------------------------------------------------------

class idear(Tool):
    """
    Tool for peak calling for MNase-seq data
    """

    def __init__(self):
        """
        Init function
        """
        print("Forge BSgenome")
        Tool.__init__(self)

    @task(
        returns=int,
        genome=FILE_IN, seed_file_param=IN,
        genome_2bit=FILE_OUT, chrom_size=FILE_OUT, seed_file=FILE_OUT, index=FILE_OUT,
        isModifier=False)
    def bsgenome_creater(
            self, genome, seed_file_param, genome_2bit, chrom_size, seed_file, index):
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
            Assembly used for teh aligned sequences
        peak_bed : str
            Location of the peak bed file

        Returns
        -------
        peak_bed : str
            Location of the collated bed file
        """

        command_line_2bit = "faToTwoBit " + genome + " " + genome_2bit
        command_line_chrom_size = "twoBitInfo " + genome_2bit + " stdout | sort -k2rn"

        args = shlex.split(command_line_2bit)
        process = subprocess.Popen(args)
        process.wait()

        args = shlex.split(command_line_chrom_size)
        with open(chrom_size, "w") as f_out:
            sub_proc = subprocess.Popen(args, shell=True, stdout=f_out)
            sub_proc.wait()

        chrom_list = []
        with open(chrom_size, "rb") as f_in:
            for line in f_in:
                line = line.split(" ")
                chrom_list.append(line[0])

        # Create the seed file
        with open(, "wb") as f_out:


        command_line = 'forge_bsgenome --file ' + genome
        logger.info("BSGENOME CMD: " + command_line)

        command_line_build = "R CMD build BSgenome."
        command_line_check = "R CMD check <pkgdir>"


        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        try:
            with open(peak_bed + ".tmp", "rb") as f_in:
                with open(peak_bed, "wb") as f_out:
                    f_out.write(f_in.read())
        except IOError:
            logger.fatal("iDEAR failed to generate peak file")
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

        sample_name = None
        background_name = None

        if "idear_sample_param" in self.configuration:
            sample_name = str(self.configuration["idear_sample_param"])
        if "idear_background_param" in self.configuration:
            background_name = str(self.configuration["idear_background_param"])

        results = self.idear_peak_calling(
            sample_name, background_name,
            input_files["bam_1"], input_files["bam_2"],
            input_files["bg_bam_1"], input_files["bg_bam_2"],
            output_files["bed"],
            input_metadata["bam"].taxon_id,
            input_metadata["bam"].meta_data["assembly"]
        )
        results = compss_wait_on(results)

        output_metadata = {
            "bed": Metadata(
                data_type=input_metadata['bam_1'].data_type,
                file_type="BED",
                file_path=output_files["bed"],
                sources=[input_metadata["bam_1"].file_path],
                taxon_id=input_metadata["bam_1"].taxon_id,
                meta_data={
                    "assembly": input_metadata["bam_1"].meta_data["assembly"],
                    "tool": "idear"
                }
            )
        }

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
