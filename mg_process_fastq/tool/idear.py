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
    # from utils.dummy_pycompss import compss_wait_on # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

from mg_process_fastq.tool.common import common

# ------------------------------------------------------------------------------


class idearTool(Tool):  # pylint: disable=invalid-name
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
        logger.info("iDEAR Peak Caller")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @staticmethod
    def _untar_files(tar_file):
        """
        Extract bam files from the archive file created by run()

        Parameters
        ----------
        tar_file : str
            tar.gz file containing the compressed bam files

        Returns
        -------
        list
            List of the locations of the bam files. Returns an empty list if
            there is an error
        """
        included_files = []
        try:
            sample_tar_file = os.path.split(tar_file)[0]

            tar = tarfile.open(tar_file)

            for tarinfo in tar:
                if tarinfo.isreg():
                    included_files.append(os.path.join(sample_tar_file, tarinfo.name))

            tar.extractall(path=sample_tar_file)
            tar.close()
        except (OSError, IOError) as error:
            logger.fatal("UNTAR: I/O error({0}): {1}".format(error.errno, error.strerror))
            return []

        return included_files

    @task(
        returns=int,
        sample_name=IN, bg_name=IN, sample_bam_tar_file=FILE_IN,
        bg_bam__tar_file=FILE_IN, species=IN, assembly=IN, bsgenome=FILE_IN,
        peak_bw=FILE_OUT, isModifier=False)
    def idear_peak_calling(  # pylint: disable=no-self-use,too-many-locals,too-many-arguments
            self, sample_name, bg_name, sample_bam_tar_file,
            bg_bam_tar_file, common_species_name, assembly,
            bsgenome, peak_bw):
        """
        Make iDamID-seq peak calls. These are saved as bed files That can then
        get displayed on genome browsers. Uses an R script that wraps teh iDEAR
        protocol.

        Parameters
        ----------
        sample_name : str
        bg_name : str
        sample_bam_tar_file : str
            Location of the aligned sequences in bam format
        bg_bam_tar_file : str
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

        logger.info("TAR SAMPLE FILE: " + sample_bam_tar_file)
        logger.info("TAR BACKGROUND FILE: " + bg_bam_tar_file)

        sample_files = self._untar_files(sample_bam_tar_file)
        bg_files = self._untar_files(bg_bam_tar_file)

        logger.info("SAMPLE FILES: " + ", ".join(sample_files))
        logger.info("BACKGROUND FILES: " + ", ".join(bg_files))

        if not sample_files or not bg_files:
            logger.fatal("iDEAR requires sample and background bam files")
            return False

        # mkdir tmp_R_lib
        # R CMD INSTALL -l ${PWD}/tmp_R_lib BSgenome.human.grch38_1.4.2.tar.gz

        rscript = os.path.join(os.path.dirname(__file__), "../scripts/idear.R")
        rlib = os.path.join(os.getcwd(), "tmp_R_lib")

        if not os.path.exists(rlib):
            os.makedirs(rlib)

        args = shlex.split("R CMD INSTALL -l tmp_R_lib " + bsgenome)
        process = subprocess.Popen(args)
        process.wait()

        args = [
            'Rscript', rscript,
            '--sample_name', sample_name,
            '--background_name', bg_name,
            '--files', ",".join(sample_files),
            '--bg_files', ",".join(bg_files),
            '--species', str(common_species_name),
            '--assembly', assembly,
            '--output', peak_bw + '.tmp',
            '--local_lib', rlib]

        if "idear_significance" in self.configuration:
            args.append("--significance")
            args.append(self.configuration["idear_significance"])

        logger.info("iDEAR CMD: " + ' '.join(args))

        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait()
        proc_out, proc_err = process.communicate()

        try:
            with open(peak_bw + ".tmp", "rb") as f_in:
                with open(peak_bw, "wb") as f_out:
                    f_out.write(f_in.read())
        except (OSError, IOError):
            logger.fatal("iDEAR failed to generate peak file")
            logger.fatal("iDEAR stdout" + proc_out)
            logger.fatal("iDEAR stderr" + proc_err)
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

        if isinstance(input_files["bam"], list):
            tmp_sample_tar_file = os.path.join(
                os.path.split(input_files["bam"][0])[0],
                "tmp_sample_bam_files.tar"
            )
            common.tar_folder(input_files["bam"], tmp_sample_tar_file, "tmp_sample")
            input_datatype = input_metadata['bam'][0].data_type
        else:
            tmp_sample_tar_file = os.path.join(
                os.path.split(input_files["bam"])[0],
                "tmp_sample_bam_files.tar"
            )
            common.tar_folder([input_files["bam"]], tmp_sample_tar_file, "tmp_sample")
            input_datatype = input_metadata['bam'].data_type

        if isinstance(input_files["bg_bam"], list):
            tmp_background_tar_file = os.path.join(
                os.path.split(input_files["bg_bam"][0])[0],
                "tmp_background_bam_files.tar"
            )
            common.tar_folder(input_files["bg_bam"], tmp_background_tar_file, "tmp_background")
        else:
            tmp_background_tar_file = os.path.join(
                os.path.split(input_files["bg_bam"])[0],
                "tmp_background_bam_files.tar"
            )
            common.tar_folder([input_files["bg_bam"]], tmp_background_tar_file, "tmp_background")

        sample_name = None
        background_name = None
        common_name = None

        if "idear_sample_param" in self.configuration:
            sample_name = str(self.configuration["idear_sample_param"])
        if "idear_background_param" in self.configuration:
            background_name = str(self.configuration["idear_background_param"])
        if "idear_common_name" in self.configuration:
            common_name = str(self.configuration["idear_common_name"])

        self.idear_peak_calling(
            sample_name, background_name,
            tmp_sample_tar_file,
            tmp_background_tar_file,
            common_name,
            input_metadata["bsgenome"].meta_data["assembly"],
            input_files["bsgenome"],
            output_files["bigwig"]
        )

        output_metadata = {
            "bigwig": Metadata(
                data_type=input_datatype,
                file_type="BIGWIG",
                file_path=output_files["bigwig"],
                sources=[
                    input_files["bam"],
                    input_files["bg_bam"],
                    input_metadata["bsgenome"].file_path
                ],
                taxon_id=input_metadata["bsgenome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["bsgenome"].meta_data["assembly"],
                    "tool": "idear"
                }
            )
        }

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
