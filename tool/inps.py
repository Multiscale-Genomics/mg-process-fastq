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

class inps(Tool):
    """
    Tool for peak calling for MNase-seq data
    """

    def __init__(self, configuration=None):
        """
        Init function
        """
        logger.info("iNPS Peak Caller")
        Tool.__init__(self)

    @task(
        returns=int,
        bam_file=FILE_IN, peak_bed=FILE_OUT, inps_params=IN,
        isModifier=False)
    def inps_peak_calling(self, bam_file, peak_bed, inps_params):
        """
        Convert Bam to Bed then make Nucleosome peak calls. These are saved as
        bed files That can then get displayed on genome browsers.

        Parameters
        ----------
        bam_file : str
            Location of the aligned sequences in bam format
        peak_bed : str
            Location of the collated bed file of nucleosome peak calls

        Returns
        -------
        peak_bed : str
            Location of the collated bed file of nucleosome peak calls
        """
        bed_file = bam_file + ".bed"
        # pyenv3 = os.path.join(os.path.expanduser("~"), "bin/py3")

        command_line_1 = 'bedtools bamtobed -i ' + bam_file
        command_line_2 = "iNPS " + " ".join(inps_params)
        command_line_2 = command_line_2 + " -i " + bed_file + " -o " + peak_bed + "_tmp"

        logger.info("iNPS - cmd1:", command_line_1)
        logger.info("iNPS - cmd2:", command_line_2)

        args = shlex.split(command_line_1)
        with open(bed_file, "w") as f_out:
            sub_proc = subprocess.Popen(args, stdout=f_out)
            sub_proc.wait()

        args = shlex.split(command_line_2)
        sub_proc = subprocess.Popen(args)
        sub_proc.wait()

        with open(peak_bed + "_tmp_Gathering.like_bed", "rb") as f_in:
            with open(peak_bed, "wb") as f_out:
                f_out.write(f_in.read())

        return 0

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

        command_params = []

        if "inps_sp_param" in self.configuration:
            command_params = command_params + [
                "--s_p", str(self.configuration["inps_sp_param"])]
        if "inps_pe_max_param" in self.configuration:
            command_params = command_params + [
                "--pe_max", str(self.configuration["inps_pe_max_param"])]
        if "inps_pe_min_param" in self.configuration:
            command_params = command_params + [
                "--pe_min", str(self.configuration["inps_pe_min_param"])]

        results = self.inps_peak_calling(
            input_files["bam"],
            output_files["bed"],
            command_params
        )
        results = compss_wait_on(results)

        output_metadata = {
            "bed": Metadata(
                data_type=input_metadata['bam'].data_type,
                file_type="BED",
                file_path=output_files["bed"],
                sources=[input_metadata["bam"].file_path],
                taxon_id=input_metadata["bam"].taxon_id,
                meta_data={
                    "assembly": input_metadata["bam"].meta_data["assembly"],
                    "tool": "inps"
                }
            )
        }

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
