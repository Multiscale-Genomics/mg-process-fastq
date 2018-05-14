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
import tarfile

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task  # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

# ------------------------------------------------------------------------------


class sleuthTool(Tool):
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
        print("Sleuth Tool")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(
        returns=bool,
        sleuth_config=IN, kallisto_tar=FILE_IN,
        save_file=FILE_OUT, isModifier=False)
    def sleuth_analysis(  # pylint: disable=no-self-use
            self, sleuth_config, kallisto_tar, save_file):
        """
        Differential analysis of kallisto peak calls.

        Parameters
        ----------
        sleuth_config : dict
            Data structure describing each of the experiments and the conditions
        kallisto_tar : str
            Location of the tar file containing each of the processed Kallisto
            outputs
        save_file : str
            Location of the Sleuth R object representing the processed data

        Returns
        -------
        peak_bed : str
            Location of the collated bed file
        """

        with open("ht_config.txt", "w") as cf_handle:
            cf_handle.write("sample\tcondition")
            for row in sleuth_config:
                cf_handle.write(row + "\t" + sleuth_config[row] + "\n")

        rscript = os.path.join(os.path.dirname(__file__), "../scripts/sleuth.R")

        data_tmp_dir = kallisto_tar.split("/")
        data_tmp_dir = "/".join(data_tmp_dir[:-1])

        try:
            tar = tarfile.open(kallisto_tar)
            tar.extractall(path=data_tmp_dir)
            tar.close()
        except IOError:
            return False

        args = [
            'Rscript', rscript,
            '--config', "ht_config.txt",
            '--data_dir', kallisto_tar,
            '--save', save_file]
        logger.info("SLEUTH CMD: " + ' '.join(args))

        process = subprocess.Popen(args)
        process.wait()

        return True

    def run(self, input_files, input_metadata, output_files):
        """
        The main function to run Sleuth over a set of RNA-seq experiments
        nalysed using Kallisto to identify differentially expressed genes.

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

        self.sleuth_analysis(
            self.configuration["kallisto_tar_config"],
            input_files["kallisto_tar"],
            output_files["sleuth_object"]
        )

        output_metadata = {
            "sleuth_object": Metadata(
                data_type=input_metadata['kallisto_tar'].data_type,
                file_type="RBIN",
                file_path=output_files["sleuth_object"],
                sources=input_metadata["kallisto_tar"].sources,
                taxon_id=input_metadata["kallisto_tar"].taxon_id,
                meta_data={
                    "assembly": input_metadata["kallisto_tar"].meta_data["assembly"],
                    "tool": "sleuth"
                }
            )
        }

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
