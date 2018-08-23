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
import shlex
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

from tool.common import common

# ------------------------------------------------------------------------------


class sleuthTool(Tool):  # pylint: disable=invalid-name
    """
    Tool for analysing gene differential expression using Sleuth
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

    @staticmethod
    def extract_kallisto_tar(data_tmp_dir, kallisto_tar):
        """
        Function to extract the kallisto tar file
        """
        try:
            tar = tarfile.open(kallisto_tar)
            tar.extractall(path=data_tmp_dir)

            for member in tar.getmembers():
                if member.isfile():
                    member_dir = member.name.split("/")
                    member_dir = data_tmp_dir + "/" + "/".join(member_dir[:-1])
                    logger.info("data_tmp_dir: " + data_tmp_dir)
                    logger.info("Member: " + member.name)
                    tar_sub = tarfile.open(data_tmp_dir + "/" + member.name)
                    tar_sub.extractall(path=member_dir)
                    tar_sub.close()

            tar.close()
        except IOError as error:
            logger.fatal(
                "IO ERROR {0}: Failed to extract all files:\n{1}".format(
                    error.errno, error.strerror
                )
            )
            return False

        return True

    @staticmethod
    def write_config_file(loc, sleuth_config):
        """
        Function to create the Sleuth configuration file descibing the
        experiments
        """
        with open(loc, "w") as cf_handle:
            dataset = next(iter(sleuth_config))
            column_keys = []
            for column in list(sleuth_config[dataset].keys()):
                if column != "sample":
                    column_keys.append(column)
            cf_handle.write("sample\t" + "\t".join(column_keys) + "\n")
            for row in sleuth_config:
                row_string = row
                for column in column_keys:
                    row_string += "\t" + sleuth_config[row][column]

                cf_handle.write(row_string + "\n")

        return True

    @task(
        returns=bool,
        sleuth_config=IN, kallisto_tar=FILE_IN,
        save_file=FILE_OUT, save_table=FILE_OUT, level=IN, isModifier=False)
    def sleuth_analysis(
            self, sleuth_config, kallisto_tar, save_file, save_table, level):
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
        save_table : str
            Location of the top most differentially expressed genes
        level : float
            Significance value for the cut off for filtering the top genes when
            saving the table

        Returns
        -------
        bool
            True if the task has completed successfully
        """

        data_tmp_dir = kallisto_tar.split("/")
        data_tmp_dir = "/".join(data_tmp_dir[:-1])

        self.extract_kallisto_tar(data_tmp_dir, kallisto_tar)
        self.write_config_file(data_tmp_dir + "/ht_config.txt", sleuth_config)

        rscript = os.path.join(os.path.dirname(__file__), "../scripts/sleuth.R")

        args = [
            'Rscript', rscript,
            '--config', data_tmp_dir + "/ht_config.txt",
            '--data_dir', kallisto_tar.replace(".tar.gz", ""),
            '--save', save_file,
            '--deg', save_table,
            '--degl', level]
        logger.info("SLEUTH CMD: " + ' '.join(args))

        process = subprocess.Popen(args)
        process.wait()

        return True

    @task(returns=bool,
          covariants=IN, tag=IN, sig_level=IN,
          sleuth_object=FILE_IN, image_tar=FILE_OUT)
    def sleuth_visualise(self, sleuth_object, image_tar, sig_level=None, covariants=None, tag=None):
        """
        Running each of the visualisation scripts to generate plots related to
        set of expression datasets

        Parameters
        ----------
        sleuth_object : str
            Location of the save_file sleuth object generated by the
            self.sleuth_analysis() function
        image_tar : str
            Location of the archive object that all of the generated images will
            be save to.

        Returns
        -------
        bool
            True if the task has completed successfully
        """

        if tag is None:
            tag = ""
        if covariants is None:
            covariants = []

        images_to_save = [
            sleuth_object + "_sample_heatmap_" + tag + ".png",
            sleuth_object + "_transcript_heatmap_" + tag + ".png",
        ]

        r_cmds = [
            [
                "Rscript", os.path.join(os.path.dirname(__file__), "../scripts/sleuth_heatmap.R"),
                "--file", sleuth_object,
                "--tag", tag,
                "--degl", sig_level
            ],
        ]

        for covar in covariants:
            r_cmds.append(
                [
                    "Rscript", os.path.join(os.path.dirname(__file__), "../scripts/sleuth_pca.R"),
                    "--file", sleuth_object,
                    "--tag", covar,
                    "--covariant", covar
                ]
            )
            images_to_save.append(sleuth_object + "_pca_" + covar + ".png")
            images_to_save.append(sleuth_object + "_pca1_" + covar + ".png")

        for cmd in r_cmds:
            logger.info(" ".join(cmd))
            process = subprocess.Popen(cmd)
            process.wait()

        tar = tarfile.open(image_tar.replace(".gz", ""), "w")
        for img in images_to_save:
            tar.add(img, arcname="results")
        tar.close()

        common.zip_file(image_tar.replace(".gz", ""), 2)

    def run(self, input_files, input_metadata, output_files):
        """
        The main function to run Sleuth over a set of RNA-seq experiments
        analysed using Kallisto to identify differentially expressed genes.

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

        dataset = next(iter(self.configuration["kallisto_tar_config"]))
        covariants = []
        for column in list(self.configuration["kallisto_tar_config"][dataset].keys()):
            if column != "sample":
                covariants.append(column)

        self.sleuth_analysis(
            self.configuration["kallisto_tar_config"],
            input_files["kallisto_tar"],
            output_files["sleuth_object"],
            output_files["sleuth_sig_genes_table"],
            str(self.configuration["sleuth_sig_level"])
        )

        logger.info("COVARIANTS: " + ", ".join(covariants))

        self.sleuth_visualise(
            output_files["sleuth_object"],
            output_files["sleuth_image_tar"],
            str(self.configuration["sleuth_sig_level"]),
            covariants,
            str(self.configuration["sleuth_tag"])
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
            ),
            "sleuth_sig_genes_table": Metadata(
                data_type=input_metadata['kallisto_tar'].data_type,
                file_type="TSV",
                file_path=output_files["sleuth_sig_genes_table"],
                sources=input_metadata["kallisto_tar"].sources,
                taxon_id=input_metadata["kallisto_tar"].taxon_id,
                meta_data={
                    "assembly": input_metadata["kallisto_tar"].meta_data["assembly"],
                    "tool": "sleuth"
                }
            ),
            "sleuth_image_tar": Metadata(
                data_type=input_metadata['kallisto_tar'].data_type,
                file_type="TSV",
                file_path=output_files["sleuth_image_tar"],
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
