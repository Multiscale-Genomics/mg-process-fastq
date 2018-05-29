#!/usr/bin/env python

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

import argparse

from basic_modules.workflow import Workflow
from utils import logger

from tool.sleuth import sleuthTool


# ------------------------------------------------------------------------------

class process_sleuth(Workflow):
    """
    Functions for processing Chip-Seq FastQ files. Files are the aligned,
    filtered and analysed for peak calling
    """

    def __init__(self, configuration=None):
        """
        Initialise the class

        Parameters
        ----------
        configuration : dict
            a dictionary containing parameters that define how the operation
            should be carried out, which are specific to each Tool.

            kallisto_tar_config : dict
                Requires a key-pair value for each kallisto dataset. Ideally
                this should be the experiment id (eg ERR030856) and the value
                should be a description of the experiment (control, thyroid,
                brain, lung, etc.)
            sleuth_tag : str
                Tag for naming the output files with.
        """
        logger.info("Processing Sleuth")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        Main run function for processing multiple outputs from Kallisto to
        identify genes that are differentially expressed.

        Parameters
        ----------
        input_files : dict
            Location of the initial input files required by the workflow

            kallisto_tar : str
                Archive of the full set of kallisto directories

        metadata : dict
            Input file meta data associated with their roles

             kallisto_tar : str

        output_files : dict
            Output file locations

            sleuth_object : str
                Location of the Sleuth R object

        Returns
        -------
        output_files : dict
            Output file locations associated with their roles, for the output

            sleuth_object : str
                Location of the Sleuth output R object
        output_metadata : dict
            Output metadata for the associated files in output_files

            sleuth_object : Metadata

        """
        sleuth_handle = sleuthTool(self.configuration)
        s_files, s_meta = sleuth_handle.run(
            {
                "kallisto_tar": input_files["kallisto_tar"],
            }, {
                "kallisto_tar": metadata["kallisto_tar"],
            }, {
                "sleuth_object": output_files["sleuth_object"],
                "sleuth_sig_genes_table": output_files["sleuth_sig_genes_table"],
                "sleuth_image_tar": output_files["sleuth_image_tar"]
            }
        )

        logger.info("Sleuth RESULTS:", s_meta)
        return s_files, s_meta


# ------------------------------------------------------------------------------

def main_json(config, in_metadata, out_metadata):
    """
    Alternative main function
    -------------

    This function launches the app using configuration written in
    two json files: config.json and input_metadata.json.
    """
    # 1. Instantiate and launch the App
    print("1. Instantiate and launch the App")
    from apps.jsonapp import JSONApp
    app = JSONApp()
    result = app.launch(process_sleuth,
                        config,
                        in_metadata,
                        out_metadata)

    # 2. The App has finished
    print("2. Execution finished; see " + out_metadata)
    print(result)

    return result


# ------------------------------------------------------------------------------

if __name__ == "__main__":

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(
        description="Sleuth differential gene expression analysis")
    PARSER.add_argument("--config", help="Configuration file")
    PARSER.add_argument(
        "--in_metadata", help="Location of input metadata file")
    PARSER.add_argument(
        "--out_metadata", help="Location of output metadata file")
    PARSER.add_argument(
        "--local", action="store_const", const=True, default=False)

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    CONFIG = ARGS.config
    IN_METADATA = ARGS.in_metadata
    OUT_METADATA = ARGS.out_metadata
    LOCAL = ARGS.local

    if LOCAL:
        import sys
        sys._run_from_cmdl = True  # pylint: disable=protected-access

    RESULTS = main_json(CONFIG, IN_METADATA, OUT_METADATA)

    print(RESULTS)
