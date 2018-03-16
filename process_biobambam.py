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

# Required for ReadTheDocs
from functools import wraps  # pylint: disable=unused-import

import argparse

from basic_modules.workflow import Workflow
from utils import logger

from tool.biobambam_filter import biobambam


# ------------------------------------------------------------------------------

class process_biobambam(Workflow):  # pylint disable=too-few-public-methods, invalid-name
    """
    Functions for filtering FastQ alignments with BioBamBam.
    """

    def __init__(self, configuration=None):
        """
        Initialise the class

        Parameters
        ----------
        configuration : dict
            a dictionary containing parameters that define how the operation
            should be carried out, which are specific to each Tool.
        """
        logger.info("Processing BioBamBam Filtering")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        Main run function for filtering FastQ aligned reads using BioBamBam.

        Parameters
        ----------
        input_files : dict
            Location of the initial input files required by the workflow

            bam : str
                Location of BAM file

        metadata : dict
            Input file meta data associated with their roles

            bam : str

        output_files : dict
            Output file locations

            filtered : str

        Returns
        -------
        output_files : dict
            Output file locations associated with their roles, for the output

            filtered : str
                Filtered version of the bam file

        output_metadata : dict
            Output metadata for the associated files in output_files

            filtered : Metadata
        """
        output_files_generated = {}
        output_metadata = {}

        # Filter the bam
        b3f = biobambam(self.configuration)
        b3f_files, b3f_meta = b3f.run(
            {"input": input_files['bam']},
            {"input": metadata['bam']},
            {"output": output_files["filtered"]}
        )

        try:
            output_files_generated["filtered"] = b3f_files["bam"]
            output_metadata["filtered"] = b3f_meta["bam"]

            tool_name = output_metadata['filtered'].meta_data['tool']
            output_metadata['filtered'].meta_data['tool_description'] = tool_name
            output_metadata['filtered'].meta_data['tool'] = "process_biobambam"
        except KeyError:
            logger.fatal("BioBamBam filtering failed")

        return output_files_generated, output_metadata


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
    result = app.launch(process_biobambam,
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
    PARSER = argparse.ArgumentParser(description="BioBamBam filtering")
    PARSER.add_argument(
        "--config", help="Configuration file")
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
