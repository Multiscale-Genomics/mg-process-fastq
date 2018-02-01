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
from utils import remap

from tool.bs_seeker_indexer import bssIndexerTool

# ------------------------------------------------------------------------------

class process_bs_seeker_index(Workflow):
    """
    Functions for aligning FastQ files with BWA
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
        logger.info("Processing BS Seeker2 Index")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        Main run function for generatigng the index files required by BS Seeker2.

        Parameters
        ----------
        input_files : dict
            List of strings for the locations of files. These should include:

            genome_fa : str
                Genome assembly in FASTA

        metadata : dict
            Input file meta data associated with their roles

            genome : str

        output_files : dict
            Output file locations

            bam : str
                Output bam file location

        Returns
        -------
        output_files : dict
            Output file locations associated with their roles, for the output

            index : str
        """
        output_results_files = {}
        output_metadata = {}

        logger.info("WGBS - BS-Seeker2 Index")
        # Build the matching WGBS genome index
        builder = bssIndexerTool(self.configuration)
        genome_idx, gidx_meta = builder.run(
            remap(input_files, "genome"),
            remap(metadata, "genome"),
            remap(output_files, "index")
        )
        output_results_files["index"] = genome_idx["index"]
        output_metadata["index"] = gidx_meta["index"]

        return output_results_files, output_metadata


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
    result = app.launch(process_bs_seeker_index,
                        config,
                        in_metadata,
                        out_metadata)

    # 2. The App has finished
    print("2. Execution finished; see " + out_metadata)

    return result

# ------------------------------------------------------------------------------

if __name__ == "__main__":

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="BS Seeker2 Indexer")
    PARSER.add_argument("--config", help="Configuration file")
    PARSER.add_argument("--in_metadata", help="Location of input metadata file")
    PARSER.add_argument("--out_metadata", help="Location of output metadata file")
    PARSER.add_argument("--local", action="store_const", const=True, default=False)

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
