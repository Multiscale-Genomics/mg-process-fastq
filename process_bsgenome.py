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
from functools import wraps # pylint: disable=unused-import

import argparse

from basic_modules.workflow import Workflow
from utils import logger

from tool.forge_bsgenome import bsgenomeTool

# ------------------------------------------------------------------------------

class process_bsgenome(Workflow):
    """
    Workflow to download and pre-index a given genome
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
        logger.info("Processing Genome")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        Main run function for the indexing of genome assembly FASTA files. The
        pipeline uses Bowtie2, BWA and GEM ready for use in pipelines that
        rely on alignment.

        Parameters
        ----------
        input_files : dict
            genome : str
                Location of the FASTA input file
        metadata : dict
            genome : dict
                Required meta data
        output_files : dict
            BSgenome : str
                Location of a the BSgenome R package

        Returns
        -------
        outputfiles : dict
            List of locations for the output index files
        output_metadata : dict
            Metadata about each of the files
        """

        output_files_generated = {}
        output_metadata = {}

        # BSgenome
        logger.info("Generating BSgenome")
        bsg = bsgenomeTool()
        bsgi, bsgm = bsg.run(input_files, metadata, {'index': output_files['index']})

        try:
            output_files_generated['index'] = bsgi['index']
            output_metadata['index'] = bsgm['index']

            tool_name = output_metadata['index'].meta_data['tool']
            output_metadata['index'].meta_data['tool_description'] = tool_name
            output_metadata['index'].meta_data['tool'] = "process_bsgenome"
        except KeyError:
            logger.fatal("BSgenome indexer failed")

        return (output_files_generated, output_metadata)

# ------------------------------------------------------------------------------

def main_json(config, in_metadata, out_metadata):
    """
    Alternative main function
    -------------

    This function launches the app using configuration written in
    two json files: config.json and input_metadata.json.
    """
    # 1. Instantiate and launch the App
    logger.info("1. Instantiate and launch the App")
    from apps.jsonapp import JSONApp
    app = JSONApp()
    result = app.launch(process_bsgenome,
                        config,
                        in_metadata,
                        out_metadata)

    # 2. The App has finished
    logger.info("2. Execution finished; see " + out_metadata)

    return result

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    sys._run_from_cmdl = True  # pylint: disable=protected-access

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="Index the genome file")
    PARSER.add_argument("--config", help="Configuration file")
    PARSER.add_argument("--in_metadata", help="Location of input metadata file")
    PARSER.add_argument("--out_metadata", help="Location of output metadata file")

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    CONFIG = ARGS.config
    IN_METADATA = ARGS.in_metadata
    OUT_METADATA = ARGS.out_metadata

    RESULTS = main_json(CONFIG, IN_METADATA, OUT_METADATA)
    print(RESULTS)
