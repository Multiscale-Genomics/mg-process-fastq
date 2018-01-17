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
from utils import remap

from tool.bwa_aligner import bwaAlignerTool

# ------------------------------------------------------------------------------

class process_bwa(Workflow):
    """
    Functions for processing Chip-Seq FastQ files. Files are the aligned,
    filtered and analysed for peak calling
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
        logger.info("Processing BWA Aligner")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        Main run function for processing ChIP-seq FastQ data. Pipeline aligns
        the FASTQ files to the genome using BWA. MACS 2 is then used for peak
        calling to identify transcription factor binding sites within the
        genome.

        Currently this can only handle a single data file and a single
        background file.

        Parameters
        ----------
        input_files : dict
            Location of the initial input files required by the workflow
            genome : str
                Genome FASTA file
            index : str
                Location of the BWA archived index files
            loc : str
                Location of the FASTQ reads files
        metadata : dict
            Input file meta data associated with their roles

            genome : str
            index : str
            loc : str
        output_files : dict
            Output file locations

            bam : str
                Output bam file location

        Returns
        -------
        output_files : dict
            Output file locations associated with their roles, for the output

            bam : str
                Aligned FASTQ short read file locations
        output_metadata : dict
            Output metadata for the associated files in output_files

            bam : Metadata
        """
        output_files_generated = {}
        output_metadata = {}

        logger.info("PROCESS ALIGNMENT - DEFINED OUTPUT:", output_files["bam"])

        bwa = bwaAlignerTool(self.configuration)
        bwa_files, bwa_meta = bwa.run(
            # ideally parameter "roles" don't change
            remap(input_files,
                  "genome", "loc", "index"),
            remap(metadata,
                  "genome", "loc", "index"),
            {"output": output_files["bam"]}
        )

        try:
            output_files_generated["bam"] = bwa_files["bam"]
            output_metadata["bam"] = bwa_meta["bam"]

            tool_name = output_metadata['bam'].meta_data['tool']
            output_metadata['bam'].meta_data['tool_description'] = tool_name
            output_metadata['bam'].meta_data['tool'] = "process_bwa"
        except KeyError:
            logger.fatal("BWA aligner failed")

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
    result = app.launch(process_bwa,
                        config,
                        in_metadata,
                        out_metadata)

    # 2. The App has finished
    print("2. Execution finished; see " + out_metadata)

    return result

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    sys._run_from_cmdl = True  # pylint: disable=protected-access

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="BWA Alignment")
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
