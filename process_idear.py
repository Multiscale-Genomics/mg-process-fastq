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

from mg_process_fastq.tool.idear import idearTool


# ------------------------------------------------------------------------------

class process_idear(Workflow):
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
        """
        logger.info("Processing DamID-Seq")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        Main run function for processing DamID-seq FastQ data. Pipeline aligns
        the FASTQ files to the genome using BWA. iDEAR is then used for peak
        calling to identify transcription factor binding sites within the
        genome.

        Currently this can only handle a single data file and a single
        background file.

        Parameters
        ----------
        input_files : dict
            Location of the initial input files required by the workflow

            bsgenome : str
                BSgenome index file

            bam_1 : str
                Location of the FASTQ aligned reads files

            bam_2 : str
                Location of the FASTQ repeat aligned reads files

            bg_bam_1 : str
                Location of the background FASTQ aligned reads files

            bg_bam_2 : str
                Location of the background FASTQ repeat aligned reads files

        metadata : dict
            Input file meta data associated with their roles

            bsgenome : str
            bam_1 : str
            bam_2 : str
            bg_bam_1 : str
            bg_bam_2 : str

        output_files : dict
            Output file locations

            bigwig : str



        Returns
        -------
        output_files : dict
            Output file locations associated with their roles, for the output

            bigwig : str
                Location of the bigwig peaks

        output_metadata : dict
            Output metadata for the associated files in output_files

            bigwig : Metadata

        """
        output_files_generated = {}
        output_metadata = {}

        # Add in BSgenome section

        logger.info("PROCESS DAMIDSEQ - DEFINED OUTPUT:", output_files)

        # iDEAR to call peaks
        idear_caller = idearTool(self.configuration)
        logger.progress("iDEAR Peak Caller", status="RUNNING")
        idear_caller.run(
            {
                "bam": input_files["bam"],
                "bg_bam": input_files["bg_bam"],
                "bsgenome": input_files["bsgenome"]
            }, {
                "bam": metadata["bam"],
                "bg_bam": metadata["bg_bam"],
                "bsgenome": metadata["bsgenome"]
            }, {
                "bigwig": output_files["bigwig"],
            }
        )
        logger.progress("iDEAR Peak Caller", status="DONE")

        print("DAMID-SEQ RESULTS:", output_metadata)
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
    result = app.launch(process_idear,
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
        description="iDEAR iDamID-seq peak calling")
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
