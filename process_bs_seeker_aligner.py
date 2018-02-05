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
from utils import remap

from tool.fastq_splitter import fastq_splitter
from tool.bs_seeker_aligner import bssAlignerTool
from tool.bs_seeker_filter import filterReadsTool
from tool.bs_seeker_indexer import bssIndexerTool
from tool.bs_seeker_methylation_caller import bssMethylationCallerTool

# ------------------------------------------------------------------------------

class process_bs_seeker_aligner(Workflow):
    """
    Functions for downloading and processing whole genome bisulfate sequencings
    (WGBS) files. Files are filtered, aligned and analysed for points of
    methylation
    """

    configuration = {}

    def __init__(self, configuration=None):
        """
        Initialise the class

        Parameters
        ----------
        configuration : dict
            a dictionary containing parameters that define how the operation
            should be carried out, which are specific to each Tool.
        """
        logger.info("Processing BS SEEKER2 - Aligning Reads")
        if configuration is None:
            configuration = {}
        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        This pipeline processes paired-end FASTQ files to identify
        methylated regions within the genome.

        Parameters
        ----------
        input_files : dict
            List of strings for the locations of files. These should include:

            genome_fa : str
                Genome assembly in FASTA

            fastq1 : str
                Location for the first FASTQ file for single or paired end reads

            fastq2 : str
                Location for the second FASTQ file if paired end reads

            index : str
                Location of the index file

        metadata : dict
            Input file meta data associated with their roles

            genome_fa : dict
            fastq1 : dict
            fastq2 : dict
            index : dict

        output_files : dict
            bam : str
            bai : str

        Returns
        -------

        fastq1_filtered|fastq1_filtered : str
            Locations of the filtered FASTQ files from which alignments were made

        bam|bai : str
            Location of the alignment bam file and the associated index

        """

        output_results_files = {}
        output_metadata = {}

        logger.info("BS-Seeker2 Aligner")
        # Handles the alignment of all of the split packets then merges them
        # back together.
        bss_aligner = bssAlignerTool(self.configuration)
        bam, bam_meta = bss_aligner.run(
            input_files,
            metadata,
            output_files
        )

        try:
            output_results_files["bam"] = bam["bam"]
            output_results_files["bai"] = bam["bai"]
            output_metadata["bam"] = bam_meta["bam"]
            output_metadata["bai"] = bam_meta["bai"]

            tool_name = output_metadata["bam"].meta_data["tool"]
            output_metadata["bam"].meta_data["tool_description"] = tool_name
            output_metadata["bam"].meta_data["tool"] = "process_wgbs"

            tool_name = output_metadata["bai"].meta_data["tool"]
            output_metadata["bai"].meta_data["tool_description"] = tool_name
            output_metadata["bai"].meta_data["tool"] = "process_bs_seeker_aligner"
        except KeyError:
            logger.fatal("BS SEEKER2 - Aligner failed")
            return {}, {}

        return (output_results_files, output_metadata)

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
    result = app.launch(process_bs_seeker_aligner,
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
    PARSER = argparse.ArgumentParser(description="BS-Seeker 2 Aligner")
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
