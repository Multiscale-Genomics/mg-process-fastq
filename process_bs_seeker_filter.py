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

from tool.bs_seeker_filter import filterReadsTool

# ------------------------------------------------------------------------------

class process_bsFilter(Workflow):
    """
    Functions for filtering FASTQ files.
    Files are filtered for removal of duplicate reads.
    Low quality reads in qseq file can also be filtered.
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
        logger.info("Processing BS Filter")
        if configuration is None:
            configuration = {}
        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        This pipeline processes FASTQ files to filter duplicate entries

        Parameters
        ----------
        input_files : dict
            List of strings for the locations of files. These should include:

            fastq1 : str
                Location for the first FASTQ file for single or paired end reads

            fastq2 : str
                Location for the second FASTQ file if paired end reads [OPTIONAL]

        metadata : dict
            Input file meta data associated with their roles

            fastq1 : str

            fastq2 : str
                [OPTIONAL]

        output_files : dict

            fastq1_filtered : str

            fastq2_filtered : str
                [OPTIONAL]

        Returns
        -------

        fastq1_filtered|fastq1_filtered : str
            Locations of the filtered FASTQ files from which alignments were made

        fastq2_filtered|fastq2_filtered : str
            Locations of the filtered FASTQ files from which alignments were made

        """

        output_results_files = {}
        output_metadata = {}

        logger.info("BS-Filter")

        frt = filterReadsTool(self.configuration)
        fastq1f, filter1_meta = frt.run(
            {"fastq": input_files["fastq1"]},
            {"fastq": metadata["fastq1"]},
            {"fastq_filtered": output_files["fastq1_filtered"]}
        )

        try:
            output_results_files["fastq1_filtered"] = fastq1f["fastq_filtered"]
            output_metadata["fastq1_filtered"] = filter1_meta["fastq_filtered"]
            tool_name = output_metadata["fastq1_filtered"].meta_data["tool"]
            output_metadata["fastq1_filtered"].meta_data["tool_description"] = tool_name
            output_metadata["fastq1_filtered"].meta_data["tool"] = "process_wgbs"
        except KeyError:
            logger.fatal("WGBS - FILTER: Error while filtering")
            return {}, {}

        if "fastq2" in input_files:
            logger.info("WGBS - Filter background")
            fastq2f, filter2_meta = frt.run(
                {"fastq": input_files["fastq2"]},
                {"fastq": metadata["fastq2"]},
                {"fastq_filtered": output_files["fastq2_filtered"]}
            )

            try:
                output_results_files["fastq2_filtered"] = fastq2f["fastq_filtered"]
                output_metadata["fastq2_filtered"] = filter2_meta["fastq_filtered"]

                tool_name = output_metadata["fastq2_filtered"].meta_data["tool"]
                output_metadata["fastq2_filtered"].meta_data["tool_description"] = tool_name
                output_metadata["fastq2_filtered"].meta_data["tool"] = "process_wgbs"
            except KeyError:
                logger.fatal("WGBS - FILTER (background): Error while filtering")
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
    result = app.launch(process_bsFilter,
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
    PARSER = argparse.ArgumentParser(description="WGBS Filter")
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
