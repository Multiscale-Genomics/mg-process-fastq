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

from mg_process_fastq.tool.trimgalore import trimgalore


# ------------------------------------------------------------------------------

class process_trim_galore(Workflow):  # pylint: disable=invalid-name,too-few-public-methods
    """
    Functions for filtering FASTQ files.
    Files are filtered for removal of duplicate reads.
    Low quality reads in qseq file can also be filtered.
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
        logger.info("Processing trim galore")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        This pipeline processes FASTQ files to trim low quality base calls and
        adapter sequences

        Parameters
        ----------
        input_files : dict
            List of strings for the locations of files. These should include:

            fastq : str
                Location for the first FASTQ file for single or paired end reads


        metadata : dict
            Input file meta data associated with their roles

        output_files : dict

            fastq_trimmed : str

        Returns
        -------

        fastq_trimmed|fastq_trimmed : str
            Locations of the filtered FASTQ files from which trimmings were made

        """

        output_results_files = {}
        output_metadata = {}

        logger.info("trim_galore")

        trimg = trimgalore(self.configuration)

        logger.progress("TrimGalore", status="RUNNING")
        if "fastq2" in input_files:
            trim_files, trim_meta = trimg.run(
                {
                    "fastq1": input_files["fastq1"],
                    "fastq2": input_files["fastq2"]
                },
                {
                    "fastq1": metadata["fastq1"],
                    "fastq2": metadata["fastq2"]
                },
                {
                    "fastq1_trimmed": output_files["fastq1_trimmed"],
                    "fastq2_trimmed": output_files["fastq2_trimmed"],
                    "fastq1_report": output_files["fastq1_report"],
                    "fastq2_report": output_files["fastq2_report"]
                }
            )

            try:
                output_results_files["fastq2_trimmed"] = trim_files["fastq2_trimmed"]
                output_metadata["fastq2_trimmed"] = trim_meta["fastq2_trimmed"]
                tool_name = output_metadata["fastq2_trimmed"].meta_data["tool"]
                output_metadata["fastq2_trimmed"].meta_data["tool_description"] = tool_name
                output_metadata["fastq2_trimmed"].meta_data["tool"] = "process_trim_galore"
            except KeyError:
                logger.fatal("Trim Galore fastq2: Error while trimming")
                return {}, {}
        else:
            trim_files, trim_meta = trimg.run(
                {"fastq1": input_files["fastq1"]},
                {"fastq1": metadata["fastq1"]},
                {
                    "fastq1_trimmed": output_files["fastq1_trimmed"],
                    "fastq1_report": output_files["fastq1_report"]
                }
            )
        logger.progress("TrimGalore", status="DONE")

        try:
            output_results_files["fastq1_trimmed"] = trim_files["fastq1_trimmed"]
            output_metadata["fastq1_trimmed"] = trim_meta["fastq1_trimmed"]
            tool_name = output_metadata["fastq1_trimmed"].meta_data["tool"]
            output_metadata["fastq1_trimmed"].meta_data["tool_description"] = tool_name
            output_metadata["fastq1_trimmed"].meta_data["tool"] = "process_trim_galore"
        except KeyError:
            logger.fatal("Trim Galore fastq1: Error while trimming")
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
    result = app.launch(process_trim_galore,
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
    PARSER = argparse.ArgumentParser(description="Trim galore")
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
