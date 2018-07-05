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

from mg_process_fastq.tool.bs_seeker_methylation_caller import bssMethylationCallerTool


# ------------------------------------------------------------------------------

class process_bs_seeker_peak_caller(Workflow):
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
        logger.info("Processing WGBS")
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

            bam : str
                Location for the FASTQ aligned reads bam file

            bai : str
                Location for the FASTQ aligned reads bam index file

            index : str
                Location of the index file

        metadata : dict
            Input file meta data associated with their roles

            index : dict
            bam : dict
            bai : dict

        output_files : dict
            wig_file : str
            cgmap_file : str
            atcgmap_file : str

        Returns
        -------

        wig_file : str
            Location of the wig file containing the methylation peak calls

        cgmap_file : str
            Location of the CGmap file generated by BS-Seeker2

        atcgmap_file : str
            Location of the ATCGmap file generated by BS-Seeker2

        """

        output_results_files = {}
        output_metadata = {}

        if "genome_public" in input_files:
            input_files["genome"] = input_files.pop("genome_public")
            metadata["genome"] = metadata.pop("genome_public")

        # Methylation peak caller
        peak_caller_handle = bssMethylationCallerTool(self.configuration)

        logger.progress("BSseeker2 Peak Caller", status="RUNNING")
        peak_files, peak_meta = peak_caller_handle.run(
            input_files,
            metadata,
            output_files
        )
        logger.progress("BSseeker2 Peak Caller", status="DONE")

        try:
            output_results_files["wig_file"] = peak_files["wig_file"]
            output_results_files["cgmap_file"] = peak_files["cgmap_file"]
            output_results_files["atcgmap_file"] = peak_files["atcgmap_file"]
            output_metadata["wig_file"] = peak_meta["wig_file"]
            output_metadata["cgmap_file"] = peak_meta["cgmap_file"]
            output_metadata["atcgmap_file"] = peak_meta["atcgmap_file"]

            output_metadata["wig_file"].meta_data["tool_description"] = output_metadata["wig_file"].meta_data["tool"]  # pylint: disable=line-too-long
            output_metadata["wig_file"].meta_data["tool"] = "process_bs_seeker_peak_caller"
            output_metadata["cgmap_file"].meta_data["tool_description"] = output_metadata["cgmap_file"].meta_data["tool"]  # pylint: disable=line-too-long
            output_metadata["cgmap_file"].meta_data["tool"] = "process_bs_seeker_peak_caller"
            output_metadata["atcgmap_file"].meta_data["tool_description"] = output_metadata["atcgmap_file"].meta_data["tool"]  # pylint: disable=line-too-long
            output_metadata["atcgmap_file"].meta_data["tool"] = "process_bs_seeker_peak_caller"
        except KeyError:
            logger.fatal("BS Seeker2 - Peak caller failed")
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
    result = app.launch(process_bs_seeker_peak_caller,
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
    PARSER = argparse.ArgumentParser(description="BS Seeker2 peak calling")
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
