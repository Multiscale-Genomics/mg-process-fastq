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

from tool.inps import inps

# ------------------------------------------------------------------------------

class process_iNPS(Workflow):
    """
    Functions for improved nucleosome positioning algorithm
    (iNPS).  Bam Files are analysed for peaks for nucleosome positioning
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
        logger.info("Processing iNPS")
        if configuration is None:
            configuration = {}
        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        This pipeline processes bam files to identify
        nucleosome regions within the genome and generates bed files.

        Parameters
        ----------
        input_files : dict
            bam_file : str
            Location of the aligned sequences in bam format

        output_files : dict
            peak_bed : str
            Location of the collated bed file of nucleosome peak calls

        Returns
        -------

        peak_bed : str
            Location of the collated bed file of nucleosome peak calls

        """

        output_results_files = {}
        output_metadata = {}

        logger.info("iNPS")

        # Process the MNAse-seq reads to find nucleosome
        logger.info("iNPS")
        inps_r = inps(self.configuration)
        bamf, bed_meta = inps_r.run(
            {"bam": input_files["bam"]},
            {"bam": metadata["bam"]},
            {"bed": output_files["bed"]}
        )

        try:
            output_results_files["bed"] = bamf["bed"]
            output_metadata["bed"] = bed_meta["bed"]
            tool_name = output_metadata["bed"].meta_data["tool"]
            output_metadata["bed"].meta_data["tool_description"] = tool_name
            output_metadata["bed"].meta_data["tool"] = "process_iNPS"
        except KeyError:
            logger.fatal("iNPS : Error while processing")
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
    result = app.launch(process_iNPS,
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
    PARSER = argparse.ArgumentParser(description="iNPS peak calling")
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
