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
from basic_modules.metadata import Metadata
from utils import logger
from utils import remap

from dmp import dmp

from tool.bwa_aligner import bwaAlignerTool
from tool.inps import inps

# ------------------------------------------------------------------------------

class process_mnaseseq(Workflow):
    """
    Functions for downloading and processing Mnase-seq FastQ files. Files are
    downloaded from the European Nucleotide Archive (ENA), then aligned,
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
        logger.info("Processing MNase-Seq")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        Main run function for processing MNase-Seq FastQ data. Pipeline aligns
        the FASTQ files to the genome using BWA. iNPS is then used for peak
        calling to identify nucleosome position sites within the genome.

        Parameters
        ----------
        files_ids : list
            List of file locations
        metadata : list
            Required meta data

        Returns
        -------
        outputfiles : list
            List of locations for the output bam, bed and tsv files
        """

        output_metadata = {}

        bwa = bwaAlignerTool()
        bwa_files, bwa_meta = bwa.run(
            remap(input_files, "genome", "loc", "index"),
            remap(metadata, "genome", "loc", "index"),
            {"output": output_files["bam"]}
        )

        output_metadata["bam"] = bwa_meta["bam"]
        tool_name = output_metadata['bed'].meta_data['tool']
        output_metadata['bed'].meta_data['tool_description'] = tool_name
        output_metadata['bed'].meta_data['tool'] = "process_mnaseseq"

        inps_tool = inps()
        out_peak_bed, out_peak_bed_meta = inps_tool.run(
            remap(bwa_files, "bam"),
            remap(bwa_meta, "bam"),
            {"bed": output_files["bed"]}
        )

        output_metadata["bed"] = out_peak_bed_meta["bed"]
        tool_name = output_metadata['bed'].meta_data['tool']
        output_metadata['bed'].meta_data['tool_description'] = tool_name
        output_metadata['bed'].meta_data['tool'] = "process_mnaseseq"

        return (output_files, output_metadata)

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
    result = app.launch(process_mnaseseq,
                        config,
                        in_metadata,
                        out_metadata)

    # 2. The App has finished
    print("2. Execution finished; see " + out_metadata)
    print(result)

    return result


# ------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    sys._run_from_cmdl = True  # pylint: disable=protected-access

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="MNase-seq peak calling")
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
