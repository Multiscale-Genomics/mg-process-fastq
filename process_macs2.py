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

from tool.macs2 import macs2

# ------------------------------------------------------------------------------

class process_macs2(Workflow):
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
        logger.info("Processing MACS2")
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

            bam : str
                Location of the aligned reads file

            bam_bg : str
                Location of the background aligned FASTQ reads file [OPTIONAL]

        metadata : dict
            Input file meta data associated with their roles

            bam : str

            bam_bg : str
                [OPTIONAL]

        output_files : dict
            Output file locations

            narrow_peak : str
            summits : str
            broad_peak : str
            gapped_peak : str

        Returns
        -------
        output_files : dict
            Output file locations associated with their roles, for the output

            narrow_peak : str
                Results files in bed4+1 format

            summits : str
                Results files in bed6+4 format

            broad_peak : str
                Results files in bed6+3 format

            gapped_peak : str
                Results files in bed12+3 format

        output_metadata : dict
            Output metadata for the associated files in output_files

            narrow_peak : Metadata
            summits : Metadata
            broad_peak : Metadata
            gapped_peak : Metadata
        """
        output_files_generated = {}
        output_metadata = {}

        ## MACS2 to call peaks
        macs_caller = macs2(self.configuration)
        macs_inputs = {"input": input_files["bam"]}
        macs_metadt = {"input": metadata['bam']}

        if "bg_loc" in input_files:
            macs_inputs["background"] = input_files["bam_bg"]
            macs_metadt["background"] = output_metadata['bam_bg']

        m_results_files, m_results_meta = macs_caller.run(
            macs_inputs, macs_metadt,
            # Outputs of the final step may match workflow outputs;
            # Extra entries in output_files will be disregarded.
            remap(output_files, 'narrow_peak', 'summits', 'broad_peak', 'gapped_peak'))

        if 'narrow_peak' in m_results_meta:
            output_files_generated['narrow_peak'] = m_results_files['narrow_peak']
            output_metadata['narrow_peak'] = m_results_meta['narrow_peak']

            tool_name = output_metadata['narrow_peak'].meta_data['tool']
            output_metadata['narrow_peak'].meta_data['tool_description'] = tool_name
            output_metadata['narrow_peak'].meta_data['tool'] = "process_chipseq"
        if 'summits' in m_results_meta:
            output_files_generated['summits'] = m_results_files['summits']
            output_metadata['summits'] = m_results_meta['summits']

            tool_name = output_metadata['summits'].meta_data['tool']
            output_metadata['summits'].meta_data['tool_description'] = tool_name
            output_metadata['summits'].meta_data['tool'] = "process_chipseq"
        if 'broad_peak' in m_results_meta:
            output_files_generated['broad_peak'] = m_results_files['broad_peak']
            output_metadata['broad_peak'] = m_results_meta['broad_peak']

            tool_name = output_metadata['broad_peak'].meta_data['tool']
            output_metadata['broad_peak'].meta_data['tool_description'] = tool_name
            output_metadata['broad_peak'].meta_data['tool'] = "process_chipseq"
        if 'gapped_peak' in m_results_meta:
            output_files_generated['gapped_peak'] = m_results_files['gapped_peak']
            output_metadata['gapped_peak'] = m_results_meta['gapped_peak']

            tool_name = output_metadata['gapped_peak'].meta_data['tool']
            output_metadata['gapped_peak'].meta_data['tool_description'] = tool_name
            output_metadata['gapped_peak'].meta_data['tool'] = "process_chipseq"

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
    result = app.launch(process_macs2,
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
    PARSER = argparse.ArgumentParser(description="MACS2 ChIP-seq peak calling")
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
