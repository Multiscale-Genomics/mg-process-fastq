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
from tool.biobambam_filter import biobambam
# from tool.idear import idear


# ------------------------------------------------------------------------------

class process_damidseq(Workflow):
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
            genome : str
                Genome FASTA file
            index : str
                Location of the BWA archived index files
            loc : str
                Location of the FASTQ reads files
            bg_loc : str
                Location of the background FASTQ reads files [OPTIONAL]
        metadata : dict
            Input file meta data associated with their roles

            genome : str
            index : str
            bg_loc : str
                [OPTIONAL]
        output_files : dict
            Output file locations

            bam [, "bam_bg"] : str
            filtered [, "filtered_bg"] : str


        Returns
        -------
        output_files : dict
            Output file locations associated with their roles, for the output

            bam [, "bam_bg"] : str
                Aligned FASTQ short read file [ and aligned background file]
                locations
            filtered [, "filtered_bg"] : str
                Filtered versions of the respective bam files

        output_metadata : dict
            Output metadata for the associated files in output_files

            bam [, "bam_bg"] : Metadata
            filtered [, "filtered_bg"] : Metadata

        """
        output_files_generated = {}
        output_metadata = {}

        logger.info("PROCESS DAMIDSEQ - DEFINED OUTPUT:", output_files["bam"])

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
            output_metadata['bam'].meta_data['tool'] = "process_chipseq"
        except KeyError:
            logger.fatal("BWA aligner failed")

        if "bg_loc" in input_files:
            # Align background files
            bwa_bg_files, bwa_bg_meta = bwa.run(
                # Small changes can be handled easily using "remap"
                remap(input_files,
                      "genome", "index", loc="bg_loc"),
                remap(metadata,
                      "genome", "index", loc="bg_loc"),
                # Intermediate outputs should be created via tempfile?
                {"output": output_files["bam_bg"]}
            )

            try:
                output_files_generated["bam_bg"] = bwa_bg_files["bam_bg"]
                output_metadata["bam_bg"] = bwa_bg_meta["bam_bg"]

                tool_name = output_metadata['bam_bg'].meta_data['tool']
                output_metadata['bam_bg'].meta_data['tool_description'] = tool_name
                output_metadata['bam_bg'].meta_data['tool'] = "process_chipseq"
            except KeyError:
                logger.fatal("Background BWA aligner failed")

        # For multiple files there will need to be merging into a single bam file

        # Filter the bams
        b3f = biobambam(self.configuration)
        b3f_files, b3f_meta = b3f.run(
            {"input": bwa_files['bam']},
            {"input": bwa_meta['bam']},
            {"output": output_files["filtered"]}
        )

        try:
            output_files_generated["filtered"] = bwa_files["filtered"]
            output_metadata["filtered"] = bwa_meta["filtered"]

            tool_name = output_metadata['filtered'].meta_data['tool']
            output_metadata['filtered'].meta_data['tool_description'] = tool_name
            output_metadata['filtered'].meta_data['tool'] = "process_chipseq"
        except KeyError:
            logger.fatal("BioBamBam filtering failed")

        if "bg_loc" in input_files:
            # Filter background aligned files
            b3f_bg_files, b3f_bg_meta = b3f.run(
                {"input": bwa_bg_files['bam']},
                {"input": bwa_bg_meta['bam']},
                {"output": output_files["filtered_bg"]}
            )

            try:
                output_files_generated["filtered_bg"] = bwa_files["filtered_bg"]
                output_metadata["filtered_bg"] = bwa_meta["filtered_bg"]

                tool_name = output_metadata['filtered_bg'].meta_data['tool']
                output_metadata['filtered_bg'].meta_data['tool_description'] = tool_name
                output_metadata['filtered_bg'].meta_data['tool'] = "process_chipseq"
            except KeyError:
                logger.fatal("Background BioBamBam filtering failed")

        ## MACS2 to call peaks
        # macs_caller = macs2(self.configuration)
        # macs_inputs = {"input": b3f_files['bam']}
        # macs_metadt = {"input": b3f_meta['bam']}

        # if "bg_loc" in input_files:
        #     macs_inputs["background"] = b3f_bg_files['bam']
        #     macs_metadt["background"] = b3f_bg_meta['bam']

        # m_results_files, m_results_meta = macs_caller.run(
        #     macs_inputs, macs_metadt,
        #     # Outputs of the final step may match workflow outputs;
        #     # Extra entries in output_files will be disregarded.
        #     remap(output_files, 'narrow_peak', 'summits', 'broad_peak', 'gapped_peak'))

        # if len(m_results_meta) == 0:
        #     logger.fatal("MACS2 peak calling failed")

        # if 'narrow_peak' in m_results_meta:
        #     output_files_generated['narrow_peak'] = m_results_files['narrow_peak']

        #     tool_name = output_metadata['narrow_peak'].meta_data['tool']
        #     output_metadata['narrow_peak'].meta_data['tool_description'] = tool_name
        #     output_metadata['narrow_peak'].meta_data['tool'] = "process_chipseq"
        # if 'summits' in m_results_meta:
        #     output_files_generated['summits'] = m_results_files['summits']

        #     tool_name = output_metadata['summits'].meta_data['tool']
        #     output_metadata['summits'].meta_data['tool_description'] = tool_name
        #     output_metadata['summits'].meta_data['tool'] = "process_chipseq"
        # if 'broad_peak' in m_results_meta:
        #     output_files_generated['broad_peak'] = m_results_files['broad_peak']

        #     tool_name = output_metadata['broad_peak'].meta_data['tool']
        #     output_metadata['broad_peak'].meta_data['tool_description'] = tool_name
        #     output_metadata['broad_peak'].meta_data['tool'] = "process_chipseq"
        # if 'gapped_peak' in m_results_meta:
        #     output_files_generated['gapped_peak'] = m_results_files['gapped_peak']

        #     tool_name = output_metadata['gapped_peak'].meta_data['tool']
        #     output_metadata['gapped_peak'].meta_data['tool_description'] = tool_name
        #     output_metadata['gapped_peak'].meta_data['tool'] = "process_chipseq"

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
    result = app.launch(process_chipseq,
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
    PARSER = argparse.ArgumentParser(description="ChIP-seq peak calling")
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
