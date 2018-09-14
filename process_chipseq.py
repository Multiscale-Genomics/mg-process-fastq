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

from mg_process_fastq.tool.bwa_aligner import bwaAlignerTool
from mg_process_fastq.tool.biobambam_filter import biobambam
from mg_process_fastq.tool.macs2 import macs2


# ------------------------------------------------------------------------------

class process_chipseq(Workflow):  # pylint: disable=invalid-name,too-few-public-methods
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
        logger.info("Processing ChIP-Seq")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):  # pylint: disable=too-many-branches,too-many-locals,too-many-statements,line-too-long
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

            fastq2 : str
                Location of the paired end FASTQ file [OPTIONAL]

            bg_loc : str
                Location of the background FASTQ reads files [OPTIONAL]

            fastq2_bg : str
                Location of the paired end background FASTQ reads files [OPTIONAL]

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
            narrow_peak : str
            summits : str
            broad_peak : str
            gapped_peak : str

        Returns
        -------
        output_files : dict
            Output file locations associated with their roles, for the output

            bam [, "bam_bg"] : str
                Aligned FASTQ short read file [ and aligned background file]
                locations
            filtered [, "filtered_bg"] : str
                Filtered versions of the respective bam files
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

            bam [, "bam_bg"] : Metadata
            filtered [, "filtered_bg"] : Metadata
            narrow_peak : Metadata
            summits : Metadata
            broad_peak : Metadata
            gapped_peak : Metadata
        """
        output_files_generated = {}
        output_metadata = {}

        logger.info("PROCESS CHIPSEQ - DEFINED OUTPUT:", output_files["bam"])

        if "genome_public" in input_files:
            align_input_files = remap(
                input_files, genome="genome_public", loc="loc", index="index_public")
            align_input_file_meta = remap(
                metadata, genome="genome_public", loc="loc", index="index_public")
        else:
            align_input_files = remap(input_files, "genome", "loc", "index")
            align_input_file_meta = remap(metadata, "genome", "loc", "index")

        if "fastq2" in input_files:
            align_input_files["fastq2"] = input_files["fastq2"]
            align_input_file_meta["fastq2"] = metadata["fastq2"]

        logger.progress("BWA Aligner", status="RUNNING")
        bwa = bwaAlignerTool(self.configuration)
        bwa_files, bwa_meta = bwa.run(
            align_input_files,
            align_input_file_meta,
            {"output": output_files["bam"]}
        )
        logger.progress("BWA Aligner", status="DONE")
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
            if "genome_public" in input_files:
                align_input_files_bg = remap(
                    input_files, genome="genome_public", index="index_public", loc="bg_loc")
                align_input_file_meta_bg = remap(
                    metadata, genome="genome_public", index="index_public", loc="bg_loc")
            else:
                align_input_files_bg = remap(input_files, "genome", "index", loc="bg_loc")
                align_input_file_meta_bg = remap(metadata, "genome", "index", loc="bg_loc")

            if "fastq2" in input_files:
                align_input_files_bg["fastq2"] = input_files["fastq2_bg"]
                align_input_file_meta_bg["fastq2"] = metadata["fastq2_bg"]

            logger.progress("BWA Aligner - Background", status="RUNNING")
            bwa_bg_files, bwa_bg_meta = bwa.run(
                align_input_files_bg,
                align_input_file_meta_bg,
                {"output": output_files["bam_bg"]}
            )
            logger.progress("BWA Aligner - Background", status="DONE")

            try:
                output_files_generated["bam_bg"] = bwa_bg_files["bam_bg"]
                output_metadata["bam_bg"] = bwa_bg_meta["bam_bg"]

                tool_name = output_metadata['bam_bg'].meta_data['tool']
                output_metadata['bam_bg'].meta_data['tool_description'] = tool_name
                output_metadata['bam_bg'].meta_data['tool'] = "process_chipseq"
            except KeyError:
                logger.fatal("Background BWA aligner failed")

        # Filter the bams
        b3f = biobambam(self.configuration)
        logger.progress("BioBamBam", status="RUNNING")
        b3f_files, b3f_meta = b3f.run(
            {"input": bwa_files['bam']},
            {"input": bwa_meta['bam']},
            {"output": output_files["filtered"]}
        )
        logger.progress("BioBamBam", status="DONE")

        try:
            output_files_generated["filtered"] = b3f_files["bam"]
            output_metadata["filtered"] = b3f_meta["bam"]

            tool_name = output_metadata['filtered'].meta_data['tool']
            output_metadata['filtered'].meta_data['tool_description'] = tool_name
            output_metadata['filtered'].meta_data['tool'] = "process_chipseq"
        except KeyError:
            logger.fatal("BioBamBam filtering failed")

        if "bg_loc" in input_files:
            # Filter background aligned files
            logger.progress("BioBamBam Background", status="RUNNING")
            b3f_bg_files, b3f_bg_meta = b3f.run(
                {"input": bwa_bg_files['bam']},
                {"input": bwa_bg_meta['bam']},
                {"output": output_files["filtered_bg"]}
            )
            logger.progress("BioBamBam Background", status="DONE")

            try:
                output_files_generated["filtered_bg"] = b3f_bg_files["bam"]
                output_metadata["filtered_bg"] = b3f_bg_meta["bam"]

                tool_name = output_metadata['filtered_bg'].meta_data['tool']
                output_metadata['filtered_bg'].meta_data['tool_description'] = tool_name
                output_metadata['filtered_bg'].meta_data['tool'] = "process_chipseq"
            except KeyError:
                logger.fatal("Background BioBamBam filtering failed")

        # MACS2 to call peaks
        # Duplicates have already been filtered so MACS2 does not need to due
        # any further filtering
        self.configuration["macs_keep-dup_param"] = "all"

        macs_caller = macs2(self.configuration)
        macs_inputs = {"bam": output_files_generated["filtered"]}
        macs_metadt = {"bam": output_metadata['filtered']}

        if "bg_loc" in input_files:
            macs_inputs["bam_bg"] = output_files_generated["filtered_bg"]
            macs_metadt["bam_bg"] = output_metadata['filtered_bg']

        logger.progress("MACS2", status="RUNNING")
        m_results_files, m_results_meta = macs_caller.run(
            macs_inputs, macs_metadt,
            # Outputs of the final step may match workflow outputs;
            # Extra entries in output_files will be disregarded.
            remap(
                output_files,
                'narrow_peak', 'summits', 'broad_peak', 'gapped_peak')
        )
        logger.progress("MACS2", status="DONE")

        if not m_results_meta:
            logger.fatal("MACS2 peak calling failed")

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
        if 'control_lambda' in m_results_meta:
            output_files_generated['control_lambda'] = m_results_files['control_lambda']
            output_metadata['control_lambda'] = m_results_meta['control_lambda']

            tool_name = output_metadata['control_lambda'].meta_data['tool']
            output_metadata['control_lambda'].meta_data['tool_description'] = tool_name
            output_metadata['control_lambda'].meta_data['tool'] = "process_chipseq"
        if 'treat_pileup' in m_results_meta:
            output_files_generated['treat_pileup'] = m_results_files['treat_pileup']
            output_metadata['treat_pileup'] = m_results_meta['treat_pileup']

            tool_name = output_metadata['treat_pileup'].meta_data['tool']
            output_metadata['treat_pileup'].meta_data['tool_description'] = tool_name
            output_metadata['treat_pileup'].meta_data['tool'] = "process_chipseq"

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

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="ChIP-seq peak calling")
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
