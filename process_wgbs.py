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

from tool.bs_seeker_aligner import bssAlignerTool
from tool.bs_seeker_filter import filterReadsTool
from tool.bs_seeker_indexer import bssIndexerTool
from tool.bs_seeker_methylation_caller import bssMethylationCallerTool


# ------------------------------------------------------------------------------

class process_wgbs(Workflow):
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

            genome_fa : str
                Genome assembly in FASTA

            fastq1 : str
                Location for the first FASTQ file for single or paired end
                reads

            fastq2 : str
                [OPTIONAL]Location for the second FASTQ file if paired end
                reads


        metadata : dict
            Input file meta data associated with their roles

            genome_fa : str
            fastq1 : str

            fastq2 : str
                [OPTIONAL]

        output_files : dict
            index : str
            fastq1_filtered : str

            fastq2_filtered : str
                [OPTIONAL]

            bam : str
            bai : str
            wig_file : str
            cgmap_file : str
            atcgmap_file : str

        Returns
        -------

        fastq1_filtered|fastq1_filtered : str
            Locations of the filtered FASTQ files from which alignments were
            made

        bam|bai : str
            Location of the alignment bam file and the associated index

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

        logger.info("WGBS - BS-Seeker2 Index")
        # Build the matching WGBS genome index
        builder = bssIndexerTool(self.configuration)
        logger.progress("BSseeker2 Indexer", status="RUNNING")
        genome_idx, gidx_meta = builder.run(
            remap(input_files, "genome"),
            remap(metadata, "genome"),
            remap(output_files, "index")
        )
        logger.progress("BSseeker2 Indexer", status="DONE")
        output_results_files["index"] = genome_idx["index"]
        output_metadata["index"] = gidx_meta["index"]

        # Filter the FASTQ reads to remove duplicates
        logger.info("WGBS - Filter")
        frt = filterReadsTool(self.configuration)
        logger.progress("BSseeker2 Filtering FASTQ 1", status="RUNNING")
        fastq1f, filter1_meta = frt.run(
            {"fastq": input_files["fastq1"]},
            {"fastq": metadata["fastq1"]},
            {"fastq_filtered": output_files["fastq1_filtered"]}
        )
        logger.progress("BSseeker2 Filter FASTQ 1", status="DONE")

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
            logger.progress("BSseeker2 Filtering FASTQ 2", status="RUNNING")
            fastq2f, filter2_meta = frt.run(
                {"fastq": input_files["fastq2"]},
                {"fastq": metadata["fastq2"]},
                {"fastq_filtered": output_files["fastq2_filtered"]}
            )
            logger.progress("BSseeker2 Filtering FASTQ 2", status="DONE")

            try:
                output_results_files["fastq2_filtered"] = fastq2f["fastq_filtered"]
                output_metadata["fastq2_filtered"] = filter2_meta["fastq_filtered"]

                tool_name = output_metadata["fastq2_filtered"].meta_data["tool"]
                output_metadata["fastq2_filtered"].meta_data["tool_description"] = tool_name
                output_metadata["fastq2_filtered"].meta_data["tool"] = "process_wgbs"
            except KeyError:
                logger.fatal(
                    "WGBS - FILTER (background): Error while filtering")
                return {}, {}

        logger.info("WGBS - BS-Seeker2 Aligner")
        # Handles the alignment of all of the split packets then merges them
        # back together.
        bss_aligner = bssAlignerTool(self.configuration)
        aligner_input_files = {
            "genome": input_files["genome"],
            "fastq1": fastq1f["fastq_filtered"]
        }
        aligner_input_files["index"] = genome_idx["index"]

        aligner_meta = {
            "genome": metadata["genome"],
            "fastq1": filter1_meta["fastq_filtered"],
            "index": output_metadata["index"]
        }
        if "fastq2" in input_files:
            aligner_input_files["fastq2"] = fastq2f["fastq_filtered"]
            aligner_meta["fastq2"] = filter2_meta["fastq_filtered"]

        logger.progress("BSseeker2 Aligner", status="RUNNING")
        bam, bam_meta = bss_aligner.run(
            aligner_input_files,
            aligner_meta,
            remap(output_files, "bam", "bai")
        )
        logger.progress("BSseeker2 Aligner", status="DONE")

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
            output_metadata["bai"].meta_data["tool"] = "process_wgbs"
        except KeyError:
            logger.fatal("WGBS - Aligner failed")
            return {}, {}

        # Methylation peak caller
        peak_caller_handle = bssMethylationCallerTool(self.configuration)
        mct_input_files = {
            "genome": input_files["genome"],
            "index": genome_idx["index"],
            "fastq1": fastq1f["fastq_filtered"],
            "bam": bam["bam"],
            "bai": bam["bai"]
        }

        mct_meta = {
            "genome": metadata["genome"],
            "index": gidx_meta["index"],
            "fastq1": filter1_meta["fastq_filtered"],
            "bam": output_metadata["bam"],
            "bai": bam_meta["bai"]
        }

        if "fastq2" in input_files:
            mct_input_files["fastq2"] = fastq2f["fastq_filtered"]
            mct_meta["fastq2"] = filter2_meta["fastq_filtered"]

        logger.progress("BSseeker2 Peak Caller", status="RUNNING")
        peak_files, peak_meta = peak_caller_handle.run(
            mct_input_files,
            mct_meta,
            {
                "wig_file": output_files["wig_file"],
                "atcgmap_file": output_files["atcgmap_file"],
                "cgmap_file": output_files["cgmap_file"]
            }
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
            output_metadata["wig_file"].meta_data["tool"] = "process_wgbs"
            output_metadata["cgmap_file"].meta_data["tool_description"] = output_metadata["cgmap_file"].meta_data["tool"]  # pylint: disable=line-too-long
            output_metadata["cgmap_file"].meta_data["tool"] = "process_wgbs"
            output_metadata["atcgmap_file"].meta_data["tool_description"] = output_metadata["atcgmap_file"].meta_data["tool"]  # pylint: disable=line-too-long
            output_metadata["atcgmap_file"].meta_data["tool"] = "process_wgbs"
        except KeyError:
            logger.fatal("WGBS - Peak caller failed")
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
    result = app.launch(process_wgbs,
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
    PARSER = argparse.ArgumentParser(description="WGBS peak calling")
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
