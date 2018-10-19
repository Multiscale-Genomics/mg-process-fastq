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
import os.path

from basic_modules.workflow import Workflow
from utils import logger
from utils import remap

from tool.forge_bsgenome import bsgenomeTool
from tool.bwa_mem_aligner import bwaAlignerMEMTool
from tool.biobambam_filter import biobambam
from tool.idear import idearTool


# ------------------------------------------------------------------------------

class process_damidseq(Workflow):
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
        logger.info("Processing DamID-Seq")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def _align_filter(self, align_input_files, align_input_file_meta, output_files):
        """
        Function for performing the alignment and filtering of fastq files.
        """

        output_files_generated = {}
        output_metadata_generated = {}

        bwa = bwaAlignerMEMTool(self.configuration)
        logger.progress("BWA MEM Aligner - " + align_input_files["loc"], status="RUNNING")
        bwa_files, bwa_meta = bwa.run(
            align_input_files, align_input_file_meta, {"output": output_files["bam"]}
        )
        logger.progress("BWA MEM Aligner - " + align_input_files["loc"], status="DONE")

        try:
            output_files_generated["bam"] = bwa_files["bam"]
            output_metadata_generated["bam"] = bwa_meta["bam"]

            tool_name = output_metadata_generated["bam"].meta_data["tool"]
            output_metadata_generated["bam"].meta_data["tool_description"] = tool_name
            output_metadata_generated["bam"].meta_data["tool"] = "process_damidseq"
        except KeyError as msg:
            logger.fatal(
                "KeyError error - BWA aligner failed: {0}\n{1}\n{2}\n{3}".format(
                    msg, output_files_generated["bam"],
                    "Available file keys: " + ", ".join(bwa_files.keys()),
                    "Available mets keys: " + ", ".join(bwa_meta.keys())
                )
            )
            return {}, {}

        # Filter the bams
        b3f = biobambam(self.configuration)
        logger.progress("BioBamBam Filtering - " + align_input_files["loc"], status="RUNNING")
        b3f_files, b3f_meta = b3f.run(
            {"input": bwa_files["bam"]},
            {"input": bwa_meta["bam"]},
            {"output": output_files["bam_filtered"]}
        )
        logger.progress("BioBamBam Filtering - " + align_input_files["loc"], status="DONE")

        try:
            output_files_generated["bam_filtered"] = b3f_files["output"]
            output_metadata_generated["bam_filtered"] = b3f_meta["output"]

            tool_name = output_metadata_generated["bam_filtered"].meta_data["tool"]
            output_metadata_generated["bam_filtered"].meta_data["tool_description"] = tool_name
            output_metadata_generated["bam_filtered"].meta_data["tool"] = "process_damidseq"
        except KeyError as msg:
            logger.fatal("KeyError error - BioBamBam filtering failed: {0}\n{1}".format(
                msg, output_files_generated["bam_filtered"]))
            return {}, {}

        return (output_files_generated, output_metadata_generated)

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

            fastq_1 : str
                Location of the FASTQ reads files

            fastq_2 : str
                Location of the FASTQ repeat reads files

            bg_fastq_1 : str
                Location of the background FASTQ reads files

            bg_fastq_2 : str
                Location of the background FASTQ repeat reads files

        metadata : dict
            Input file meta data associated with their roles

            genome : str
            index : str
            fastq_1 : str
            fastq_2 : str
            bg_fastq_1 : str
            bg_fastq_2 : str

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
            bigwig : str
                Location of the bigwig peaks

        output_metadata : dict
            Output metadata for the associated files in output_files

            bam [, "bam_bg"] : Metadata
            filtered [, "filtered_bg"] : Metadata
            bigwig : Metadata

        """
        output_files_generated = {
            "bam": [],
            "bam_filtered": []
        }
        output_metadata = {
            "bam": [],
            "bam_filtered": []
        }

        # BSgenome
        logger.info("Generating BSgenome")

        if "genome_public" in input_files:
            genome_input_file = {"genome": input_files["genome_public"]}
            genome_input_meta = {"genome": metadata["genome_public"]}
        else:
            genome_input_file = {"genome": input_files["genome"]}
            genome_input_meta = {"genome": metadata["genome"]}

        bsg = bsgenomeTool(self.configuration)
        logger.progress("BSgenome Indexer", status="RUNNING")
        bsgi, bsgm = bsg.run(
            genome_input_file,
            genome_input_meta,
            {
                "bsgenome": output_files["bsgenome"],
                "chrom_size": output_files["chrom_size"],
                "genome_2bit": output_files["genome_2bit"],
                "seed_file": output_files["seed_file"]
            }
        )
        logger.progress("BSgenome Indexer", status="DONE")

        try:
            file_keys = ["bsgenome", "chrom_size", "genome_2bit", "seed_file"]
            for file_key in file_keys:
                output_files_generated[file_key] = bsgi[file_key]
                output_metadata[file_key] = bsgm[file_key]

                tool_name = output_metadata[file_key].meta_data['tool']
                output_metadata[file_key].meta_data['tool_description'] = tool_name
                output_metadata[file_key].meta_data['tool'] = "process_damidseq"
        except KeyError:
            logger.fatal("BSgenome indexer failed")
            return {}, {}

        # Align and filter reads
        for prefix in ["", "bg_"]:
            for i, aln in enumerate(input_files[prefix + "fastq_1"]):
                logger.info("BWA MEM Aligning and filtering of " + aln)
                if "genome_public" in input_files:
                    align_input_files = remap(
                        input_files, genome="genome_public", index="index_public",
                        loc=input_files[prefix + "fastq_1"][i])
                    align_input_file_meta = remap(
                        metadata, genome="genome_public", index="index_public",
                        loc=input_files[prefix + "fastq_1"][i])
                else:
                    align_input_files = remap(
                        input_files, genome="genome", index="index",
                        loc=input_files[prefix + "fastq_1"][i])
                    align_input_file_meta = remap(
                        metadata, genome="genome", index="index",
                        loc=input_files[prefix + "fastq_1"][i])

                if prefix + "fastq_2" in input_files:
                    align_input_files["fastq_2"] = input_files[prefix + "fastq_2"][i]
                    align_input_file_meta["fastq_2"] = metadata[prefix + "fastq_2"][i]

                fastq_in = os.path.split(input_files[prefix + "fastq_1"][i])
                fastq_suffix = fastq_in[1].split(".")[-1]
                align_output_files = {
                    prefix + "bam": os.path.join(
                        self.configuration["execution"],
                        fastq_in[1].replace(fastq_suffix, "bam")
                    ),
                    prefix + "bam_filtered": os.path.join(
                        self.configuration["execution"],
                        fastq_in[1].replace(fastq_suffix, "filtered.bam")
                    ),
                }

                bwa_files, bwa_meta = self._align_filter(
                    align_input_files, align_input_file_meta, align_output_files)

                try:
                    output_files_generated[prefix + "bam"].append(bwa_files["bam"])
                    output_metadata[prefix + "bam"].append(bwa_meta["bam"])

                    output_files_generated[prefix + "bam_filtered"].append(
                        bwa_files["bam_filtered"])
                    output_metadata[prefix + "bam_filtered"].append(bwa_meta["bam"])
                except KeyError as msg:
                    logger.fatal("Error aligning and filtering input FASTQ files")
                    return {}, {}

        # iDEAR to call peaks
        idear_caller = idearTool(self.configuration)
        logger.progress("iDEAR Peak Caller", status="RUNNING")
        idear_files, idear_meta = idear_caller.run(
            {
                "bam": output_files_generated["bam_filtered"],
                "bg_bam": output_files_generated["bg_bam_filtered"],
                "bsgenome": input_files["bsgenome"]
            }, {
                "bam": output_metadata["bam_filtered"],
                "bg_bam": output_metadata["bg_bam_filtered"],
                "bsgenome": metadata["bsgenome"]
            }, {
                "bigwig": output_files["bigwig"],
            }
        )
        logger.progress("iDEAR Peak Caller", status="DONE")

        try:
            output_files_generated["bigwig"] = idear_files["bigwig"]
            output_metadata["bigwig"] = idear_meta["bigwig"]

            tool_name = output_metadata["bigwig"].meta_data["tool"]
            output_metadata["bigwig"].meta_data["tool_description"] = tool_name
            output_metadata["bigwig"].meta_data["tool"] = "process_damidseq"
        except KeyError as msg:
            logger.fatal("KeyError error - iDEAR filtering failed: {0}\n{1}".format(
                msg, "bigwig"))

            return {}, {}

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
    result = app.launch(process_damidseq,
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
    PARSER = argparse.ArgumentParser(description="iDamID-seq peak calling")
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
