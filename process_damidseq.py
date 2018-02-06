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

from tool.forge_bsgenome import bsgenomeTool
from tool.bwa_aligner import bwaAlignerTool
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
        output_files_generated = {}
        output_metadata = {}

        # Add in BSgenome section

        logger.info("PROCESS DAMIDSEQ - DEFINED OUTPUT:", output_files)

        alignment_set = [
            ["fastq_1", "bam_1", "bam_1_filtered"],
            ["fastq_2", "bam_2", "bam_2_filtered"],
            ["bg_fastq_1", "bg_bam_1", "bg_bam_1_filtered"],
            ["bg_fastq_2", "bg_bam_2", "bg_bam_2_filtered"],
        ]

        # BSgenome
        logger.info("Generating BSgenome")
        bsg = bsgenomeTool(self.configuration)
        bsgi, bsgm = bsg.run(
            {"genome": input_files["genome"]},
            {"genome": metadata["genome"]},
            {
                "bsgenome": output_files["bsgenome"],
                "chrom_size": output_files["chrom_size"],
                "genome_2bit": output_files["genome_2bit"],
                "seed_file": output_files["seed_file"]
            }
        )

        try:
            for file_key in ["bsgenome", "chrom_size", "genome_2bit", "seed_file"]:
                output_files_generated[file_key] = bsgi[file_key]
                output_metadata[file_key] = bsgm[file_key]

                tool_name = output_metadata[file_key].meta_data['tool']
                output_metadata[file_key].meta_data['tool_description'] = tool_name
                output_metadata[file_key].meta_data['tool'] = "process_damidseq"
        except KeyError:
            logger.fatal("BSgenome indexer failed")

        # Align and filter reads
        for aln in alignment_set:
            bwa = bwaAlignerTool(self.configuration)
            bwa_files, bwa_meta = bwa.run(
                remap(input_files,
                      "genome", "index", loc=aln[0]),
                remap(metadata,
                      "genome", "index", loc=aln[0]),
                {"output": output_files[aln[1]]}
            )

            try:
                output_files_generated[aln[1]] = bwa_files["bam"]
                output_metadata[aln[1]] = bwa_meta["bam"]

                tool_name = output_metadata[aln[1]].meta_data["tool"]
                output_metadata[aln[1]].meta_data["tool_description"] = tool_name
                output_metadata[aln[1]].meta_data["tool"] = "process_damidseq"
            except KeyError as msg:
                logger.fatal(
                    "KeyError error - BWA aligner failed: {0}\n{1}\n{2}\n{3}".format(
                        msg, aln[1],
                        "Available file keys: " + ", ".join(bwa_files.keys()),
                        "Available mets keys: " + ", ".join(bwa_meta.keys())
                    )
                )

            # Filter the bams
            b3f = biobambam(self.configuration)
            b3f_files, b3f_meta = b3f.run(
                {"input": bwa_files["bam"]},
                {"input": bwa_meta["bam"]},
                {"output": output_files[aln[2]]}
            )

            try:
                output_files_generated[aln[2]] = b3f_files["bam"]
                output_metadata[aln[2]] = b3f_meta["bam"]

                tool_name = output_metadata[aln[2]].meta_data["tool"]
                output_metadata[aln[2]].meta_data["tool_description"] = tool_name
                output_metadata[aln[2]].meta_data["tool"] = "process_damidseq"
            except KeyError as msg:
                logger.fatal("KeyError error - BioBamBam filtering failed: {0}\n{1}".format(
                    msg, aln[2]))

                return {}, {}

        ## iDEAR to call peaks
        idear_caller = idearTool(self.configuration)
        idear_files, idear_meta = idear_caller.run(
            {
                "bam_1" : output_files_generated["bam_1_filtered"],
                "bam_2" : output_files_generated["bam_2_filtered"],
                "bg_bam_1" : output_files_generated["bg_bam_1_filtered"],
                "bg_bam_2" : output_files_generated["bg_bam_2_filtered"],
                "bsgenome" : input_files["bsgenome"]
            }, {
                "bam_1" : output_metadata["bam_1_filtered"],
                "bam_2" : output_metadata["bam_2_filtered"],
                "bg_bam_1" : output_metadata["bg_bam_1_filtered"],
                "bg_bam_2" : output_metadata["bg_bam_2_filtered"],
                "bsgenome" : metadata["bsgenome"]
            }, {
                "bigwig" : output_files["bigwig"],
            }
        )

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
