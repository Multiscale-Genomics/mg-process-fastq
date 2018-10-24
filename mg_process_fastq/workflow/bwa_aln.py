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

from basic_modules.workflow import Workflow
from utils import logger

from mg_process_fastq.tool.bwa_aligner import bwaAlignerTool


# ------------------------------------------------------------------------------

class process_bwa(Workflow):
    """
    Functions for aligning FastQ files with BWA ALN
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
        logger.info("Processing BWA Aligner")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        Main run function for aligning FastQ reads with BWA ALN.

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
                [OPTIONAL] Location of the FASTQ reads file for paired end data

        metadata : dict
            Input file meta data associated with their roles

            genome : str
            index : str
            loc : str
            fastq2 : str
        output_files : dict
            Output file locations

            bam : str
                Output bam file location

        Returns
        -------
        output_files : dict
            Output file locations associated with their roles, for the output

            bam : str
                Aligned FASTQ short read file locations
        output_metadata : dict
            Output metadata for the associated files in output_files

            bam : Metadata
        """
        output_files_generated = {}
        output_metadata = {}

        logger.info("PROCESS ALIGNMENT - DEFINED OUTPUT:", output_files["bam"])

        if "genome_public" in input_files:
            input_files["genome"] = input_files.pop("genome_public")
            metadata["genome"] = metadata.pop("genome_public")

        if "index_public" in input_files:
            input_files["index"] = input_files.pop("index_public")
            metadata["index"] = metadata.pop("index_public")

        bwa = bwaAlignerTool(self.configuration)

        logger.progress("BWA ALN Aligner", status="RUNNING")
        bwa_files, bwa_meta = bwa.run(
            input_files, metadata,
            {"output": output_files["bam"], "bai": output_files["bai"]}
        )
        logger.progress("BWA ALN Aligner", status="DONE")

        try:
            output_files_generated["bam"] = bwa_files["bam"]
            output_metadata["bam"] = bwa_meta["bam"]

            tool_name = output_metadata['bam'].meta_data['tool']
            output_metadata['bam'].meta_data['tool_description'] = tool_name
            output_metadata['bam'].meta_data['tool'] = "process_bwa_aln"

            output_files_generated["bai"] = bwa_files["bai"]
            output_metadata["bai"] = bwa_meta["bai"]

            tool_name = output_metadata['bai'].meta_data['tool']
            output_metadata['bai'].meta_data['tool_description'] = tool_name
            output_metadata['bai'].meta_data['tool'] = "process_bwa_aln"
        except KeyError:
            logger.fatal("BWA aligner failed")

        return output_files_generated, output_metadata
