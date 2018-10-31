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
from utils import remap

from mg_process_fastq.tool.bwa_aligner import bwaAlignerTool
from mg_process_fastq.tool.inps import inps


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

        if "genome_public" in input_files:
            align_input_files = remap(
                input_files, genome="genome_public", loc="loc", index="index_public")
            align_input_file_meta = remap(
                metadata, genome="genome_public", loc="loc", index="index_public")
        else:
            align_input_files = remap(input_files, "genome", "loc", "index")
            align_input_file_meta = remap(metadata, "genome", "loc", "index")

        bwa = bwaAlignerTool()
        logger.progress("BWA ALN Aligner", status="RUNNING")
        bwa_files, bwa_meta = bwa.run(
            align_input_files, align_input_file_meta,
            {"output": output_files["bam"], "bai": output_files["bai"]}
        )
        logger.progress("BWA ALN Aligner", status="DONE")

        output_files_generated = {}

        try:
            output_files_generated["bam"] = bwa_files["bam"]
            output_metadata["bam"] = bwa_meta["bam"]

            tool_name = output_metadata['bam'].meta_data['tool']
            output_metadata['bam'].meta_data['tool_description'] = tool_name
            output_metadata['bam'].meta_data['tool'] = "process_mnaseseq"

            output_files_generated["bai"] = bwa_files["bai"]
            output_metadata["bai"] = bwa_meta["bai"]

            tool_name = output_metadata['bai'].meta_data['tool']
            output_metadata['bai'].meta_data['tool_description'] = tool_name
            output_metadata['bai'].meta_data['tool'] = "process_mnaseseq"
        except KeyError:
            logger.fatal("BWA Alignment failed")

        inps_tool = inps()
        logger.progress("iNPS Peak Caller", status="RUNNING")
        inps_files, inps_meta = inps_tool.run(
            remap(bwa_files, "bam"),
            remap(bwa_meta, "bam"),
            {"bed": output_files["bed"]}
        )
        logger.progress("iNPS Peak Caller", status="DONE")

        try:
            output_files_generated["bed"] = inps_files["bed"]
            output_metadata["bed"] = inps_meta["bed"]

            tool_name = output_metadata['bed'].meta_data['tool']
            output_metadata['bed'].meta_data['tool_description'] = tool_name
            output_metadata['bed'].meta_data['tool'] = "process_mnaseseq"
        except KeyError:
            logger.fatal("BWA Alignment failed")

        print("MNASESEQ RESULTS:", output_metadata)
        return output_files, output_metadata
