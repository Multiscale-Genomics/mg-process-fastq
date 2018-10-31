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

from mg_process_fastq.tool.bowtie_indexer import bowtieIndexerTool
from mg_process_fastq.tool.bwa_indexer import bwaIndexerTool
from mg_process_fastq.tool.gem_indexer import gemIndexerTool


# ------------------------------------------------------------------------------

class process_genome(Workflow):
    """
    Workflow to download and pre-index a given genome
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
        logger.info("Processing Genomes")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        Main run function for the indexing of genome assembly FASTA files. The
        pipeline uses Bowtie2, BWA and GEM ready for use in pipelines that
        rely on alignment.

        Parameters
        ----------
        input_files : dict
            genome : str
                List of file locations
        metadata : dict
            genome : dict
                Required meta data
        output_files : dict
            bwa_index : str
                Location of the BWA index archive files
            bwt_index : str
                Location of the Bowtie2 index archive file
            gem_index : str
                Location of the GEM index file
            genome_gem : str
                Location of a the FASTA file generated for the GEM indexing
                step

        Returns
        -------
        outputfiles : dict
            List of locations for the output index files
        output_metadata : dict
            Metadata about each of the files
        """

        output_files_generated = {}
        output_metadata = {}

        if "genome_public" in input_files:
            genome_input_file = {"genome": input_files["genome_public"]}
            genome_input_meta = {"genome": metadata["genome_public"]}
        else:
            genome_input_file = {"genome": input_files["genome"]}
            genome_input_meta = {"genome": metadata["genome"]}

        # Bowtie2 Indexer
        logger.info("Generating indexes for Bowtie2")
        bowtie2 = bowtieIndexerTool()
        logger.progress("Bowtie2 Indexer", status="RUNNING")
        bti, btm = bowtie2.run(
            genome_input_file,
            genome_input_meta,
            {'index': output_files['bwt_index']}
        )
        logger.progress("Bowtie2 Indexer", status="DONE")

        try:
            output_files_generated['bwt_index'] = bti["index"]
            output_metadata['bwt_index'] = btm['index']

            tool_name = output_metadata['bwt_index'].meta_data['tool']
            output_metadata['bwt_index'].meta_data['tool_description'] = tool_name
            output_metadata['bwt_index'].meta_data['tool'] = "process_genome"
        except KeyError:
            logger.fatal("Bowtie2 indexer failed")

        # BWA Indexer
        logger.info("Generating indexes for BWA")
        bwa = bwaIndexerTool()
        logger.progress("BWA Indexer", status="RUNNING")
        bwai, bwam = bwa.run(
            genome_input_file,
            genome_input_meta,
            {'index': output_files['bwa_index']}
        )
        logger.progress("BWA Indexer", status="DONE")

        try:
            output_files_generated['bwa_index'] = bwai['index']
            output_metadata['bwa_index'] = bwam['index']

            tool_name = output_metadata['bwa_index'].meta_data['tool']
            output_metadata['bwa_index'].meta_data['tool_description'] = tool_name
            output_metadata['bwa_index'].meta_data['tool'] = "process_genome"
        except KeyError:
            logger.fatal("BWA indexer failed")

        # GEM Indexer
        logger.info("Generating indexes for GEM")
        gem = gemIndexerTool()
        logger.progress("GEM Indexer", status="RUNNING")
        gemi, gemm = gem.run(
            genome_input_file, genome_input_meta,
            {
                'index': output_files['gem_index']
            }
        )
        logger.progress("GEM Indexer", status="DONE")

        try:
            output_files_generated['gem_index'] = gemi['index']
            output_metadata['gem_index'] = gemm['index']

            tool_name = output_metadata['gem_index'].meta_data['tool']
            output_metadata['gem_index'].meta_data['tool_description'] = tool_name
            output_metadata['gem_index'].meta_data['tool'] = "process_genome"
        except KeyError:
            logger.fatal("GEM indexer failed")

        return (output_files_generated, output_metadata)
