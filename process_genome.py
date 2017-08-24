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

from basic_modules.tool import Tool
from basic_modules.workflow import Workflow
from basic_modules.metadata import Metadata

from dmp import dmp

from tool.bowtie_indexer import bowtieIndexerTool
from tool.bwa_indexer import bwaIndexerTool
from tool.gem_indexer import gemIndexerTool

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
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        The downloading can be done using the current common.py functions. These
        should be prevented from running the indexing step as this will be done
        as part of this workflow.

        Parameters
        ----------
        input_files : list
            List of file locations
        metadata : list
            Required meta data
        output_files : list
            List of output file locations

        Returns
        -------
        outputfiles : list
            List of locations for the output index files
        """

        genome_fa = input_files[0]
        output_metadata = {}
        output_metadata['genome_idx'] = {}

        # Bowtie2 Indexer
        bowtie2 = bowtieIndexerTool()
        bti, btm = bowtie2.run([genome_fa], [], {})
        output_metadata['genome_idx']['bowtie'] = btm

        # BWA Indexer
        bwa = bwaIndexerTool()
        bwai, bwam = bwa.run([genome_fa], [], {})
        output_metadata['genome_idx']['bwa'] = bwam

        # GEM Indexer
        gem = gemIndexerTool()
        gemi, gemm = gem.run([genome_fa], [], {})
        output_metadata['genome_idx']['gem'] = gemm

        # Need to get the file_id for the FASTA file to use as source_id
        # Also need to get the accession of the assembly from the parent
        output_metadata['output_files'] = [
            Metadata("index", "Assembly"),
            Metadata("index", "Assembly"),
            Metadata("index", "Assembly"),
        ]

        return (bti + bwai + gemi, output_metadata)

# ------------------------------------------------------------------------------

def main(input_files, output_files, input_metadata):
    """
    Main function
    -------------

    This function launches the app.
    """

    # import pprint  # Pretty print - module for dictionary fancy printing

    # 1. Instantiate and launch the App
    print("1. Instantiate and launch the App")
    from apps.workflowapp import WorkflowApp
    app = WorkflowApp()
    result = app.launch(process_genome, input_files, input_metadata,
                        output_files, {})

    # 2. The App has finished
    print("2. Execution finished")
    print(result)
    return result

def prepare_files(
        dm_handler, taxon_id, genome_fa, assembly):
    """
    Function to load the DM API with the required files and prepare the
    parameters passed from teh command line ready for use in the main function
    """
    print(dm_handler.get_files_by_user("test"))

    genome_file = dm_handler.set_file(
        "test", genome_fa, "fasta", "Assembly", taxon_id, None, [],
        meta_data={"assembly" : assembly})

    # Maybe it is necessary to prepare a metadata parser from json file
    # when building the Metadata objects.
    metadata = [
        Metadata("fasta", "Assembly", None, {'assembly' : assembly}, genome_file),
    ]

    files = [
        genome_fa,
    ]

    files_out = []

    return [files, files_out, metadata]

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    sys._run_from_cmdl = True

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="Index the genome file")
    PARSER.add_argument("--taxon_id", help="Species (9606)")
    PARSER.add_argument("--genome", help="Genome FASTA file")
    PARSER.add_argument("--assembly", help="Assembly ID")

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    GENOME_FA = ARGS.genome
    ASSEMBLY = ARGS.assembly
    TAXON_ID = ARGS.taxon_id

    #
    # MuG Tool Steps
    # --------------
    #
    # 1. Create data files
    DM_HANDLER = dmp(test=True)

    # Get the assembly

    #2. Register the data with the DMP
    PARAMS = prepare_files(DM_HANDLER, TAXON_ID, GENOME_FA, ASSEMBLY)

    RESULTS = main(PARAMS[0], PARAMS[1], PARAMS[2])

    print(DM_HANDLER.get_files_by_user("test"))
