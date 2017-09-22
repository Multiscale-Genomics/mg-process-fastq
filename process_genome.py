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
import os

from basic_modules.workflow import Workflow
from basic_modules.metadata import Metadata
from utils import remap

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

        genome_fa = input_files['genome']
        output_metadata = {}
        output_metadata = {}

        # Bowtie2 Indexer
        bowtie2 = bowtieIndexerTool()
        bti, btm = bowtie2.run(input_files, metadata, {'index': output_files['bwt_index']})
        output_metadata['bwt_index'] = btm['index']

        # BWA Indexer
        bwa = bwaIndexerTool()
        bwai, bwam = bwa.run(input_files, metadata, {'index': output_files['bwa_index']})
        output_metadata['bwa_index'] = bwam['index']

        # GEM Indexer
        gem = gemIndexerTool()
        gemi, gemm = gem.run(
            input_files, metadata,
            {
                'index': output_files['gem_index'],
                'genome_gem': output_files['genome_gem']
            }
        )
        output_metadata['gem_index'] = gemm['index']
        output_metadata['genome_gem'] = gemm['genome_gem']

        return (output_files, output_metadata)

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

def main_json():
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
    root_path = os.path.dirname(__file__)
    result = app.launch(process_genome,
                        root_path,
                        "tests/json/config_genome_indexer.json",
                        "tests/json/input_genome_indexer_metadata.json")

    # 2. The App has finished
    print("2. Execution finished; see " + root_path + "/results.json")
    print(result)

    return result

def prepare_files(
        dm_handler, taxon_id, genome_fa, assembly):
    """
    Function to load the DM API with the required files and prepare the
    parameters passed from teh command line ready for use in the main function
    """
    print(dm_handler.get_files_by_user("test"))

    root_name = genome_fa.split("/")
    parent_dir = '/'.join(root_name[0:-1])

    genome_file = dm_handler.set_file(
        "test", genome_fa, "file", "fasta", 64000, parent_dir, "Assembly",
        taxon_id, None, None, meta_data={"assembly" : assembly})

    # Maybe it is necessary to prepare a metadata parser from json file
    # when building the Metadata objects.
    metadata = {
        'genome': Metadata("Assembly", "fasta", genome_fa, None,
            {'assembly' : assembly}, genome_file),
    }

    files = {
        'genome': genome_fa,
    }


    files_out = {
        'bwa_index': genome_fa + '.bwa.tar.gz',
        'bwt_index': genome_fa + '.bwt.tar.gz',
        'gem_index': genome_fa.replace('.fasta', '_gem.fasta') + '.gem.gz',
        'genome_gem': genome_fa.replace('.fasta', '_gem.fasta')
    }

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
    PARSER.add_argument("--json",
                        help="Use defined JSON config files",
                        action='store_const', const=True, default=False)

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    GENOME_FA = ARGS.genome
    ASSEMBLY = ARGS.assembly
    TAXON_ID = ARGS.taxon_id
    JSON_CONFIG = ARGS.json

    if JSON_CONFIG is True:
        RESULTS = main_json()
    else:
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

    print(RESULTS)
