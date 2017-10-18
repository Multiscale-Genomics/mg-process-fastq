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

# -*- coding: utf-8 -*-

from __future__ import print_function

import argparse
import sys

# Required for ReadTheDocs
from functools import wraps # pylint: disable=unused-import

from basic_modules.workflow import Workflow
from dmp import dmp

from tool.tb_full_mapping import tbFullMappingTool
from tool.tb_parse_mapping import tbParseMappingTool
from tool.tb_filter import tbFilterTool
from tool.tb_generate_tads import tbGenerateTADsTool
from tool.tb_save_hdf5_matrix import tbSaveAdjacencyHDF5Tool

# ------------------------------------------------------------------------------
class process_hic(Workflow):
    """
    Functions for downloading and processing Mnase-seq FastQ files. Files are
    downloaded from the European Nucleotide Archive (ENA), then aligned,
    filtered and analysed for peak calling
    """

    configuration = {}

    def __init__(self, configuration=None):
        """
        Initialise the tool with its configuration.


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
        Main run function for processing MNase-Seq FastQ data. Pipeline aligns
        the FASTQ files to the genome using BWA. iNPS is then used for peak
        calling to identify nucleosome position sites within the genome.

        Parameters
        ----------
        files_ids : list
            List of file locations
        metadata : list
            Required meta data
        output_files : list
            List of output file locations

        Returns
        -------
        outputfiles : list
            List of locations for the output bam, bed and tsv files
        """

        genome_fa = input_files[0]
        genome_gem = input_files[1]
        assembly = metadata['assembly']
        fastq_file_1 = input_files[2]
        fastq_file_2 = input_files[3]
        enzyme_name = metadata['enzyme_name']
        resolutions = metadata['resolutions']
        window_type = metadata['window_type']
        windows1 = metadata['windows1']
        windows2 = metadata['windows2']
        normalized = metadata['normalized']
        saveas_hdf5 = metadata['hdf5']
        expt_name = metadata['expt_name']

        print("HIC - metadata:", metadata)

        input_metadata_mapping1 = {
            'windows' : windows1,
        }
        input_metadata_mapping2 = {
            'windows' : windows2,
        }

        if window_type == 'frag':
            input_metadata_mapping1['windows'] = None
            input_metadata_mapping2['windows'] = None

        if enzyme_name is not None:
            input_metadata_mapping1['enzyme_name'] = enzyme_name
            input_metadata_mapping2['enzyme_name'] = enzyme_name

        tfm1 = tbFullMappingTool()
        tfm1_files, tfm1_meta = tfm1.run([genome_gem, fastq_file_1], [], input_metadata_mapping1)

        tfm2 = tbFullMappingTool()
        tfm2_files, tfm2_meta = tfm2.run([genome_gem, fastq_file_2], [], input_metadata_mapping2)

        tpm = tbParseMappingTool()
        files = [genome_fa] + tfm1_files + tfm2_files
        input_metadata_parser = {
            'enzyme_name' : enzyme_name,
            'mapping' : [tfm1_meta['func'], tfm2_meta['func']],
            'expt_name' : expt_name
        }
        print("TB MAPPED FILES:", files)
        print("TB PARSE METADATA:", input_metadata_parser)
        tpm_files, tpm_meta = tpm.run(files, [], input_metadata_parser)

        print("TB PARSED FILES:", tpm_files)

        tbf = tbFilterTool()
        tf_files, tf_meta = tbf.run(tpm_files, [], {'conservative' : True, 'expt_name' : expt_name})

        #adjlist_loc = f2a.save_hic_data()

        print("TB FILTER FILES:", tf_files[0])

        tgt = tbGenerateTADsTool()
        tgt_meta_in = {
            'resolutions' : resolutions,
            'normalized' : False,
            'expt_name' : expt_name
        }
        tgt_files, tgt_meta = tgt.run([tf_files[0]], [], tgt_meta_in)

        # Generate the HDF5 and meta data required for the RESTful API.
        # - Chromosome meta is from the tb_parse_mapping step

        hdf5_file = None
        if saveas_hdf5 is True:
            th5 = tbSaveAdjacencyHDF5Tool()
            th5_files_in = [tf_files[0], genome_fa]
            th5_meta_in = {
                'assembly' : assembly,
                'resolutions' : resolutions,
                'normalized' : normalized,
                'chromosomes_meta' : tpm_meta['chromosomes']
            }
            th5_files, th5_meta = th5.run(th5_files_in, [], th5_meta_in)
            hdf5_file = th5_files[0]

        # List of files to get saved
        return ([tfm1_files[0], tfm2_files[0], tpm_files[0], tf_files[0], hdf5_file], [])

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
    result = app.launch(process_hic, input_files, input_metadata, output_files,
                        {})

    # 2. The App has finished
    print("2. Execution finished")
    print(result)
    return result

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    sys._run_from_cmdl = True

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="Generate adjacency files")
    PARSER.add_argument("--genome", help="Genome assembly FASTA file")
    PARSER.add_argument("--genome_gem", help="Genome assembly GEM file")
    PARSER.add_argument("--taxon_id", help="Species (9606)")
    PARSER.add_argument("--assembly", help="Assembly (GRCh38)")
    PARSER.add_argument("--file1", help="Location of FASTQ file 1")
    PARSER.add_argument("--file2", help="Location of FASTQ file 2")
    PARSER.add_argument(
        "--resolutions",
        help="CSV string of the resolutions to be computed for the models")
    PARSER.add_argument("--enzyme_name", help="Enzyme used to digest the DNA")
    PARSER.add_argument("--window_type", help="Windowing type [frag, iter]", default="frag")
    PARSER.add_argument(
        "--windows1",
        help="FASTQ windowing - start locations",
        default="1,25,50,75,100")
    PARSER.add_argument(
        "--windows2",
        help="FASTQ windowing - paired end locations",
        default="1,25,50,75,100")
    PARSER.add_argument("--normalized", help="Normalize the alignments", default=False)
    PARSER.add_argument("--tag", help="tag", default='test_name')


    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    # Assumes that there are 2 fastq files for the paired ends
    #windows1 = ((1,25), (1,50), (1,75),(1,100))
    #windows2 = ((1,25), (1,50), (1,75),(1,100))
    #windows2 = ((101,125), (101,150), (101,175),(101,200))

    GENOME_FA = ARGS.genome
    GENOME_GEM = ARGS.genome_gem
    TAXON_ID = ARGS.taxon_id
    ASSEMBLY = ARGS.assembly
    FASTQ_01_FILE = ARGS.file1
    FASTQ_02_FILE = ARGS.file2
    ENZYME_NAME = ARGS.enzyme_name
    RESOLUTIONS = ARGS.resolutions
    WINDOW_TYPE = ARGS.window_type
    WINDOWS1ARG = ARGS.windows1
    WINDOWS2ARG = ARGS.windows1
    NORMALIZED = ARGS.normalized
    EXPT_NAME = ARGS.tag

    if WINDOWS1ARG is not None:
        W1 = [int(i) for i in WINDOWS1ARG.split(",")]
        WINDOWS1 = [[W1[0], j] for j in W1[1:]]
    if WINDOWS2ARG is not None:
        W2 = [int(i) for i in WINDOWS2ARG.split(",")]
        WINDOWS2 = [[W2[0], j] for j in W2[1:]]

    print("WINDOWS1:", WINDOWS1ARG, WINDOWS1)
    print("WINDOWS2:", WINDOWS2ARG, WINDOWS2)
    print("ENZYME_NAME:", ENZYME_NAME)

    if RESOLUTIONS is None:
        # RESOLUTIONS = [
        #     1000, 2500, 5000, 10000, 25000, 50000, 100000, 250000, 500000,
        #     1000000, 10000000
        # ]
        RESOLUTIONS = [1000000, 10000000]
    else:
        RESOLUTIONS = RESOLUTIONS.split(',')

    METADATA = {
        'user_id' : 'test',
        'assembly' : ASSEMBLY,
        'resolutions' : RESOLUTIONS,
        'enzyme_name' : ENZYME_NAME,
        'windows1' : WINDOWS1,
        'windows2' : WINDOWS2,
        'normalized' : NORMALIZED,
        'hdf5' : True,
        'expt_name' : EXPT_NAME,
        'window_type' : WINDOW_TYPE
    }

    #
    # MuG Tool Steps
    # --------------
    #
    # 1. Create data files
    DM_HANDLER = dmp(test=True)

    #2. Register the data with the DMP
    genome_file = DM_HANDLER.set_file(
        "test", GENOME_FA, "fasta", "Assembly", TAXON_ID,
        meta_data={'assembly' : ASSEMBLY})
    genome_idx = DM_HANDLER.set_file(
        "test", GENOME_GEM, "gem", "Assembly Index", TAXON_ID,
        meta_data={'assembly' : ASSEMBLY})
    fastq_01_file_in = DM_HANDLER.set_file(
        "test", FASTQ_01_FILE, "fastq", "Hi-C", TAXON_ID,
        meta_data=METADATA)
    fastq_02_file_in = DM_HANDLER.set_file(
        "test", FASTQ_02_FILE, "fastq", "Hi-C", TAXON_ID,
        meta_data=METADATA)

    FILES = [
        GENOME_FA,
        GENOME_GEM,
        FASTQ_01_FILE,
        FASTQ_02_FILE
    ]

    # 3. Instantiate and launch the App
    RESULTS = main(FILES, [], METADATA)

    print(RESULTS)
    print(DM_HANDLER.get_files_by_user("test"))
