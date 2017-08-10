#!/usr/bin/env python

"""
.. Copyright 2017 EMBL-European Bioinformatics Institute

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
from basic_modules.metadata import Metadata

from dmp import dmp

from tool.bwa_aligner import bwaAlignerTool
from tool.inps import inps

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
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, file_ids, metadata, output_files):
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

        genome_file = file_ids[0]
        bwa_amb = file_ids[1]
        bwa_ann = file_ids[2]
        bwa_bwt = file_ids[3]
        bwa_pac = file_ids[4]
        bwa_sa = file_ids[5]
        file_loc = file_ids[6]

        bwa = bwaAlignerTool()
        out_file_bam, out_bam_meta = bwa.run(
            [genome_file, file_loc, bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_sa],
            []
        )

        # Needs moving to its own tool
        inps_tool = inps()
        out_peak_bed, out_meta = inps_tool.run([out_file_bam[0]], [])

        return ([out_file_bam[0], out_peak_bed], [out_bam_meta, out_meta])

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
    result = app.launch(process_mnaseseq, input_files, input_metadata,
                        output_files, {})

    # 2. The App has finished
    print("2. Execution finished")
    print(result)
    return result

def prepare_files( # pylint: disable=too-many-arguments
        dm_handler, taxon_id, genome_fa, assembly, file_loc):
    """
    Function to load the DM API with the required files and prepare the
    parameters passed from teh command line ready for use in the main function
    """
    print(dm_handler.get_files_by_user("test"))

    genome_file = dm_handler.set_file(
        "test", genome_fa, "fasta", "Assembly", taxon_id, None, [],
        meta_data={"assembly" : assembly})
    dm_handler.set_file(
        "test", genome_fa + ".amb", "amb", "Assembly", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly})
    dm_handler.set_file(
        "test", genome_fa + ".ann", "ann", "Assembly", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly})
    dm_handler.set_file(
        "test", genome_fa + ".bwt", "bwt", "Assembly", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly})
    dm_handler.set_file(
        "test", genome_fa + ".pac", "pac", "Assembly", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly})
    dm_handler.set_file(
        "test", genome_fa + ".sa", "sa", "Assembly", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly})

    in_files = [
        genome_fa,
        genome_fa + ".amb",
        genome_fa + ".ann",
        genome_fa + ".bwt",
        genome_fa + ".pac",
        genome_fa + ".sa",
        file_loc
    ]

    return [in_files, [], {}]


# ------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    sys._run_from_cmdl = True

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="Mnase-seq peak calling")
    PARSER.add_argument("--assembly", help="Genome assembly ID (GCA_000001635.2)")
    PARSER.add_argument("--taxon_id", help="Taxon ID (10090)")
    PARSER.add_argument("--genome", help="Genome assembly FASTA file")
    PARSER.add_argument("--file", help="Location of FASTQ file")

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    GENOME_FA = ARGS.genome
    TAXON_ID = ARGS.taxon_id
    ASSEMBLY = ARGS.assembly
    FILE_LOC = ARGS.file

    #
    # MuG Tool Steps
    # --------------
    #
    # 1. Create data files
    DM_HANDLER = dmp(test=True)

    #2. Register the data with the DMP
    PARAMS = prepare_files(DM_HANDLER, TAXON_ID, GENOME_FA, ASSEMBLY, FILE_LOC)

    # 3. Instantiate and launch the App
    RESULTS = main(PARAMS[0], [], [PARAMS[2]])

    print(RESULTS)
    print(DM_HANDLER.get_files_by_user("test"))
