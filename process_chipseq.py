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

# Required for ReadTheDocs
from functools import wraps # pylint: disable=unused-import

import argparse

from basic_modules.workflow import Workflow
from basic_modules.metadata import Metadata

from dmp import dmp

from tool.bwa_aligner import bwaAlignerTool
from tool.biobambam_filter import biobambam
from tool.macs2 import macs2

# ------------------------------------------------------------------------------

class process_chipseq(Workflow):
    """
    Functions for processing Chip-Seq FastQ files. Files are the aligned,
    filtered and analysed for peak calling
    """

    #configuration = {}

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


    def run(self, file_ids, metadata, output_files):
        """
        Main run function for processing ChIP-seq FastQ data. Pipeline aligns
        the FASTQ files to the genome using BWA. MACS 2 is then used for peak
        calling to identify transcription factor binding sites within the
        genome.

        Currently this can only handle a single data file and a single
        background file.

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

        run_genome_fa = file_ids[0]
        bwa_amb = file_ids[1]
        bwa_ann = file_ids[2]
        bwa_bwt = file_ids[3]
        bwa_pac = file_ids[4]
        bwa_sa = file_ids[5]
        file_loc = file_ids[6]
        file_bgd_loc = file_ids[7]

        out_bam = file_loc.replace(".fastq", '.bam')

        bwa = bwaAlignerTool(self.configuration)
        bwa_results = bwa.run(
            [run_genome_fa, file_loc, bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_sa],
            []
        )
        out_bam = bwa_results[0][0]
        
        out_bgd_bam = None
        if file_bgd_loc != None:
            bwa_results_bgd = bwa.run(
                [run_genome_fa, file_bgd_loc, bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_sa],
                []
            )
            out_bgd_bam = bwa_results_bgd[0][0]

        # For multiple files there will need to be merging into a single bam file

        # Filter the bams
        b3f = biobambam(self.configuration)
        b3f_results = b3f.run(
            [out_bam],
            []
        )
        out_filtered_bam = b3f_results[0][0]

        out_filtered_bgd_bam = None
        if file_bgd_loc != None:
            b3f_results_bgd = b3f.run(
                [out_bgd_bam],
                []
            )
            out_filtered_bgd_bam = b3f_results_bgd[0][0]

        ## MACS2 to call peaks
        macs_caller = macs2(self.configuration)

        if file_bgd_loc != None:
            m_results_files, m_results_meta = macs_caller.run([out_filtered_bam, out_filtered_bgd_bam], [])
        else:
            m_results_files, m_results_meta = macs_caller.run([out_filtered_bam], [])

        return (
            [out_bam, out_filtered_bam] + m_results_files,
            [{}, m_results_meta]
        )

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
    result = app.launch(process_chipseq, input_files, input_metadata, output_files, {})

    # 2. The App has finished
    print("2. Execution finished")
    print(result)
    return result

def prepare_files(
        dm_handler, taxon_id, genome_fa, assembly, file_loc, file_bg_loc=None):
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

    # Maybe it is necessary to prepare a metadata parser from json file
    # when building the Metadata objects.
    metadata = [
        Metadata("fasta", "Assembly"),
        Metadata("index", "Assembly"),
        Metadata("index", "Assembly"),
        Metadata("index", "Assembly"),
        Metadata("index", "Assembly"),
        Metadata("index", "Assembly"),
        Metadata("fasta", "ChIP-seq")
    ]

    dm_handler.set_file(
        "test", file_loc, "fasta", "ChIP-seq", taxon_id, None, [],
        meta_data={'assembly' : assembly})
    if file_bg_loc:
        dm_handler.set_file(
            "test", file_bg_loc, "fasta", "ChIP-seq", taxon_id, None, [],
            meta_data={'assembly' : assembly})
        metadata.append(Metadata("fasta", "ChIP-seq"))

    print(dm_handler.get_files_by_user("test"))

    files = [
        genome_fa,
        genome_fa + ".amb",
        genome_fa + ".ann",
        genome_fa + ".bwt",
        genome_fa + ".pac",
        genome_fa + ".sa",
        file_loc,
        file_bg_loc
    ]

    out_bam = file_loc.replace(".fastq", '.filtered.bam')
    out_peaks_narrow = file_loc.replace(".fastq", '_peaks.narrowPeak')
    out_peaks_xls = file_loc.replace(".fastq", '_peaks.xls')
    out_peaks_broad = file_loc.replace(".fastq", '_summits.bed')
    out_peaks_gapped = file_loc.replace(".fastq", '_summits.bed')
    out_summits = file_loc.replace(".fastq", '_summits.bed')


    files_out = [
        out_bam,
        out_peaks_narrow,
        out_peaks_xls,
        out_peaks_broad,
        out_peaks_gapped,
        out_summits
    ]
    return [files, files_out, metadata]

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    sys._run_from_cmdl = True

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="ChIP-seq peak calling")
    PARSER.add_argument("--taxon_id", help="Taxon_ID (9606)")
    PARSER.add_argument("--genome", help="Genome FASTA file")
    PARSER.add_argument("--assembly", help="Genome assembly ID (GCA_000001405.25)")
    PARSER.add_argument("--file", help="Location of FASTQ input file")
    PARSER.add_argument("--bgd_file", help="Location of FASTQ background file", default=None)

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    TAXON_ID = ARGS.taxon_id
    GENOME_FA = ARGS.genome
    ASSEMBLY = ARGS.assembly
    FILE_LOC = ARGS.file
    FILE_BG_LOC = ARGS.bgd_file

    #
    # MuG Tool Steps
    # --------------
    #
    # 1. Create data files
    DM_HANDLER = dmp(test=True)

    #2. Register the data with the DMP
    PARAMS = prepare_files(DM_HANDLER, TAXON_ID, GENOME_FA, ASSEMBLY, FILE_LOC, FILE_BG_LOC)

    # 3. Instantiate and launch the App
    RESULTS = main(PARAMS[0], PARAMS[1], PARAMS[2])

    print(RESULTS)
    print(DM_HANDLER.get_files_by_user("test"))
