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


from tool.kallisto_indexer import kallistoIndexerTool
from tool.kallisto_quant import kallistoQuantificationTool

# ------------------------------------------------------------------------------

class process_rnaseq(Workflow):
    """
    Functions for downloading and processing RNA-seq FastQ files. Files are
    downloaded from the European Nucleotide Archive (ENA), then they are mapped
    to quantify the amount of cDNA
    """

    def __init__(self, configuration={}):
        """
        Initialise the tool with its configuration.


        Parameters
        ----------
        configuration : dict
            a dictionary containing parameters that define how the operation
            should be carried out, which are specific to each Tool.
        """
        self.configuration.update(configuration)


    def run(self, file_ids, metadata, output_files):
        """
        Main run function for processing RNA-Seq FastQ data. Pipeline aligns
        the FASTQ files to the genome using Kallisto. Kallisto is then also
        used for peak calling to identify levels of expression.

        Parameters
        ----------
        files_ids : list
            List of file locations (genome FASTA, FASTQ_01, FASTQ_02 (for
            paired ends))
        metadata : list
            Required meta data
        output_files : list
            List of output file locations


        Returns
        -------
        outputfiles : list
            List of locations for the output bam, bed and tsv files
        """

        genome_fa = file_ids[0]

        # Index the cDNA
        # This could get moved to the general tools section
        k_index = kallistoIndexerTool()
        genome_idx_loc = genome_fa.replace('.fasta', '.idx')
        genome_idx_loc = genome_idx_loc.replace('.fa', '.idx')
        gi_out = k_index.run([genome_fa, genome_idx_loc], metadata)

        # Quantification
        k_quant = kallistoQuantificationTool()

        if len(file_ids) == 2:
            results = k_quant.run([genome_idx_loc, file_ids[1]], metadata)
        elif len(file_ids) == 3:
            results = k_quant.run([genome_idx_loc, file_ids[1], file_ids[2]], metadata)

        return results


# -----------------------------------------------------------------------------

def main(inputFiles, inputMetadata, outputFiles):
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
    result = app.launch(process_rnaseq, inputFiles, inputMetadata,
                        outputFiles, {})

    # 2. The App has finished
    print("2. Execution finished")

# ------------------------------------------------------------------------------


if __name__ == "__main__":
    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="Parse RNA-seq for expression analysis")
    PARSER.add_argument("--assembly", help="Genome assembly ID (GCA_000001405.25)")
    PARSER.add_argument("--taxon_id", help="Taxon_ID (9606)")
    PARSER.add_argument("--genome", help="Location of the genome cDNA FASTA file")
    PARSER.add_argument("--file", help="Location of the FASTQ file")
    PARSER.add_argument(
        "--file2",
        help="[OPTIONAL] Location of the paired end FASTQ file",
        default=None)

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    TAXON_ID = ARGS.taxon_id
    GENOME_FA = ARGS.genome
    ASSEMBLY = ARGS.assembly
    FILE_LOC = ARGS.file
    PAIRED_FILE = ARGS.file2

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

    ###########################################################

    assembly = ARGS.assembly
    taxon_id = ARGS.taxon_id
    genome_fa = ARGS.genome
    file_loc = ARGS.file
    paired_file = ARGS.file2

    da = dmp(test=True)

    print(da.get_files_by_user("test"))

    genome_file = da.set_file("test", genome_fa, "fasta", "cDNA", taxon_id, meta_data={'assembly' : assembly})
    file_in = da.set_file("test", file_loc, "fasta", "RNA-seq", taxon_id, meta_data={'assembly' : assembly})

    print(da.get_files_by_user("test"))

    file_loc_split = file_loc.split("\t")
    abundance_h5_file = "/".join(file_loc_split[0:-1]) + "/abundance.h5"
    abundance_tsv_file = "/".join(file_loc_split[0:-1]) + "/abundance.tsv"
    run_info_file = "/".join(file_loc[0:-1]) + "/run_info.json"

    file_out_0 = da.set_file("test_abundance_h5", abundance_h5_file, "hdf5", "RNA-seq", taxon_id, None, [file_in], meta_data={'assembly' : assembly})
    file_out_1 = da.set_file("test_abundance", abundance_tsv_file, "tsv", "RNA-seq", taxon_id, None, [file_in], meta_data={'assembly' : assembly})
    file_out_2 = da.set_file("test_abundance_log", run_info_file, "log", "RNA-seq", taxon_id, None, [file_in], meta_data={'assembly' : assembly})

    files_in = [
        genome_fa,
        file_loc
    ]

    files_out = [
        abundance_h5_file,
        abundance_tsv_file,
        run_info_file
    ]

    # Read metadata file and build a dictionary with the metadata:
    from basic_modules.metadata import Metadata
    # Maybe it is necessary to prepare a metadata parser from json file
    # when building the Metadata objects.
    inputMetadata = [Metadata("fasta", "RNA-seq")]

    if paired_file is not None:
        files_in.append(paired_file)
        file2_in = da.set_file("test", paired_file, "fasta", "RNA-seq", taxon_id, None, [], meta_data={'assembly' : assembly})
        inputMetadata.append(Metadata("fasta", "RNA-seq"))



    #pr = process_rnaseq()
    #results = pr.run(files, {'user_id' : 'test'})

    main(files_in, inputMetadata, files_out)

    print(results)

    print(da.get_files_by_user("test"))


