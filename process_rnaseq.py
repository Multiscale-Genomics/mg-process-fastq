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

import argparse
import os.path

from basic_modules.workflow import Workflow
from basic_modules.metadata import Metadata
from utils import remap

from tool.kallisto_indexer import kallistoIndexerTool
from tool.kallisto_quant import kallistoQuantificationTool

# ------------------------------------------------------------------------------

class process_rnaseq(Workflow):
    """
    Functions for downloading and processing RNA-seq FastQ files. Files are
    downloaded from the European Nucleotide Archive (ENA), then they are mapped
    to quantify the amount of cDNA
    """

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

        # Index the cDNA
        # This could get moved to the general tools section
        k_index = kallistoIndexerTool()
        k_out, k_meta = k_index.run(
            remap(input_files, "cdna"),
            remap(metadata, "cdna"),
            remap(output_files, "index"),
        )

        # Quantification
        k_quant = kallistoQuantificationTool()

        if "fastq2" not in input_files:
            kq_input_files = {
                "cdna": input_files["cdna"],
                "fastq1": input_files["fastq1"],
                "index": k_out["index"]
            }
            kq_input_meta = {
                "cdna": metadata["cdna"],
                "fastq1": metadata["fastq1"],
                "index": k_meta["index"]
            }

            kq_files, kq_meta = k_quant.run(
                kq_input_files,
                kq_input_meta,
                remap(output_files, "abundance_h5_file", "abundance_tsv_file", "run_info_file")
            )
        elif "fastq2" in input_files:
            kq_input_files = {
                "cdna": input_files["cdna"],
                "fastq1": input_files["fastq1"],
                "fastq2": input_files["fastq2"],
                "index": k_out["index"]
            }
            kq_input_meta = {
                "cdna": metadata["cdna"],
                "fastq1": metadata["fastq1"],
                "fastq2": metadata["fastq2"],
                "index": k_meta["index"]
            }

            kq_files, kq_meta = k_quant.run(
                kq_input_files,
                kq_input_meta,
                remap(output_files, "abundance_h5_file", "abundance_tsv_file", "run_info_file")
            )

        kq_files["index"] = k_out["index"]
        kq_meta["index"] = k_meta["index"]

        kq_meta['index'].meta_data['tool_description'] = kq_meta['index'].meta_data['tool']
        kq_meta['index'].meta_data['tool'] = "process_rnaseq"
        kq_meta['abundance_h5_file'].meta_data['tool_description'] = kq_meta['abundance_h5_file'].meta_data['tool']
        kq_meta['abundance_h5_file'].meta_data['tool'] = "process_rnaseq"
        kq_meta['abundance_tsv_file'].meta_data['tool_description'] = kq_meta['abundance_tsv_file'].meta_data['tool']
        kq_meta['abundance_tsv_file'].meta_data['tool'] = "process_rnaseq"
        kq_meta['run_info_file'].meta_data['tool_description'] = kq_meta['run_info_file'].meta_data['tool']
        kq_meta['run_info_file'].meta_data['tool'] = "process_rnaseq"

        return (kq_files, kq_meta)

# -----------------------------------------------------------------------------

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
    result = app.launch(process_rnaseq,
                        config,
                        in_metadata,
                        out_metadata)

    # 2. The App has finished
    print("2. Execution finished; see " + out_metadata)
    print(result)

    return result

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    sys._run_from_cmdl = True  # pylint: disable=protected-access

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="Parse RNA-seq for expression analysis")
    PARSER.add_argument("--config", help="Configuration file")
    PARSER.add_argument("--in_metadata", help="Location of input metadata file")
    PARSER.add_argument("--out_metadata", help="Location of output metadata file")

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    CONFIG = ARGS.config
    IN_METADATA = ARGS.in_metadata
    OUT_METADATA = ARGS.out_metadata

    RESULTS = main_json(CONFIG, IN_METADATA, OUT_METADATA)

    print(RESULTS)
