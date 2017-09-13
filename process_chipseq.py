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
import os.path

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

    def run(self, input_files, metadata, output_files):
        """
        Main run function for processing ChIP-seq FastQ data. Pipeline aligns
        the FASTQ files to the genome using BWA. MACS 2 is then used for peak
        calling to identify transcription factor binding sites within the
        genome.

        Currently this can only handle a single data file and a single
        background file.

        Parameters
        ----------
        input_files : dict
            Location of the initial input files required by the workflow
            genome : str
                Genome FASTA file
            amb : str
                BWA index file
            ann : str
                BWA index file
            bwt : str
                BWA index file
            pac : str
                BWA index file
            sa : str
                BWA index file
            loc : str
                Location of the FASTQ reads files
            bg_loc : str
                Location of the background FASTQ reads files [OPTIONAL]
        metadata : dict
            Input file meta data associated with their roles

            genome : str
            amb : str
            ann : str
            bwt : str
            pac : str
            sa : str
            loc : str
            bg_loc : str
                [OPTIONAL]
        output_files : dict
            Output file locations

            bam [, "bam_bg"] : str
            filtered [, "filtered_bg"] : str
            narrow_peak : str
            summits : str
            broad_peak : str
            gapped_peak : str

        Returns
        -------
        output_files : dict
            Output file locations associated with their roles, for the output

            bam [, "bam_bg"] : str
                Aligned FASTQ short read file [ and aligned background file]
                locations
            filtered [, "filtered_bg"] : str
                Filtered versions of the respective bam files
            narrow_peak : str
                Results files in bed4+1 format
            summits : str
                Results files in bed6+4 format
            broad_peak : str
                Results files in bed6+3 format
            gapped_peak : str
                Results files in bed12+3 format
        output_metadata : dict
            Output metadata for the associated files in output_files

            bam [, "bam_bg"] : Metadata
            filtered [, "filtered_bg"] : Metadata
            narrow_peak : Metadata
            summits : Metadata
            broad_peak : Metadata
            gapped_peak : Metadata
        """
        print("INPUT RUN METADATA:", metadata)

        bwa = bwaAlignerTool(self.configuration)

        print(
            "PROCESS CHIPSEQ - FILES PASSED TO TOOLS:",
            remap(input_files, "genome", "loc", "amb", "ann", "bwt", "pac", "sa")
        )
        print("PROCESS CHIPSEQ - DEFINED OUTPUT:", output_files["bam"])
        bwa_files, bwa_meta = bwa.run(
            # ideally parameter "roles" don't change
            remap(input_files,
                  "genome", "loc", "amb", "ann", "bwt", "pac", "sa"),
            remap(metadata,
                  "genome", "loc", "amb", "ann", "bwt", "pac", "sa"),
            {"output": output_files["bam"]}
        )

        if "bg_loc" in input_files:
            bwa_bg_files, bwa_bg_meta = bwa.run(
                # Small changes can be handled easily using "remap"
                remap(input_files,
                      "genome", "amb", "ann", "bwt", "pac", "sa", bg_loc="loc"),
                remap(metadata,
                      "genome", "amb", "ann", "bwt", "pac", "sa", bg_loc="loc"),
                # Intermediate outputs should be created via tempfile?
                {"output": output_files["bam_bg"]}
            )

        # For multiple files there will need to be merging into a single bam file

        # Filter the bams
        b3f = biobambam(self.configuration)
        b3f_files, b3f_meta = b3f.run(
            # Alternatively, we can rebuild the dict
            {"input": bwa_files["output"]['bam']},
            {"input": bwa_meta["output"]['bam']},
            # Intermediate outputs should be created via tempfile?
            {"output": output_files["filtered"]}
        )

        if "bg_loc" in input_files:
            b3f_bg_files, b3f_bg_meta = b3f.run(
                {"input": bwa_bg_files["output"]['bam']},
                {"input": bwa_bg_meta["output"]['bam']},
                # Intermediate outputs should be created via tempfile?
                {"output": output_files["filtered_bg"]}
            )

        ## MACS2 to call peaks
        macs_caller = macs2(self.configuration)
        macs_inputs = {"input": b3f_files["output"]['bam']}
        macs_metadt = {"input": b3f_meta["output"]['bam']}

        if "bg_loc" in input_files:
            # The dicts can be built incrementally
            macs_inputs["background"] = b3f_bg_files["output"]['bam']
            macs_metadt["background"] = b3f_bg_meta["output"]['bam']

        m_results_files, m_results_meta = macs_caller.run(
            macs_inputs, macs_metadt,
            # Outputs of the final step may match workflow outputs;
            # Extra entries in output_files will be disregarded.
            remap(output_files, 'narrow_peak', 'summits', 'broad_peak', 'gapped_peak'))

        # Outputs are collected with some name changes
        m_results_files["bam"] = bwa_files["output"]["bam"]
        m_results_files["filtered"] = b3f_files["output"]["bam"]

        # Equivalent meta data is collected
        m_results_meta["bam"] = bwa_meta["output"]["bam"]
        m_results_meta["filtered"] = b3f_meta["output"]["bam"]

        if "bg_loc" in input_files:
            m_results_files["bam_bg"] = bwa_bg_files["output"]["bam"]
            m_results_files["filtered_bg"] = b3f_bg_files["output"]["bam"]

            m_results_meta["bam_bg"] = bwa_bg_meta["output"]["bam"]
            m_results_meta["filtered_bg"] = b3f_bg_meta["output"]["bam"]

        print("CHIPSEQ RESULTS:", m_results_meta)
        return m_results_files, m_results_meta

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
    root_path = os.path.dirname(os.path.abspath(__file__))
    result = app.launch(process_chipseq,
                        root_path,
                        "tests/json/config_chipseq.json",
                        "tests/json/input_chipseq_metadata.json")

    # 2. The App has finished
    print("2. Execution finished; see " + root_path + "/results.json")
    print(result)

    return result

def prepare_files(
        dm_handler, taxon_id, genome_fa, assembly, file_loc, file_bg_loc=None):
    """
    Function to load the DM API with the required files and prepare the
    parameters passed from teh command line ready for use in the main function
    """
    # print(dm_handler.get_files_by_user("test"))

    root_name = file_loc.split("/")
    parent_dir = '/'.join(root_name[0:-1])

    genome_file = dm_handler.set_file(
        "test", genome_fa, "file", "fasta", 64000, parent_dir,
        "Assembly", taxon_id, None, None,
        meta_data={"assembly" : assembly, 'tool': 'bwa_indexer'})
    amb_file = dm_handler.set_file(
        "test", genome_fa + ".amb", "file", "amb", 64000, parent_dir,
        "Index", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly, 'tool': 'bwa_indexer'})
    ann_file = dm_handler.set_file(
        "test", genome_fa + ".ann", "file", "ann", 64000, parent_dir,
        "Index", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly, 'tool': 'bwa_indexer'})
    bwt_file = dm_handler.set_file(
        "test", genome_fa + ".bwt", "file", "bwt", 64000, parent_dir,
        "Index", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly, 'tool': 'bwa_indexer'})
    pac_file = dm_handler.set_file(
        "test", genome_fa + ".pac", "file", "pac", 64000, parent_dir,
        "Index", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly, 'tool': 'bwa_indexer'})
    sa_file = dm_handler.set_file(
        "test", genome_fa + ".sa", "file", "sa", 64000, parent_dir,
        "Index", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly, 'tool': 'bwa_indexer'})

    # Maybe it is necessary to prepare a metadata parser from json file
    # when building the Metadata objects.
    metadata = {
        "genome": Metadata(
            "Assembly", "fasta", genome_fa, None,
            {'assembly' : assembly}, genome_file),
        "amb": Metadata(
            "Index", "amb", genome_fa + ".amb", [genome_file],
            {'assembly': assembly, "tool": "bwa_indexer"}, amb_file),
        "ann": Metadata(
            "Index", "ann", genome_fa + ".ann", [genome_file],
            {'assembly': assembly, "tool": "bwa_indexer"}, ann_file),
        "bwt": Metadata(
            "Index", "bwt", genome_fa + ".bwt", [genome_file],
            {'assembly': assembly, "tool": "bwa_indexer"}, bwt_file),
        "pac": Metadata(
            "Index", "pac", genome_fa + ".pac", [genome_file],
            {'assembly': assembly, "tool": "bwa_indexer"}, pac_file),
        "sa": Metadata(
            "Index", "sa", genome_fa + ".sa", [genome_file],
            {'assembly': assembly, "tool": "bwa_indexer"}, sa_file),
    }

    fq1_file = dm_handler.set_file(
        "test", file_loc, "file", "fastq", 64000, parent_dir, "ChIP-seq",
        taxon_id, None, None, meta_data={'assembly' : assembly})
    metadata["loc"] = Metadata(
        "data_chip_seq", "fastq", file_loc, None,
        {'assembly' : assembly}, fq1_file
    )

    files = {
        "genome": genome_fa,
        "amb": genome_fa + ".amb",
        "ann": genome_fa + ".ann",
        "bwt": genome_fa + ".bwt",
        "pac": genome_fa + ".pac",
        "sa": genome_fa + ".sa",
        "loc": file_loc
    }

    root_name = file_loc.split("/")
    root_name[-1] = root_name[-1].replace('.fastq', '')

    files_out = {
        "bam": file_loc.replace(".fastq", ".bam"),
        "filtered": file_loc.replace(".fastq", "_filtered.bam"),
        "output": file_loc.replace(".fastq", ".tsv"),
        'narrow_peak': '/'.join(root_name) + '_filtered_peaks.narrowPeak',
        'summits': '/'.join(root_name) + '_filtered_summits.bed',
        'broad_peak': '/'.join(root_name) + '_filtered_peaks.broadPeak',
        'gapped_peak': '/'.join(root_name) + '_filtered_peaks.gappedPeak'
    }

    if file_bg_loc:
        fq2_file = dm_handler.set_file(
            "test", file_bg_loc, "fastq", "file", "ChIP-seq", taxon_id, None, None,
            meta_data={'assembly' : assembly})

        metadata["bg_loc"] = Metadata(
            "data_chip_seq", "fastq", file_bg_loc, None,
            {'assembly': assembly, 'main_expt': fq1_file, 'background': True},
            fq2_file
        )

        files["bg_loc"] = file_bg_loc

        files_out["bam_bg"] = file_bg_loc.replace(".fastq", ".bam"),
        files_out["filtered_bg"] = file_bg_loc.replace(".fastq", "_filtered.bam"),


    # print(dm_handler.get_files_by_user("test"))

    return [files, files_out, metadata]

# ------------------------------------------------------------------------------

def remap(indict, *args, **kwargs):
    """
    Re-map keys of indict using information from arguments.

    Non-keyword arguments are keys of input dictionary that are passed
    unchanged to the output. Keyword arguments must be in the form

    old="new"

    and act as a translation table for new key names.
    """
    outdict = {role: indict[role] for role in args}
    outdict.update(
        {new: indict[old] for old, new in kwargs.items()}
    )
    return outdict

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
    PARSER.add_argument("--json", help="Location of FASTQ background file", default=None)

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    TAXON_ID = ARGS.taxon_id
    GENOME_FA = ARGS.genome
    ASSEMBLY = ARGS.assembly
    FILE_LOC = ARGS.file
    FILE_BG_LOC = ARGS.bgd_file
    JSON_CONFIG = ARGS.json

    if JSON_CONFIG is not None:
        RESULTS = main_json()
    else:
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
        print(DM_HANDLER.get_files_by_user("test"))

    print(RESULTS)
