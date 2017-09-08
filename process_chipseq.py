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
            Input file locations associated with their roles
        metadata : dict
            Input file meta data associated with their roles
        output_files : dict
             Output file locations associated with their roles

        Returns
        -------
        output_files : dict
             Output file locations associated with their roles, for the output
             bam, bed and tsv files
        """

        # XXX input_files.keys()  # also metadata.keys()
        # "genome", "amb", "ann", "bwt", "pac", "sa", "loc" [, "bg_loc"]

        # XXX output_files.keys()
        # "bam", "filtered", "output"

        bwa = bwaAlignerTool(self.configuration)
        bwa_files, bwa_meta = bwa.run(
            # ideally parameter "roles" don't change
            remap(input_files,
                  "genome", "loc", "amb", "ann", "bwt", "pac"),
            remap(metadata,
                  "genome", "loc", "amb", "ann", "bwt", "pac"),
            {"output": output_files["bam"]}
        )

        if "bg_loc" in input_files:
            bwa_bg_files, bwa_bg_meta = bwa.run(
                # XXX: small changes can be handled easily using "remap"
                remap(input_files,
                      "genome", "amb", "ann", "bwt", "pac", bg_loc="loc"),
                remap(metadata,
                      "genome", "amb", "ann", "bwt", "pac", bg_loc="loc"),
                # XXX intermediate outputs should be created via tempfile?
                {"output": "bg.bam"}
            )

        # For multiple files there will need to be merging into a single bam file

        # Filter the bams
        b3f = biobambam(self.configuration)
        b3f_files, b3f_meta = b3f.run(
            # XXX alternatively, we can rebuild the dict
            {"input": bwa_files["output"]['bam']},
            {"input": bwa_meta["output"]['bam']},
            # XXX intermediate outputs should be created via tempfile?
            {"output": "filtered.bam"}
        )

        if "bg_loc" in input_files:
            b3f_bg_files, b3f_bg_meta = b3f.run(
                {"input": bwa_bg_files["output"]['bam']},
                {"input": bwa_bg_meta["output"]['bam']},
                # XXX intermediate outputs should be created via tempfile?
                {"output": "bg_filtered.bam"}
            )

        ## MACS2 to call peaks
        macs_caller = macs2(self.configuration)
        macs_inputs = {"input": b3f_files["output"]['bam']}
        macs_metadt = {"input": b3f_meta["output"]['bam']}

        if "bg_loc" in input_files:
            # XXX the dicts can be built incrementally
            macs_inputs["background"] = b3f_bg_files["output"]['bam']
            macs_metadt["background"] = b3f_bg_meta["output"]['bam']

        m_results_files, m_results_meta = macs_caller.run(
            macs_inputs, macs_metadt,
            # XXX outputs of the final step may match workflow outputs;
            # XXX extra entries in output_files will be disregarded.
            output_files)

        # XXX outputs are collected with some name changes
        m_results_files.update(
            remap(bwa_files, output="bam")
        ).update(
            remap(b3f_files, output="filtered")
        )

        # XXX or equivalently
        m_results_meta["bam"] = bwa_meta["output"]
        m_results_meta["filtered"] = b3f_meta["output"]

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
    amb_file = dm_handler.set_file(
        "test", genome_fa + ".amb", "amb", "Assembly", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly})
    ann_file = dm_handler.set_file(
        "test", genome_fa + ".ann", "ann", "Assembly", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly})
    bwt_file = dm_handler.set_file(
        "test", genome_fa + ".bwt", "bwt", "Assembly", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly})
    pac_file = dm_handler.set_file(
        "test", genome_fa + ".pac", "pac", "Assembly", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly})
    sa_file = dm_handler.set_file(
        "test", genome_fa + ".sa", "sa", "Assembly", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly})

    # Maybe it is necessary to prepare a metadata parser from json file
    # when building the Metadata objects.
    metadata = {
        "genome": Metadata("fasta", "Assembly", None, {'assembly' : assembly}, genome_file),
        "amb": Metadata("index", "Assembly", [genome_file], {'assembly' : assembly}, amb_file),
        "ann": Metadata("index", "Assembly", [genome_file], {'assembly' : assembly}, ann_file),
        "bwt": Metadata("index", "Assembly", [genome_file], {'assembly' : assembly}, bwt_file),
        "pac": Metadata("index", "Assembly", [genome_file], {'assembly' : assembly}, pac_file),
        "sa": Metadata("index", "Assembly", [genome_file], {'assembly' : assembly}, sa_file),
    }

    fq1_file = dm_handler.set_file(
        "test", file_loc, "fastq", "ChIP-seq", taxon_id, None, [],
        meta_data={'assembly' : assembly})
    metadata["loc"] = Metadata(
        "fastq", "ChIP-seq", None,
        {'assembly' : assembly},
        fq1_file
    )
    if file_bg_loc:
        fq2_file = dm_handler.set_file(
            "test", file_bg_loc, "fastq", "ChIP-seq", taxon_id, None, [],
            meta_data={'assembly' : assembly})
        metadata["bg_loc"] = Metadata(
            "fastq", "ChIP-seq", None,
            {'assembly' : assembly, 'background' : True},
            fq2_file
        )

    print(dm_handler.get_files_by_user("test"))

    files = {
        "genome": genome_fa,
        "amb": genome_fa + ".amb",
        "ann": genome_fa + ".ann",
        "bwt": genome_fa + ".bwt",
        "pac": genome_fa + ".pac",
        "sa": genome_fa + ".sa",
        "loc": file_loc
    }
    if file_bg_loc:
        files["bg_loc"] = file_bg_loc

    files_out = {
        "bam": genome_fa + "_out.bam",
        "filtered": genome_fa + "_filtered.bam",
        "output": genome_fa + ".tsv"
    }

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
