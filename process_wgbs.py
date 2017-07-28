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
import pysam

from basic_modules.workflow import Workflow
from basic_modules.metadata import Metadata

from dmp import dmp

from tool.bs_seeker_aligner import bssAlignerTool
from tool.bs_seeker_filter import filterReadsTool
from tool.bs_seeker_indexer import bssIndexerTool
from tool.bs_seeker_methylation_caller import bssMethylationCallerTool

from fastqreader import fastqreader

# ------------------------------------------------------------------------------

class process_wgbs(Workflow):
    """
    Functions for downloading and processing whole genome bisulfate sequencings
    (WGBS) files. Files are filtered, aligned and analysed for points of
    methylation
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
        if configuration is None : 
            configuration = {}
        self.configuration.update(configuration)

    def single_splitter(self, in_file1, tag='tmp'):
        """
        Function to divide the FastQ files into separte sub files of 1000000
        sequences so that the aligner can run in parallel.

        Parameters
        ----------
        in_file1 : str
            Location of first FASTQ file
        tag : str
            DEFAULT = tmp
            Tag used to identify the files. Useful if this is getting run
            manually on a single machine multiple times to prevent collisions of
            file names


        Returns
        -------
        Returns: Returns a list of the files that have been generated.
                 Each sub list containing the two paired end files for that
                 subset.
        paired_files : list
            List of lists of pair end files. Each sub list containing the two
            paired end files for that subset.
        """

        fqr = fastqreader()
        fqr.openFastQ(in_file1)
        fqr.createOutputFiles(tag)

        record1 = fqr.next(1)

        count_r3 = 0

        file_loc_1 = fqr.fastq1.split("/")
        file_loc_1[-1] = file_loc_1[-1].replace(
            ".fastq",
            "." + str(fqr.output_tag) + "_" + str(fqr.output_file_count) + ".fastq")
        file_loc_1.insert(-1, tag)

        files_out = [["/".join(file_loc_1)]]

        while fqr.eof(1) is False:
            fqr.writeOutput(record1, 1)
            record1 = fqr.next(1)
            count_r3 += 1

            if count_r3 % 1000000 == 0:
                fqr.incrementOutputFiles()
                file_loc_1 = fqr.fastq1.split("/")
                file_loc_1[-1] = file_loc_1[-1].replace(
                    ".fastq",
                    "." + str(fqr.output_tag) + "_" + str(fqr.output_file_count) + ".fastq")
                file_loc_1.insert(-1, "tmp")

                files_out.append(["/".join(file_loc_1)])

        fqr.closeFastQ()
        fqr.closeOutputFiles()

        return files_out

    def paired_splitter(self, in_file1, in_file2, tag='tmp'):
        """
        Function to divide the FastQ files into separte sub files of 1000000
        sequences so that the aligner can run in parallel.

        Parameters
        ----------
        in_file1 : str
            Location of first paired end FASTQ file
        in_file2 : str
            Location of second paired end FASTQ file
        tag : str
            DEFAULT = tmp
            Tag used to identify the files. Useful if this is getting run
            manually on a single machine multiple times to prevent collisions of
            file names


        Returns
        -------
        Returns: Returns a list of lists of the files that have been generated.
                 Each sub list containing the two paired end files for that
                 subset.
        paired_files : list
            List of lists of pair end files. Each sub list containing the two
            paired end files for that subset.
        """

        fqr = fastqreader()
        fqr.openFastQ(in_file1, in_file2)
        fqr.createOutputFiles(tag)

        record1 = fqr.next(1)
        record2 = fqr.next(2)

        count_r1 = 0
        count_r2 = 0
        count_r3 = 0

        file_loc_1 = fqr.fastq1.split("/")
        file_loc_1[-1] = file_loc_1[-1].replace(
            ".fastq",
            "." + str(fqr.output_tag) + "_" + str(fqr.output_file_count) + ".fastq")
        file_loc_1.insert(-1, tag)

        file_loc_2 = fqr.fastq2.split("/")
        file_loc_2[-1] = file_loc_2[-1].replace(
            ".fastq",
            "." + str(fqr.output_tag) + "_" + str(fqr.output_file_count) + ".fastq")
        file_loc_2.insert(-1, tag)
        files_out = [["/".join(file_loc_1), "/".join(file_loc_2)]]

        while fqr.eof(1) is False and fqr.eof(2) is False:
            r1_id = record1["id"].split(" ")
            r2_id = record2["id"].split(" ")

            if r1_id[0] == r2_id[0]:
                fqr.writeOutput(record1, 1)
                fqr.writeOutput(record2, 2)

                record1 = fqr.next(1)
                record2 = fqr.next(2)

                count_r1 += 1
                count_r2 += 1
                count_r3 += 1
            elif r1_id[0] < r2_id[0]:
                record1 = fqr.next(1)
                count_r1 += 1
            else:
                record2 = fqr.next(2)
                count_r2 += 1

            if count_r3 % 1000000 == 0:
                fqr.incrementOutputFiles()
                file_loc_1 = fqr.fastq1.split("/")
                file_loc_1[-1] = file_loc_1[-1].replace(
                    ".fastq",
                    "." + str(fqr.output_tag) + "_" + str(fqr.output_file_count) + ".fastq")
                file_loc_1.insert(-1, "tmp")

                file_loc_2 = fqr.fastq2.split("/")
                file_loc_2[-1] = file_loc_2[-1].replace(
                    ".fastq",
                    "." + str(fqr.output_tag) + "_" + str(fqr.output_file_count) + ".fastq")
                file_loc_2.insert(-1, "tmp")

                files_out.append(["/".join(file_loc_1), "/".join(file_loc_2)])

        fqr.closeFastQ()
        fqr.closeOutputFiles()

        return files_out


    def run(self, input_files, output_files, metadata):
        """
        This pipeline processes paired-end FASTQ files to identify
        methylated regions within the genome.

        Parameters
        ----------
        input_files : list
            List of strings for the locations of files. These should include:

            genome_fa : str
                Genome assembly in FASTA
            fastq1 : str
                FASTQ file for the first pair end file
            fastq2 : str
                FASTQ file for the second pair end file

        Returns
        -------
        wig_file : str
            Location of the wig file containing the methylation peak calls
        cgmap_file : str
            Location of the CGmap file generated by BS-Seeker2
        atcgmap_file : str
            Location of the ATCGmap file generated by BS-Seeker2
        """
        genome_fa = input_files[0]
        fastq1 = input_files[1]
        fastq2 = input_files[2]
        bt2_1 = input_files[3]
        bt2_2 = input_files[4]
        bt2_3 = input_files[5]
        bt2_4 = input_files[6]
        bt2_rev_1 = input_files[7]
        bt2_rev_2 = input_files[8]

        output_metadata = {}

        # Filter the FASTQ reads to remove duplicates
        frt = filterReadsTool()
        fastq1f, filter1_meta = frt.run([fastq1], metadata)
        fastq2f, filter2_meta = frt.run([fastq2], metadata)

        output_metadata['fastq1'] = filter1_meta
        output_metadata['fastq2'] = filter2_meta

        # Build the matching WGBS genome index
        builder = bssIndexerTool()
        genome_idx, gidx_meta = builder.run(
            [genome_fa],
            metadata
        )
        output_metadata['genome_idx'] = gidx_meta

        # Split the FASTQ files into smaller, easier to align packets
        if fastq2 is not None:
            tmp_fastq = self.paired_splitter(fastq1f[0], fastq2f[0], 'tmp')
        else:
            tmp_fastq = self.single_splitter(fastq1f[0], 'tmp')

        bam_sort_files = []
        bam_merge_files = []
        fastq_for_alignment = []
        for bams in tmp_fastq:
            bam_root = bams[0] + "_bspe.bam"
            tmp = bams
            tmp.append(bam_root)

            fastq_for_alignment.append(tmp)

            bam_sort_files.append([bam_root, bam_root + ".sorted.bam"])
            bam_merge_files.append(bam_root + ".sorted.bam")

        # Run the bs_seeker2-align.py steps on the split up fastq files
        for ffa in fastq_for_alignment:
            bss_aligner = bssAlignerTool()
            if fastq2 is not None:
                input_files = [
                    genome_fa, ffa[0], ffa[1],
                    bt2_1, bt2_2, bt2_3, bt2_4, bt2_rev_1, bt2_rev_2]
                output_files = [ffa[2]]
            else:
                input_files = [
                    genome_fa, ffa[0], None,
                    bt2_1, bt2_2, bt2_3, bt2_4, bt2_rev_1, bt2_rev_2]
                output_files = [ffa[1]]

            bam, bam_meta = bss_aligner.run(input_files, output_files, metadata)

            if 'alignment' in output_metadata:
                output_metadata['alignment'].append(bam_meta)
            else:
                output_metadata['alignment'] = [bam_meta]

        # Sort and merge the aligned bam files
        # Pre-sort the original input bam files
        for bfs in bam_sort_files:
            pysam.sort("-o", bfs[1], bfs[0])

        f_bam = fastq1.split("/")
        f_bam[-1] = f_bam[-1].replace(".fastq", ".sorted.bam")
        out_bam_file = "/".join(f_bam)

        pysam.merge(out_bam_file, *bam_merge_files)

        pysam.sort("-o", out_bam_file, "-T", out_bam_file + "_sort", out_bam_file)

        pysam.index(out_bam_file)

        # Methylation peak caller
        peak_caller_handle = bssMethylationCallerTool()

        metadata['index_path'] = genome_fa + '_bowtie2'
        peak_files, peak_meta = peak_caller_handle.run(
            [out_bam_file],
            [],
            metadata
        )
        output_metadata['peak_calling'] = peak_meta




        return (peak_files, output_metadata)


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
    result = app.launch(process_wgbs, input_files, output_files, input_metadata,
                        {})

    # 2. The App has finished
    print("2. Execution finished")
    print(result)
    return result

def prepare_files(
        dm_handler, taxon_id, genome_fa, assembly, file_1_loc, file_2_loc):
    """outes and prepare the
    parameters passed from teh command line ready for use in the main function
    """
    print(dm_handler.get_files_by_user("test"))

    genome_file = dm_handler.set_file(
        "test", genome_fa, "fasta", "Assembly", taxon_id, None, [],
        meta_data={"assembly" : assembly})
    dm_handler.set_file(
        "test", genome_fa + ".1.bt2", "bt2", "Index", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly})
    dm_handler.set_file(
        "test", genome_fa + ".2.bt2", "bt2", "Index", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly})
    dm_handler.set_file(
        "test", genome_fa + ".3.bt2", "bt2", "Index", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly})
    dm_handler.set_file(
        "test", genome_fa + ".4.bt2", "bt2", "Index", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly})
    dm_handler.set_file(
        "test", genome_fa + ".rev.1.bt2", "bt2", "Index", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly})
    dm_handler.set_file(
        "test", genome_fa + ".rev.2.bt2", "bt2", "Index", taxon_id, None, [genome_file],
        meta_data={'assembly' : assembly})


    # Maybe it is necessary to prepare a metadata parser from json file
    # when building the Metadata objects.
    metadata = {'files' : [
        Metadata("fasta", "Assembly"),
        Metadata("fastq", "WGBS"),
        Metadata("fastq", "WGBS"),
        Metadata("index", "Assembly"),
        Metadata("index", "Assembly"),
        Metadata("index", "Assembly"),
        Metadata("index", "Assembly"),
        Metadata("index", "Assembly"),
        Metadata("index", "Assembly"),
    ]}

    fastq1_id = dm_handler.set_file(
        "test", file_1_loc, "fasta", "WGBS", taxon_id, None, [],
        meta_data={'assembly' : assembly})

    fastq2_id = dm_handler.set_file(
        "test", file_2_loc, "fasta", "WGBS", taxon_id, None, [],
        meta_data={
            'assembly' : assembly,
            'paired_end' : fastq1_id
        })

    dm_handler.add_file_metadata(fastq1_id, 'paired_end', fastq2_id)

    files_input = [
        genome_fa,
        file_1_loc,
        file_2_loc,
        genome_fa + ".1.bt2",
        genome_fa + ".2.bt2",
        genome_fa + ".3.bt2",
        genome_fa + ".4.bt2",
        genome_fa + ".rev.1.bt2",
        genome_fa + ".rev.2.bt2",
    ]

    return (
        files_input,
        [],
        metadata
    )

# ------------------------------------------------------------------------------


if __name__ == "__main__":
    import sys
    sys._run_from_cmdl = True

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="Parse WGBS data")
    PARSER.add_argument("--fastq1", help="Location of first paired end FASTQ")
    PARSER.add_argument("--fastq2", help="Location of second paired end FASTQ", default=None)
    PARSER.add_argument("--taxon_id", help="Taxon_ID (10090)")
    PARSER.add_argument("--assembly", help="Assembly (GRCm38)")
    PARSER.add_argument("--genome", help="Genome assembly FASTA file")
    PARSER.add_argument("--aligner", help="Aligner to use (eg bowtie2)")
    PARSER.add_argument("--aligner_path", help="Directory for the aligner program")
    PARSER.add_argument("--bss_path", help="Directory for the BS-Seeker2 program")

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    TAXON_ID = ARGS.taxon_id
    GENOME_FA = ARGS.genome
    ASSEMBLY = ARGS.assembly
    FASTQ1 = ARGS.fastq1
    FASTQ2 = ARGS.fastq2
    ALIGNER = ARGS.aligner
    ALIGNER_PATH = ARGS.aligner_path
    BSS_PATH = ARGS.bss_path

    #
    # MuG Tool Steps
    # --------------
    #
    # 1. Create data files
    DM_HANDLER = dmp(test=True)

    #2. Register the data with the DMP
    PARAMS = prepare_files(DM_HANDLER, TAXON_ID, GENOME_FA, ASSEMBLY, FASTQ1, FASTQ2)

    METADATA = PARAMS[2]
    METADATA['user_id'] = 'test'
    METADATA['assembly'] = ASSEMBLY
    METADATA['aligner'] = ALIGNER
    METADATA['aligner_path'] = ALIGNER_PATH
    METADATA['bss_path'] = BSS_PATH

    # 3. Instantiate and launch the App
    RESULTS = main(PARAMS[0], PARAMS[1], PARAMS[2])

    print(RESULTS)
    print(DM_HANDLER.get_files_by_user("test"))
