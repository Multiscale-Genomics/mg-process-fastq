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
import os
import sys

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT
    from dummy_pycompss import task
    from dummy_pycompss import compss_wait_on

#from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

from tool.common import common

# ------------------------------------------------------------------------------

class bwaAlignerTool(Tool):
    """
    Tool for aligning sequence reads to a genome using BWA
    """

    def __init__(self, configuration=None):
        """
        Init function
        """
        print("BWA Aligner")
        Tool.__init__(self)

    @task(returns=int, genome_file_loc=FILE_IN, read_file_loc=FILE_IN,
          bam_loc=FILE_OUT, amb_loc=FILE_IN, ann_loc=FILE_IN, bwt_loc=FILE_IN,
          pac_loc=FILE_IN, sa_loc=FILE_IN, isModifier=False)
    def bwa_aligner(self, genome_file_loc, read_file_loc, bam_loc, amb_loc, # pylint: disable=unused-argument,too-many-arguments
                    ann_loc, bwt_loc, pac_loc, sa_loc): # pylint: disable=unused-argument
        """
        BWA Aligner

        Parameters
        ----------
        genome_file_loc : str
            Location of the genomic fasta
        read_file_loc : str
            Location of the FASTQ file

        Returns
        -------
        bam_loc : str
            Location of the output file
        """

        # Rename all files for the genome and indexes so that they are in the
        # expected format by BWA - Might not be an issue with PyCOMPSs v2.1
        # Care needs to be taken when running this locally and the renaming
        # should not happen
        os.rename(amb_loc, genome_file_loc + ".amb")
        os.rename(ann_loc, genome_file_loc + ".ann")
        os.rename(bwt_loc, genome_file_loc + ".bwt")
        os.rename(pac_loc, genome_file_loc + ".pac")
        os.rename(sa_loc, genome_file_loc + ".sa")

        od_list = genome_file_loc.split("/")
        #od_list = narrowpeak.split("/")
        output_dir = "/".join(od_list[0:-1])

        print("Input files:", os.listdir(output_dir))

        common_handle = common()
        bam_loc1 = common_handle.bwa_align_reads(genome_file_loc, read_file_loc, bam_loc)

        print("Output files:", os.listdir(output_dir))
        print("BWA ALIGNER TASK:", os.path.isfile(bam_loc), os.path.getsize(bam_loc))
        print("Files to return:", bam_loc, bam_loc1)

        return 0

    def run(self, input_files, output_files, metadata=None):
        """
        The main function to align bam files to a genome using BWA

        Parameters
        ----------
        input_files : list
            File 0 is the genome file location, file 1 is the FASTQ file

        Returns
        -------
        output : list
            First element is a list of output_bam_files, second element is the
            matching meta data
        """

        output_metadata = {}
        out_bam = input_files[1].replace(".fastq", '.bam')

        results = self.bwa_aligner(
            input_files[0], input_files[1], out_bam, input_files[2],
            input_files[3], input_files[4], input_files[5], input_files[6])

        results = compss_wait_on(results)

        print("BWA ALIGNER:", os.path.isfile(out_bam))

        return ([out_bam], [output_metadata])

# ------------------------------------------------------------------------------
