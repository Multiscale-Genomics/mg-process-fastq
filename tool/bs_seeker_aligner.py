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

import shlex
import subprocess
import sys
import tarfile

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT, IN
    from dummy_pycompss import task

from basic_modules.tool import Tool

# ------------------------------------------------------------------------------

class bssAlignerTool(Tool):
    """
    Script from BS-Seeker2 for building the index for alignment. In this case
    it uses Bowtie2.
    """

    def __init__(self):
        """
        Init function
        """
        print("BS-Seeker Aligner")
        Tool.__init__(self)

    @task(
        input_fastq1=FILE_IN, input_fastq2=FILE_IN,
        aligner=IN, aligner_path=IN, bss_path=IN,
        genome_fasta=FILE_IN, genome_idx=FILE_IN, bam_out=FILE_OUT)
    def bs_seeker_aligner(
            self, input_fastq1, input_fastq2, aligner, aligner_path, bss_path,
            genome_fasta, genome_idx, bam_out):
        """
        Alignment of the paired ends to the reference genome

        Generates bam files for the alignments

        This is performed by running the external program rather than
        reimplementing the code from the main function to make it easier when
        it comes to updating the changes in BS-Seeker2

        Parameters
        ----------
        input_fastq1 : str
            Location of paired end FASTQ file 1
        input_fastq2 : str
            Location of paired end FASTQ file 2
        aligner : str
            Aligner to use
        aligner_path : str
            Location of the aligner
        genome_fasta : str
            Location of the genome FASTA file
        genome_idx : str
            Location of the tar.gz genome index file
        bam_out : str
            Location of the aligned bam file

        Returns
        -------
        bam_out : file
            Location of the BAM file generated during the alignment.
        """
        g_dir = genome_fasta.split("/")
        g_dir = "/".join(g_dir[0:-1])

        tar = tarfile.open(genome_idx)
        tar.extractall()
        tar.close()

        command_line = (
            "python " + bss_path + "/bs_seeker2-align.py"
            " --input_1 " + input_fastq1 + " --input_2 " + input_fastq2 + ""
            " --aligner " + aligner + " --path " + aligner_path + ""
            " --genome " + genome_fasta + " -d " + g_dir + ""
            " --bt2-p 4 -o " + bam_out
        ).format()
        print ("command for aligner : ", command_line)
        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        return True


    def run(self, input_files, output_files, metadata=None):
        """
        Tool for indexing the genome assembly using BS-Seeker2. In this case it
        is using Bowtie2

        Parameters
        ----------
        input_files : list
            FASTQ file
        output_files : list
            Results files.
        metadata : list

        Returns
        -------
        array : list
            Location of the filtered FASTQ file
        """

        genome_fasta = input_files[0]
        genome_idx = input_files[1]
        fastq_file_1 = input_files[2]
        fastq_file_2 = input_files[3]

        aligner = metadata['aligner']
        aligner_path = metadata['aligner_path']
        bss_path = metadata['bss_path']

        # input and output share most metadata
        output_metadata = {}

        output_file = fastq_file_1.replace('_1.fastq', '.bam')

        results = self.bs_seeker_aligner(
            fastq_file_1, fastq_file_2,
            aligner, aligner_path, bss_path,
            genome_fasta, genome_idx,
            output_file
        )

        if results is False:
            pass

        return ([output_file], output_metadata)

# ------------------------------------------------------------------------------
