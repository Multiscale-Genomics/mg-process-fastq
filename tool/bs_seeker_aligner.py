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

import shlex
import subprocess
import sys
import os.path
import tarfile

import pysam

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_INOUT, FILE_OUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on, compss_open, barrier
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_INOUT, FILE_OUT, IN
    from dummy_pycompss import task
    from dummy_pycompss import compss_wait_on, compss_open, barrier

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

from tool.fastq_splitter import fastq_splitter

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

    @task(bam_file=FILE_INOUT)
    def bam_sort(self, bam_file):
        """
        Wrapper for the pysam SAMtools sort function

        Parameters
        ----------
        bam_file : str
            Location of the bam file to sort
        """
        pysam.sort("-o", bam_file, "-T", bam_file + "_sort", bam_file)
        return True

    @task(bam_file_1=FILE_INOUT, bam_file_2=FILE_IN)
    def bam_merge(self, bam_file_1, bam_file_2):
        """
        Wrapper for the pysam SAMtools merge function

        Parameters
        ----------
        bam_file_1 : str
            Location of the bam file to merge into
        bam_file_2 : str
            Location of the bam file that is to get merged into bam_file_1
        """
        pysam.merge(bam_file_1 + "_merge.bam", bam_file_1, bam_file_2)

        with open(bam_file_1 + "_merge.bam", "rb") as f_in:
            with open(bam_file_1, "wb") as f_out:
                f_out.write(f_in.read())

        return True

    @task(bam_in=FILE_IN, bam_out=FILE_OUT)
    def bam_copy(self, bam_in, bam_out):
        """
        Wrapper function to copy from one bam file to another

        Parameters
        ----------
        bam_in : str
            Location of the input bam file
        bam_out : str
            Location of the output bam file
        """
        with open(bam_in, "rb") as f_in:
            with open(bam_out, "wb") as f_out:
                f_out.write(f_in.read())

        return True

    @task(bam_file=FILE_IN, bam_idx_file=FILE_OUT)
    def bam_index(self, bam_file, bam_idx_file):
        """
        Wrapper for the pysam SAMtools merge function

        Parameters
        ----------
        bam_file : str
            Location of the bam file that is to be indexed
        bam_idx_file : str
            Location of the bam index file (.bai)
        """
        pysam.index(bam_file, bam_file + "_tmp.bai")

        with open(bam_file + "_tmp.bai", "rb") as f_in:
            with open(bam_idx_file, "wb") as f_out:
                f_out.write(f_in.read())

        return True

    @task(
        returns=bool, isModifier=False,
        input_fastq_1=FILE_IN, input_fastq_2=FILE_IN,
        aligner=IN, aligner_path=IN, bss_path=IN,
        genome_fasta=FILE_IN, genome_idx=FILE_IN, bam_out=FILE_OUT
    )
    def bs_seeker_aligner(
            self, input_fastq_1, input_fastq_2, aligner, aligner_path, bss_path,
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
        g_dir = genome_idx.split("/")
        g_dir = "/".join(g_dir[:-1])

        tar = tarfile.open(genome_idx)
        tar.extractall(path=g_dir)
        tar.close()

        command_line = (
            "python " + bss_path + "/bs_seeker2-align.py"
            " --input_1 " + input_fastq_1 + ""
            " --input_2 " + input_fastq_2 + ""
            " --aligner " + aligner + " --path " + aligner_path + ""
            " --genome " + genome_fasta + " -d " + g_dir + ""
            " --bt2-p 4 -o " + bam_out + "_tmp.bam"
        ).format()
        print("command for aligner : ", command_line)
        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        pysam.sort("-o", bam_out + "_tmp.bam", "-T", bam_out + "_tmp.bam" + "_sort", bam_out + "_tmp.bam")

        with open(bam_out + "_tmp.bam", "rb") as f_in:
            with open(bam_out, "wb") as f_out:
                f_out.write(f_in.read())

        return True

    def run(self, input_files, metadata, output_files):
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

        genome_fasta = input_files["genome"]
        genome_idx = input_files["index"]

        fqs = fastq_splitter()

        fastq1 = input_files["fastq1"]
        fastq_file_gz = fastq1 + ".tar.gz"
        if "fastq2" in input_files:
            fastq2 = input_files["fastq2"]
            fastq_file_list = fqs.paired_splitter(
                fastq1, fastq2, fastq1 + ".tar.gz"
            )
        else:
            fastq_file_list = fqs.single_splitter(
                fastq1, fastq1 + ".tar.gz"
            )

        fastq_file_list = compss_wait_on(fastq_file_list)

        with compss_open(fastq_file_gz, "rb") as f_in:
            with open(fastq_file_gz, "wb") as f_out:
                f_out.write(f_in.read())

        gz_data_path = fastq_file_gz.split("/")
        gz_data_path = "/".join(gz_data_path[:-1])

        tar = tarfile.open(fastq_file_gz)
        tar.extractall(path=gz_data_path)
        tar.close()

        aligner = metadata['aligner']
        aligner_path = metadata['aligner_path']
        bss_path = metadata['bss_path']

        # input and output share most metadata
        output_metadata = {}

        output_bam_file = output_files["bam"]
        output_bai_file = output_files["bai"]

        output_bam_list = []
        for fastq_file_pair in fastq_file_list:
            tmp_fq1 = gz_data_path + "/tmp/" + fastq_file_pair[0]
            tmp_fq2 = gz_data_path + "/tmp/" + fastq_file_pair[1]
            print("WGBS - gz_data_path:", gz_data_path)
            print("WGBS - fastq_file_pair - 0:", tmp_fq1,
                  os.path.isfile(tmp_fq1), os.path.getsize(tmp_fq1))
            print("WGBS - fastq_file_pair - 1:", tmp_fq2,
                  os.path.isfile(tmp_fq2), os.path.getsize(tmp_fq2))

            output_bam_file_tmp = tmp_fq1 + ".bam"
            output_bam_list.append(output_bam_file_tmp)

            print(
                "FILES:", tmp_fq1, tmp_fq2,
                aligner, aligner_path, bss_path,
                genome_fasta, genome_idx,
                output_bam_file_tmp)

            results = self.bs_seeker_aligner(
                tmp_fq1, tmp_fq2,
                aligner, aligner_path, bss_path,
                genome_fasta, genome_idx,
                output_bam_file_tmp
            )

        barrier()

        results = self.bam_copy(output_bam_list.pop(0), output_bam_file)
        results = compss_wait_on(results)

        while True:
            if len(output_bam_list) == 0:
                break
            results = self.bam_merge(output_bam_file, output_bam_list.pop(0))
            results = compss_wait_on(results)

        results = self.bam_sort(output_bam_file)
        results = compss_wait_on(results)

        results = self.bam_index(output_bam_file, output_bai_file)
        results = compss_wait_on(results)

        output_metadata = {
            "bam": Metadata(
                "wgbs", "bam", [metadata["genome"].file_path],
                {
                    "assembly": metadata["genome"].meta_data["assembly"],
                    "tool": "bs_seeker_aligner"
                }
            ),
            "bai": Metadata(
                "wgbs", "bai", [metadata["genome"].file_path],
                {
                    "assembly": metadata["genome"].meta_data["assembly"],
                    "tool": "bs_seeker_aligner"
                }
            )
        }

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
