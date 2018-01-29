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
import tarfile

import pysam

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on, compss_open, barrier
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import compss_wait_on, compss_open, barrier # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

from tool.fastq_splitter import fastq_splitter
from tool.bam_utils import bamUtilsTask

# ------------------------------------------------------------------------------

class bssAlignerTool(Tool):
    """
    Script from BS-Seeker2 for building the index for alignment. In this case
    it uses Bowtie2.
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
        print("BS-Seeker Aligner")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

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
        script = bss_path + "/bs_seeker2-align.py"
        params = [
            "--input_1", input_fastq_1,
            "--input_2", input_fastq_2,
            "--aligner", aligner,
            "--path", aligner_path,
            "--genome", genome_fasta,
            "--bt2-p", "4",
            "-o", bam_out + "_tmp.bam",
            "-f", "bam"
        ]

        results = self.run_aligner(genome_idx, bam_out, script, params)

        return results

    @task(
        returns=bool, isModifier=False,
        input_fastq=FILE_IN, aligner=IN, aligner_path=IN, bss_path=IN,
        genome_fasta=FILE_IN, genome_idx=FILE_IN, bam_out=FILE_OUT
    )
    def bs_seeker_aligner_single(
            self, input_fastq, aligner, aligner_path, bss_path,
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
        script = bss_path + "/bs_seeker2-align.py"
        params = [
            "-i", input_fastq,
            "--aligner", aligner,
            "--path", aligner_path,
            "--genome", genome_fasta,
            "--bt2-p", "4",
            "-o", bam_out + "_tmp.bam",
            "-f", "bam"
        ]

        results = self.run_aligner(genome_idx, bam_out, script, params)

        return results

    @staticmethod
    def run_aligner(genome_idx, bam_out, script, params):
        """
        Run the aligner

        Parameters
        ----------
        genome_idx : str
            Location of the genome index archive
        bam_out : str
            Location of the output bam file
        script : str
            Location of the BS Seeker2 aligner script
        params : list
            Parameter list for the aligner

        Returns
        -------
        bool
            True if the function completed successfully
        """
        g_dir = genome_idx.split("/")
        g_dir = "/".join(g_dir[:-1])
        params += ["-d", g_dir]

        try:
            tar = tarfile.open(genome_idx)
            tar.extractall(path=g_dir)
            tar.close()
        except IOError:
            logger.fatal("WGBS - BS SEEKER2: Missing index archive")
            return False

        command_line = (
            "python " + script + " " + " ".join(params)
        ).format()
        print("command for aligner : ", command_line)
        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        pysam.sort("-o", bam_out + "_tmp.bam",
                   "-T", bam_out + "_tmp.bam" + "_sort", bam_out + "_tmp.bam")

        try:
            with open(bam_out + "_tmp.bam", "rb") as f_in:
                with open(bam_out, "wb") as f_out:
                    f_out.write(f_in.read())
        except IOError:
            logger.fatal("WGBS - BS SEEKER2: Failed sorting")
            return False

        return True

    def run(self, input_files, input_metadata, output_files):
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

        try:
            if "bss_path" in self.configuration:
                bss_path = self.configuration["bss_path"]
            else:
                raise KeyError
            if "aligner_path" in self.configuration:
                aligner_path = self.configuration["aligner_path"]
            else:
                raise KeyError
            if "aligner" in self.configuration:
                aligner = self.configuration["aligner"]
            else:
                raise KeyError
        except KeyError:
            logger.fatal("WGBS - BS SEEKER2: Unassigned configuration variables")

        genome_fasta = input_files["genome"]
        genome_idx = input_files["index"]

        sources = [input_files["genome"]]

        fqs = fastq_splitter()

        fastq1 = input_files["fastq1"]
        sources.append(input_files["fastq1"])

        fastq_file_gz = fastq1 + ".tar.gz"
        if "fastq2" in input_files:
            fastq2 = input_files["fastq2"]
            sources.append(input_files["fastq2"])
            fastq_file_list = fqs.paired_splitter(
                fastq1, fastq2, fastq_file_gz
            )
        else:
            fastq_file_list = fqs.single_splitter(
                fastq1, fastq_file_gz
            )

        fastq_file_list = compss_wait_on(fastq_file_list)
        if not fastq_file_list:
            logger.fatal("FASTQ SPLITTER: run failed")
            return {}, {}

        if hasattr(sys, '_run_from_cmdl') is True:
            pass
        else:
            with compss_open(fastq_file_gz, "rb") as f_in:
                with open(fastq_file_gz, "wb") as f_out:
                    f_out.write(f_in.read())

        gz_data_path = fastq_file_gz.split("/")
        gz_data_path = "/".join(gz_data_path[:-1])

        try:
            tar = tarfile.open(fastq_file_gz)
            tar.extractall(path=gz_data_path)
            tar.close()
        except tarfile.TarError:
            logger.fatal("Split FASTQ files: Malformed tar file")
            return {}, {}

        # input and output share most metadata
        output_metadata = {}

        output_bam_file = output_files["bam"]
        output_bai_file = output_files["bai"]

        output_bam_list = []
        for fastq_file_pair in fastq_file_list:
            if "fastq2" in input_files:
                tmp_fq1 = gz_data_path + "/tmp/" + fastq_file_pair[0]
                tmp_fq2 = gz_data_path + "/tmp/" + fastq_file_pair[1]
                output_bam_file_tmp = tmp_fq1 + ".bam"
                output_bam_list.append(output_bam_file_tmp)

                results = self.bs_seeker_aligner(
                    tmp_fq1, tmp_fq2,
                    aligner, aligner_path, bss_path,
                    genome_fasta, genome_idx,
                    output_bam_file_tmp
                )
            else:
                tmp_fq = gz_data_path + "/tmp/" + fastq_file_pair[0]
                output_bam_file_tmp = tmp_fq + ".bam"
                output_bam_list.append(output_bam_file_tmp)

                results = self.bs_seeker_aligner_single(
                    tmp_fq,
                    aligner, aligner_path, bss_path,
                    genome_fasta, genome_idx,
                    output_bam_file_tmp
                )

        barrier()

        bam_handler = bamUtilsTask()

        results = bam_handler.bam_copy(output_bam_list.pop(0), output_bam_file)
        results = compss_wait_on(results)

        if results is False:
            logger.fatal("BS SEEKER2 Aligner: Bam copy failed")
            return {}, {}

        while True:
            if len(output_bam_list) == 0:
                break
            results = bam_handler.bam_merge(output_bam_file, output_bam_list.pop(0))
            results = compss_wait_on(results)

            if results is False:
                logger.fatal("BS SEEKER2 Aligner: Bam merging failed")
                return {}, {}

        results = bam_handler.bam_sort(output_bam_file)
        results = compss_wait_on(results)

        if results is False:
            logger.fatal("BS SEEKER2 Aligner: Bam sorting failed")
            return {}, {}

        results = bam_handler.bam_index(output_bam_file, output_bai_file)
        results = compss_wait_on(results)

        if results is False:
            logger.fatal("BS SEEKER2 Aligner: Bam indexing failed")
            return {}, {}

        output_metadata = {
            "bam": Metadata(
                data_type="data_wgbs",
                file_type="BAM",
                file_path=output_bam_file,
                sources=sources,
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "bwa_indexer"
                }
            ),
            "bai": Metadata(
                data_type="data_wgbs",
                file_type="BAI",
                file_path=output_bai_file,
                sources=[input_metadata["genome"].file_path],
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "bwa_indexer"
                }
            )
        }

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
