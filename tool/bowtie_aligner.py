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
import os
import sys
import shutil
import tarfile

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import IN, FILE_IN, FILE_OUT, FILE_INOUT
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on, compss_open, barrier
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import IN, FILE_IN, FILE_OUT, FILE_INOUT # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task
    from utils.dummy_pycompss import compss_wait_on, compss_open, barrier

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

from tool.fastq_splitter import fastq_splitter
from tool.aligner_utils import alignerUtils
from tool.bam_utils import bamUtilsTask

# ------------------------------------------------------------------------------


class bowtie2AlignerTool(Tool):
    """
    Tool for aligning sequence reads to a genome using BWA
    """

    def __init__(self, configuration=None):
        """
        Init function
        """
        logger.info("Bowtie2 Aligner")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(returns=bool, genome_file_loc=FILE_IN, read_file_loc=FILE_IN,
          bam_loc=FILE_OUT, genome_idx=FILE_IN, aln_params=IN, isModifier=False)
    def bowtie2_aligner_single(  # pylint: disable=too-many-arguments
            self, genome_file_loc, read_file_loc, bam_loc, genome_idx, aln_params):  # pylint: disable=unused-argument
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
        g_dir = genome_idx.split("/")
        g_dir = "/".join(g_dir[:-1])

        try:
            tar = tarfile.open(genome_idx)
            tar.extractall(path=g_dir)
            tar.close()
        except IOError:
            return False

        gfl = genome_file_loc.split("/")
        genome_fa_ln = genome_idx.replace('.tar.gz', '/') + gfl[-1]
        shutil.copy(genome_file_loc, genome_fa_ln)

        if (
                os.path.isfile(genome_fa_ln) is False or
                os.path.getsize(genome_fa_ln) == 0):
            return False
        if (
                os.path.isfile(read_file_loc) is False or
                os.path.getsize(read_file_loc) == 0):
            return False

        out_bam = read_file_loc + '.out.bam'

        au_handle = alignerUtils()
        logger.info(
            "BOWTIE2 FINISHED: " + str(au_handle.bowtie2_align_reads(
                genome_fa_ln, read_file_loc, out_bam, aln_params))
        )

        try:
            with open(bam_loc, "wb") as f_out:
                with open(out_bam, "rb") as f_in:
                    f_out.write(f_in.read())
        except IOError:
            return False

        # shutil.rmtree(g_dir)

        return True

    @task(returns=bool, genome_file_loc=FILE_IN, read_file_loc1=FILE_IN,
          read_file_loc2=FILE_IN, bam_loc=FILE_OUT, genome_idx=FILE_IN,
          aln_params=IN, isModifier=False)
    def bowtie2_aligner_paired(  # pylint: disable=too-many-arguments
            self, genome_file_loc, read_file_loc1, read_file_loc2, bam_loc,
            genome_idx, aln_params):  # pylint: disable=unused-argument
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
        g_dir = genome_idx.split("/")
        g_dir = "/".join(g_dir[:-1])

        try:
            tar = tarfile.open(genome_idx)
            tar.extractall(path=g_dir)
            tar.close()
        except IOError:
            return False

        gfl = genome_file_loc.split("/")
        genome_fa_ln = genome_idx.replace('.tar.gz', '/') + gfl[-1]
        shutil.copy(genome_file_loc, genome_fa_ln)

        out_bam = read_file_loc1 + '.out.bam'
        au_handle = alignerUtils()
        logger.info(
            "BOWTIE2 FINISHED: " + str(au_handle.bowtie2_align_reads(
                genome_fa_ln, read_file_loc1, read_file_loc2, out_bam, aln_params))
        )

        try:
            with open(bam_loc, "wb") as f_out:
                with open(out_bam, "rb") as f_in:
                    f_out.write(f_in.read())
        except IOError:
            return False

        # shutil.rmtree(g_dir)

        return True

    def get_aln_params(self, params):
        """

        """
        command_params = []

        # Input Options
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]

        # Alignment Options
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]

        # Scoring Options
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]

        # Reporting Options
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]

        # Effort Options
        if "bowtie2_seed-extension-attempts_param" in params:
            command_params = command_params + [
                "-D", str(params["bowtie2_seed-extension-attempts_param"])]
        if "bowtie2_reseed_param" in params:
            command_params = command_params + [
                "-R", str(params["bowtie2_reseed_param"])]

        # Paired-end Options
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]
        if "bowtie2__param" in params:
            command_params = command_params + [
                "", str(params[""])]

        # SAM Options
        if "bowtie2_no-unal_param" in params:
            command_params = command_params.append("--no-unal")
        if "bowtie2_no-hd_param" in params:
            command_params = command_params.append("--no-hd")
        if "bowtie2_no-sq_param" in params:
            command_params = command_params.append("--no-dq")
        if "bowtie2_rg-id_param" in params:
            command_params = command_params + [
                "--rg-id", str(params["bowtie2_rg-id_param"])]
        if "bowtie2_rg_param" in params:
            command_params = command_params + [
                "--rg", str(params["bowtie2_rg_param"])]
        if "bowtie2_omit-sec-seq_param" in params:
            command_params = command_params.append("--omit-sec-seq")
        if "bowtie2_soft-clipped-unmapped-tlen_param" in params:
            command_params = command_params.append("--soft-clipped-unmapped-tlen")
        if "bowtie2_sam-no-qname-trunc_param" in params:
            command_params = command_params.append("--sam-no-qname-trunc")
        if "bowtie2_xeq_param" in params:
            command_params = command_params.append("--xeq")


        return command_params

    def run(self, input_files, input_metadata, output_files):
        """
        The main function to align bam files to a genome using BWA

        Parameters
        ----------
        input_files : dict
            File 0 is the genome file location, file 1 is the FASTQ file
        metadata : dict
        output_files : dict

        Returns
        -------
        output_files : dict
            First element is a list of output_bam_files, second element is the
            matching meta data
        output_metadata : dict
        """

        sources = [input_files["genome"]]

        fqs = fastq_splitter()

        fastq1 = input_files["loc"]
        sources.append(input_files["loc"])

        fastq_file_gz = str(fastq1 + ".tar.gz")
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
            logger.info("Getting the tar file")
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

        output_bam_file = output_files["output"]
        #output_bai_file = output_files["bai"]

        logger.info("BWA ALIGNER: Aligning sequence reads to the genome")

        output_bam_list = []
        for fastq_file_pair in fastq_file_list:
            if "fastq2" in input_files:
                tmp_fq1 = gz_data_path + "/tmp/" + fastq_file_pair[0]
                tmp_fq2 = gz_data_path + "/tmp/" + fastq_file_pair[1]
                output_bam_file_tmp = tmp_fq1 + ".bam"
                output_bam_list.append(output_bam_file_tmp)

                # print("FILES:", tmp_fq1, tmp_fq2, output_bam_file_tmp)

                results = self.bowtie2_aligner_paired(
                    str(input_files["genome"]), tmp_fq1, tmp_fq2, output_bam_file_tmp,
                    str(input_files["index"]), self.get_aln_params(self.configuration)
                )
            else:
                tmp_fq = gz_data_path + "/tmp/" + fastq_file_pair[0]
                output_bam_file_tmp = tmp_fq + ".bam"
                output_bam_list.append(output_bam_file_tmp)

                logger.info("BWAL ALN FILES:" + tmp_fq)
                results = self.bowtie2_aligner_single(
                    str(input_files["genome"]), tmp_fq, output_bam_file_tmp,
                    str(input_files["index"]), self.get_aln_params(self.configuration)
                )
        barrier()

        bam_handle = bamUtilsTask()

        results = bam_handle.bam_copy(output_bam_list.pop(0), output_bam_file)
        results = compss_wait_on(results)

        if results is False:
            logger.fatal("BWA Aligner: Bam copy failed")
            return {}, {}

        while True:
            if output_bam_list:
                results = bam_handle.bam_merge(output_bam_file, output_bam_list.pop(0))
                results = compss_wait_on(results)

                if results is False:
                    logger.fatal("BWA Aligner: Bam merging failed")
                    return {}, {}
            else:
                break

        results = bam_handle.bam_sort(output_bam_file)
        results = compss_wait_on(results)

        if results is False:
            logger.fatal("BWA Aligner: Bam sorting failed")
            return {}, {}

        # results = bam_handle.bam_index(output_bam_file, output_bai_file)
        # results = compss_wait_on(results)

        # if results is False:
        #     logger.fatal("BWA Aligner: Bam indexing failed")
        #     return {}, {}

        logger.info("BWA ALIGNER: Alignments complete")

        output_metadata = {
            "bam": Metadata(
                data_type=input_metadata['loc'].data_type,
                file_type="BAM",
                file_path=output_files["output"],
                sources=[input_metadata["genome"].file_path, input_metadata['loc'].file_path],
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "bwa_aligner"
                }
            )
        }

        return ({"bam": output_files["output"]}, output_metadata)

# ------------------------------------------------------------------------------
