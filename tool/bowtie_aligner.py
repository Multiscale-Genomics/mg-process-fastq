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
    from pycompss.api.parameter import IN, FILE_IN, FILE_OUT
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on, compss_open, barrier
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import IN, FILE_IN, FILE_OUT # pylint: disable=ungrouped-imports
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
        Bowtie2 Aligner - Single End

        Parameters
        ----------
        genome_file_loc : str
            Location of the genomic fasta
        read_file_loc1 : str
            Location of the FASTQ file
        bam_loc : str
            Location of the output aligned bam file
        genome_idx : idx
            Location of the Bowtie2 index file
        aln_params : dict
            Alignment parameters

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
                genome_fa_ln, out_bam, aln_params, read_file_loc))
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
        Bowtie2 Aligner - Paired End

        Parameters
        ----------
        genome_file_loc : str
            Location of the genomic fasta
        read_file_loc1 : str
            Location of the FASTQ file
        read_file_loc2 : str
            Location of the FASTQ file
        bam_loc : str
            Location of the output aligned bam file
        genome_idx : idx
            Location of the Bowtie2 index file
        aln_params : dict
            Alignment parameters

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
                genome_fa_ln, out_bam, aln_params, read_file_loc1, read_file_loc2))
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
        Function to handle to extraction of commandline parameters and formatting
        them for use in the aligner for Bowtie2

        Parameters
        ----------
        params : dict

        Returns
        -------
        list
        """
        command_params = ["-q"]

        command_parameters = {
            # Input Options - 11
            "bowtie2_interleaved_param" : ["--interleaved", False],
            "bowtie2_tab5_param" : ["--tab5", False],
            "bowtie2_tab6_param" : ["--tab6", False],
            "bowtie2_qseq_param" : ["--qseq", False],
            "bowtie2_read_only_param" : ["-r", True],
            "bowtie2_skip_1st_n_reads_param" : ["-s", True],
            "bowtie2_aln_1st_n_reads_param" : ["-u", True],
            "bowtie2_trim5_param" : ["-5", True],
            "bowtie2_trim3_param" : ["-3", True],
            "bowtie2_phred33_param" : ["--phred33", False],
            "bowtie2_phre64_param" : ["--phred64", False],
            # Alignment Options - 12
            "bowtie2_num_mismatch_param" : ["-N", True],
            "bowtie2_seed_len_param" : ["-L", True],
            "bowtie2_seed_func_param" : ["-i", True],
            "bowtie2_ambg_char_func_param" : ["--n-cell", True],
            "bowtie2_dpads_param" : ["--dpad", True],
            "bowtie2_gbar_param" : ["--gbar", True],
            "bowtie2_ignore_quals_param" : ["--ignore-quals", False],
            "bowtie2_nofw_param" : ["--nofw", False],
            "bowtie2_norc_param" : ["--norc", False],
            "bowtie2_no_1mm_upfront_param" : ["--no-1mm-upfront", False],
            "bowtie2_end_to_end_param" : ["--end-to-end", False],
            "bowtie2_local_param" : ["--local", False],
            # Effort Options - 2
            "bowtie2_seed-extension-attempts_param" : ["-D", True],
            "bowtie2_reseed_param" : ["-R", True],
            # Paired-end Options - 10
            "bowtie2_min_frag_len_param" : ["-I", True],
            "bowtie2_max_frag_len_param" : ["-X", True],
            "bowtie2_fr_param" : ["--fr", False],
            "bowtie2_rf_param" : ["--rf", False],
            "bowtie2_ff_param" : ["--ff", False],
            "bowtie2_no_mixed_param" : ["--no-mixed", False],
            "bowtie2_no_discordant_param" : ["--no-discordant", False],
            "bowtie2_dovetail_param" : ["--dovetail", False],
            "bowtie2_no_contain_param" : ["--no-contain", False],
            "bowtie2_no_overlap_param" : ["--no-overlap", False],
            # SAM Options - 9
            "bowtie2_no_unal_param" : ["--no-unal", False],
            "bowtie2_no_hd_param" : ["--no-hd", False],
            "bowtie2_no_sq_param" : ["--no-dq", False],
            "bowtie2_rg_id_param" : ["--rg-id", True],
            "bowtie2_rg_param" : ["--rg", True],
            "bowtie2_omit_sec_seq_param" : ["--omit-sec-seq", False],
            "bowtie2_soft_clipped_unmapped_tlen_param" : ["--soft-clipped-unmapped-tlen", False],
            "bowtie2_sam_no_qname_trunc_param" : ["--sam-no-qname-trunc", False],
            "bowtie2_xeq_param" : ["--xeq", False],
        }

        for param in params:
            if param in command_parameters:
                if command_parameters[param][1]:
                    command_params = command_params + [command_parameters[param][0], params[param]]
                else:
                    command_params.append(command_parameters[param][0])

        # Scoring Options - 8
        if "bowtie2_ma_param" in params:
            command_params = command_params + [
                "--ma_", str(params["bowtie2_ma_param"])]
        if "bowtie2_mp-mx_param" in params and "bowtie2_mp-mn_param" in params:
            command_params = command_params + [
                "--mp",
                str(params["bowtie2_mp-mx_param"]) + "," + str(params["bowtie2_mp-mn_param"])]
        if "bowtie2_np_param" in params:
            command_params = command_params + [
                "--np", str(params["bowtie2_np_param"])]
        if "bowtie2_rdg-o_param" in params and "bowtie2_rdg-e_param" in params:
            command_params = command_params + [
                "--rdg",
                str(params["bowtie2_rdg-o_param"]) + "," + str(params["bowtie2_rdg-e_param"])]
        if "bowtie2_rfg_param" in params:
            command_params = command_params + [
                "--rfg",
                str(params["bowtie2_rfg-o_param"]) + "," + str(params["bowtie2_rfg-e_param"])]
        # if "bowtie2_score-min_param" in params:
        #     command_params = command_params + [
        #         "--score-min", str(params["bowtie2_score-min_param"])]

        # Reporting Options
        # if "bowtie2_reporting-k_param" in params:
        #     command_params = command_params + [
        #         "-k", str(params["bowtie2_reporting-k_param"])]
        # if "bowtie2_reporting-a_param" in params:
        #     command_params = command_params.append("-a")



        return command_params

    def run(self, input_files, input_metadata, output_files):
        """
        The main function to align bam files to a genome using Bowtie2

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

        logger.info("BOWTIE2 ALIGNER: Aligning sequence reads to the genome")

        output_bam_list = []
        for fastq_file_pair in fastq_file_list:
            if "fastq2" in input_files:
                tmp_fq1 = gz_data_path + "/tmp/" + fastq_file_pair[0]
                tmp_fq2 = gz_data_path + "/tmp/" + fastq_file_pair[1]
                output_bam_file_tmp = tmp_fq1 + ".bam"
                output_bam_list.append(output_bam_file_tmp)

                results = self.bowtie2_aligner_paired(
                    str(input_files["genome"]), tmp_fq1, tmp_fq2, output_bam_file_tmp,
                    str(input_files["index"]), self.get_aln_params(self.configuration)
                )
            else:
                tmp_fq = gz_data_path + "/tmp/" + fastq_file_pair[0]
                output_bam_file_tmp = tmp_fq + ".bam"
                output_bam_list.append(output_bam_file_tmp)

                logger.info("BOWTIE2 ALN FILES:" + tmp_fq)
                results = self.bowtie2_aligner_single(
                    str(input_files["genome"]), tmp_fq, output_bam_file_tmp,
                    str(input_files["index"]), self.get_aln_params(self.configuration)
                )
        barrier()

        bam_handle = bamUtilsTask()

        results = bam_handle.bam_copy(output_bam_list.pop(0), output_bam_file)
        results = compss_wait_on(results)

        if results is False:
            logger.fatal("BOWTIE2 Aligner: Bam copy failed")
            return {}, {}

        while True:
            if output_bam_list:
                results = bam_handle.bam_merge(output_bam_file, output_bam_list.pop(0))
                results = compss_wait_on(results)

                if results is False:
                    logger.fatal("BOWTIE2 Aligner: Bam merging failed")
                    return {}, {}
            else:
                break

        results = bam_handle.bam_sort(output_bam_file)
        results = compss_wait_on(results)

        if results is False:
            logger.fatal("BOWTIE2 Aligner: Bam sorting failed")
            return {}, {}

        # results = bam_handle.bam_index(output_bam_file, output_bai_file)
        # results = compss_wait_on(results)

        # if results is False:
        #     logger.fatal("BOWTIE2 Aligner: Bam indexing failed")
        #     return {}, {}

        logger.info("BOWTIE2 ALIGNER: Alignments complete")

        output_metadata = {
            "bam": Metadata(
                data_type=input_metadata['loc'].data_type,
                file_type="BAM",
                file_path=output_files["output"],
                sources=[input_metadata["genome"].file_path, input_metadata['loc'].file_path],
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "bowtie_aligner"
                }
            )
        }

        return ({"bam": output_files["output"]}, output_metadata)

# ------------------------------------------------------------------------------
