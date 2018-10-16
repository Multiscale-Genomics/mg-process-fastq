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
import tarfile
import shutil

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import IN, FILE_IN, FILE_OUT
    from pycompss.api.task import task
    from pycompss.api.constraint import constraint
    from pycompss.api.api import barrier, compss_wait_on, compss_open, compss_delete_file
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import IN, FILE_IN, FILE_OUT  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task, constraint  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import barrier, compss_wait_on, compss_open, compss_delete_file  # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

from tool.fastq_splitter import fastq_splitter
from tool.aligner_utils import alignerUtils
from tool.bam_utils import bamUtilsTask
from tool.common import common

# ------------------------------------------------------------------------------


class bowtie2AlignerTool(Tool):  # pylint: disable=invalid-name
    """
    Tool for aligning sequence reads to a genome using BWA
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
        logger.info("Bowtie2 Aligner")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(returns=bool, genome_file_name=IN, genome_idx=FILE_IN,
          bt2_1_file=FILE_OUT, bt2_2_file=FILE_OUT, bt2_3_file=FILE_OUT,
          bt2_4_file=FILE_OUT, bt2_rev1_file=FILE_OUT, bt2_rev2_file=FILE_OUT)
    def untar_index(  # pylint: disable=too-many-locals,too-many-arguments
            self, genome_file_name, genome_idx,
            bt2_1_file, bt2_2_file, bt2_3_file, bt2_4_file,
            bt2_rev1_file, bt2_rev2_file):
        """
        Extracts the Bowtie2 index files from the genome index tar file.

        Parameters
        ----------
        genome_file_name : str
            Location string of the genome fasta file
        genome_idx : str
            Location of the Bowtie2 index file
        bt2_1_file : str
            Location of the <genome>.1.bt2 index file
        bt2_2_file : str
            Location of the <genome>.2.bt2 index file
        bt2_3_file : str
            Location of the <genome>.3.bt2 index file
        bt2_4_file : str
            Location of the <genome>.4.bt2 index file
        bt2_rev1_file : str
            Location of the <genome>.rev.1.bt2 index file
        bt2_rev2_file : str
            Location of the <genome>.rev.2.bt2 index file

        Returns
        -------
        bool
            Boolean indicating if the task was successful
        """
        if "no-untar" in self.configuration and self.configuration["no-untar"] is True:
            return True

        gfl = genome_file_name.split("/")
        au_handle = alignerUtils()
        au_handle.bowtie2_untar_index(
            gfl[-1], genome_idx,
            bt2_1_file, bt2_2_file, bt2_3_file, bt2_4_file,
            bt2_rev1_file, bt2_rev2_file)

        return True

    @constraint(ComputingUnits="4")
    @task(returns=bool, genome_file_loc=FILE_IN, read_file_loc=FILE_IN,
          bam_loc=FILE_OUT, bt2_1_file=FILE_IN, bt2_2_file=FILE_IN,
          bt2_3_file=FILE_IN, bt2_4_file=FILE_IN, bt2_rev1_file=FILE_IN,
          bt2_rev2_file=FILE_IN, aln_params=IN, isModifier=False)
    def bowtie2_aligner_single(  # pylint: disable=too-many-arguments, no-self-use
            self, genome_file_loc, read_file_loc, bam_loc,
            bt2_1_file, bt2_2_file, bt2_3_file, bt2_4_file,  # pylint: disable=unused-argument
            bt2_rev1_file, bt2_rev2_file, aln_params):  # pylint: disable=unused-argument
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
        bt2_1_file : str
            Location of the <genome>.1.bt2 index file
        bt2_2_file : str
            Location of the <genome>.2.bt2 index file
        bt2_3_file : str
            Location of the <genome>.3.bt2 index file
        bt2_4_file : str
            Location of the <genome>.4.bt2 index file
        bt2_rev1_file : str
            Location of the <genome>.rev.1.bt2 index file
        bt2_rev2_file : str
            Location of the <genome>.rev.2.bt2 index file
        aln_params : dict
            Alignment parameters

        Returns
        -------
        bam_loc : str
            Location of the output file
        """
        out_bam = read_file_loc + '.out.bam'

        au_handle = alignerUtils()
        logger.info(
            "BOWTIE2 FINISHED: " + str(au_handle.bowtie2_align_reads(
                genome_file_loc, out_bam, aln_params, read_file_loc))
        )

        common_handle = common()
        return_val = common_handle.to_output_file(out_bam, bam_loc, False)

        if return_val is False:
            logger.fatal("IO Error: Missing file - {}".format(out_bam))

        os.remove(out_bam)

        return return_val

    @constraint(ComputingUnits="4")
    @task(returns=bool, genome_file_loc=FILE_IN, read_file_loc1=FILE_IN,
          read_file_loc2=FILE_IN, bam_loc=FILE_OUT, bt2_1_file=FILE_IN,
          bt2_2_file=FILE_IN, bt2_3_file=FILE_IN, bt2_4_file=FILE_IN,
          bt2_rev1_file=FILE_IN, bt2_rev2_file=FILE_IN,
          aln_params=IN, isModifier=False)
    def bowtie2_aligner_paired(  # pylint: disable=too-many-arguments, no-self-use
            self, genome_file_loc, read_file_loc1, read_file_loc2, bam_loc,
            bt2_1_file, bt2_2_file, bt2_3_file, bt2_4_file,  # pylint: disable=unused-argument
            bt2_rev1_file, bt2_rev2_file, aln_params):  # pylint: disable=unused-argument
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
        bt2_1_file : str
            Location of the <genome>.1.bt2 index file
        bt2_2_file : str
            Location of the <genome>.2.bt2 index file
        bt2_3_file : str
            Location of the <genome>.3.bt2 index file
        bt2_4_file : str
            Location of the <genome>.4.bt2 index file
        bt2_rev1_file : str
            Location of the <genome>.rev.1.bt2 index file
        bt2_rev2_file : str
            Location of the <genome>.rev.2.bt2 index file
        aln_params : dict
            Alignment parameters

        Returns
        -------
        bam_loc : str
            Location of the output file
        """
        out_bam = read_file_loc1 + '.out.bam'
        au_handle = alignerUtils()
        logger.info(
            "BOWTIE2 FINISHED: " + str(au_handle.bowtie2_align_reads(
                genome_file_loc, out_bam, aln_params, read_file_loc1, read_file_loc2))
        )

        try:
            with open(bam_loc, "wb") as f_out:
                with open(out_bam, "rb") as f_in:
                    f_out.write(f_in.read())
        except (OSError, IOError):
            return False

        os.remove(out_bam)

        return True

    @staticmethod
    def get_aln_params(params, paired=False):
        """
        Function to handle to extraction of commandline parameters and formatting
        them for use in the aligner for Bowtie2

        Parameters
        ----------
        params : dict
        paired : bool
            Indicate if the parameters are paired-end specific. [DEFAULT=False]

        Returns
        -------
        list
        """
        command_params = ["-q"]

        command_parameters = {
            # Input Options - 11
            "bowtie2_interleaved_param": ["--interleaved", False],
            "bowtie2_tab5_param": ["--tab5", False],
            "bowtie2_tab6_param": ["--tab6", False],
            # "bowtie2_qseq_param": ["--qseq", False],
            # "bowtie2_read_only_param": ["-r", True],
            "bowtie2_skip_1st_n_reads_param": ["-s", True],
            "bowtie2_aln_1st_n_reads_param": ["-u", True],
            "bowtie2_trim5_param": ["-5", True],
            "bowtie2_trim3_param": ["-3", True],
            "bowtie2_phred33_param": ["--phred33", False],
            "bowtie2_phre64_param": ["--phred64", False],
            # Alignment Options - 12
            "bowtie2_num_mismatch_param": ["-N", True],
            "bowtie2_seed_len_param": ["-L", True],
            # "bowtie2_seed_func_param": ["-i", True],
            # "bowtie2_ambg_char_func_param": ["--n-cell", True],
            "bowtie2_dpads_param": ["--dpad", True],
            "bowtie2_gbar_param": ["--gbar", True],
            "bowtie2_ignore_quals_param": ["--ignore-quals", False],
            "bowtie2_nofw_param": ["--nofw", False],
            "bowtie2_norc_param": ["--norc", False],
            "bowtie2_no_1mm_upfront_param": ["--no-1mm-upfront", False],
            "bowtie2_end_to_end_param": ["--end-to-end", False],
            "bowtie2_local_param": ["--local", False],
            # Effort Options - 2
            "bowtie2_seed_extension_attempts_param": ["-D", True],
            "bowtie2_reseed_param": ["-R", True],
            # SAM Options - 9
            "bowtie2_no_unal_param": ["--no-unal", False],
            "bowtie2_no_hd_param": ["--no-hd", False],
            "bowtie2_no_sq_param": ["--no-dq", False],
            "bowtie2_rg_id_param": ["--rg-id", True],
            "bowtie2_rg_param": ["--rg", True],
            "bowtie2_omit_sec_seq_param": ["--omit-sec-seq", False],
            "bowtie2_soft_clipped_unmapped_tlen_param": ["--soft-clipped-unmapped-tlen", False],
            "bowtie2_sam_no_qname_trunc_param": ["--sam-no-qname-trunc", False],
            "bowtie2_xeq_param": ["--xeq", False],
        }

        if paired:
            # Paired-end Options - 10
            command_parameters["bowtie2_min_frag_len_param"] = ["-I", True]
            command_parameters["bowtie2_max_frag_len_param"] = ["-X", True]
            command_parameters["bowtie2_frrfff_param"] = ["", True]
            command_parameters["bowtie2_no_mixed_param"] = ["--no-mixed", False]
            command_parameters["bowtie2_no_discordant_param"] = ["--no-discordant", False]
            command_parameters["bowtie2_dovetail_param"] = ["--dovetail", False]
            command_parameters["bowtie2_no_contain_param"] = ["--no-contain", False]
            command_parameters["bowtie2_no_overlap_param"] = ["--no-overlap", False]

        for param in params:
            if param in command_parameters:
                if command_parameters[param][1] and params[param] != "":
                    command_params = command_params + [command_parameters[param][0], params[param]]
                else:
                    if command_parameters[param][0] and params[param] is not False:
                        command_params.append(command_parameters[param][0])

        # Scoring Options - 8
        if "bowtie2_ma_param" in params and params["bowtie2_ma_param"] != "":
            command_params = command_params + [
                "--ma", str(params["bowtie2_ma_param"])]
        if ("bowtie2_mp_mx_param" in params and "bowtie2_mp_mn_param" in params and
                params["bowtie2_mp_mx_param"] != "" and params["bowtie2_mp_mn_param"] != ""):
            command_params = command_params + [
                "--mp",
                str(params["bowtie2_mp_mx_param"]) + "," + str(params["bowtie2_mp_mn_param"])]
        if "bowtie2_np_param" in params and params["bowtie2_np_param"] != "":
            command_params = command_params + [
                "--np", str(params["bowtie2_np_param"])]
        if ("bowtie2_rdg_o_param" in params and "bowtie2_rdg_e_param" in params and
                params["bowtie2_rdg_o_param"] != "" and params["bowtie2_rdg_e_param"] != ""):
            command_params = command_params + [
                "--rdg",
                str(params["bowtie2_rdg_o_param"]) + "," + str(params["bowtie2_rdg_e_param"])]
        if ("bowtie2_rfg_o_param" in params and "bowtie2_rfg_e_param" in params and
                params["bowtie2_rfg_o_param"] != "" and params["bowtie2_rfg_e_param"] != ""):
            command_params = command_params + [
                "--rfg",
                str(params["bowtie2_rfg_o_param"]) + "," + str(params["bowtie2_rfg_e_param"])]
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

    def run(self, input_files, input_metadata, output_files):  # pylint: disable=too-many-branches,too-many-locals,too-many-statements
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

        tasks_done = 0
        task_count = 6

        untar_idx = True
        if "no-untar" in self.configuration and self.configuration["no-untar"] is True:
            untar_idx = False
            task_count = 5

        index_files = {
            "1.bt2": input_files["genome"] + ".1.bt2",
            "2.bt2": input_files["genome"] + ".2.bt2",
            "3.bt2": input_files["genome"] + ".3.bt2",
            "4.bt2": input_files["genome"] + ".4.bt2",
            "rev.1.bt2": input_files["genome"] + ".rev.1.bt2",
            "rev.2.bt2": input_files["genome"] + ".rev.2.bt2"
        }

        if untar_idx:
            logger.progress("Untar Index", task_id=tasks_done, total=task_count)
            self.untar_index(
                input_files["genome"],
                input_files["index"],
                index_files["1.bt2"],
                index_files["2.bt2"],
                index_files["3.bt2"],
                index_files["4.bt2"],
                index_files["rev.1.bt2"],
                index_files["rev.2.bt2"]
            )
            tasks_done += 1
            logger.progress("Untar Index", task_id=tasks_done, total=task_count)

        sources = [input_files["genome"]]

        fqs = fastq_splitter()

        fastq1 = input_files["loc"]
        sources.append(input_files["loc"])

        logger.progress("FASTQ Splitter", task_id=tasks_done, total=task_count)
        fastq_file_gz = os.path.join(
            self.configuration["execution"], os.path.split(fastq1)[1] + ".tar.gz")
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

        # Required to prevent iterating over the future objects
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

        gz_data_path = os.path.split(fastq_file_gz)[0]

        try:
            tar = tarfile.open(fastq_file_gz)
            tar.extractall(path=gz_data_path)
            tar.close()
            compss_delete_file(fastq_file_gz)
            try:
                os.remove(fastq_file_gz)
            except (OSError, IOError) as msg:
                logger.warn(
                    "Unable to remove file I/O error({0}): {1}".format(
                        msg.errno, msg.strerror
                    )
                )
        except tarfile.TarError:
            logger.fatal("Split FASTQ files: Malformed tar file")
            return {}, {}

        tasks_done += 1
        logger.progress("FASTQ Splitter", task_id=tasks_done, total=task_count)

        # input and output share most metadata
        output_metadata = {}

        output_bam_file = output_files["output"]

        logger.info("BOWTIE2 ALIGNER: Aligning sequence reads to the genome")
        logger.progress("ALIGNER - jobs = " + str(len(fastq_file_list)),
                        task_id=tasks_done, total=task_count)

        output_bam_list = []
        for fastq_file_pair in fastq_file_list:
            if "fastq2" in input_files:
                command_params = self.get_aln_params(self.configuration, True)
                tmp_fq1 = os.path.join(gz_data_path, "tmp", fastq_file_pair[0])
                tmp_fq2 = os.path.join(gz_data_path, "tmp", fastq_file_pair[1])
                output_bam_file_tmp = tmp_fq1 + ".bam"
                output_bam_list.append(output_bam_file_tmp)

                logger.info("BOWTIE2 ALIGN (PAIRED) FILES:\n\t{}\n\t{}".format(tmp_fq1, tmp_fq2))

                self.bowtie2_aligner_paired(
                    str(input_files["genome"]), tmp_fq1, tmp_fq2,
                    output_bam_file_tmp,
                    index_files["1.bt2"],
                    index_files["2.bt2"],
                    index_files["3.bt2"],
                    index_files["4.bt2"],
                    index_files["rev.1.bt2"],
                    index_files["rev.2.bt2"],
                    command_params
                )
            else:
                command_params = self.get_aln_params(self.configuration)
                tmp_fq = os.path.join(gz_data_path, "tmp", fastq_file_pair[0])
                output_bam_file_tmp = tmp_fq + ".bam"
                output_bam_list.append(output_bam_file_tmp)

                logger.info("BOWTIE2 ALIGN (SINGLE) FILES:" + tmp_fq)
                self.bowtie2_aligner_single(
                    str(input_files["genome"]), tmp_fq, output_bam_file_tmp,
                    index_files["1.bt2"],
                    index_files["2.bt2"],
                    index_files["3.bt2"],
                    index_files["4.bt2"],
                    index_files["rev.1.bt2"],
                    index_files["rev.2.bt2"],
                    command_params
                )

        barrier()

        # Remove all tmp fastq files now that the reads have been aligned
        if untar_idx:
            for idx_file in index_files:
                compss_delete_file(index_files[idx_file])

        if hasattr(sys, '_run_from_cmdl') is True:
            pass
        else:
            for fastq_file_pair in fastq_file_list:
                tmp_fq = os.path.join(gz_data_path, "tmp", fastq_file_pair[0])
                compss_delete_file(tmp_fq)
                try:
                    os.remove(tmp_fq)
                except (OSError, IOError) as msg:
                    logger.warn(
                        "Unable to remove file I/O error({0}): {1}".format(
                            msg.errno, msg.strerror
                        )
                    )
                if "fastq2" in input_files:
                    tmp_fq = os.path.join(gz_data_path, "tmp", fastq_file_pair[1])
                    compss_delete_file(tmp_fq)
                    try:
                        os.remove(tmp_fq)
                    except (OSError, IOError) as msg:
                        logger.warn(
                            "Unable to remove file I/O error({0}): {1}".format(
                                msg.errno, msg.strerror
                            )
                        )

        tasks_done += 1
        logger.progress("ALIGNER", task_id=tasks_done, total=task_count)

        bam_handle = bamUtilsTask()

        logger.progress("Merging bam files", task_id=tasks_done, total=task_count)
        bam_handle.bam_merge(output_bam_list)
        tasks_done += 1
        logger.progress("Merging bam files", task_id=tasks_done, total=task_count)

        # Remove all bam files that are not the final file
        for i in output_bam_list[1:len(output_bam_list)]:
            try:
                compss_delete_file(i)
                os.remove(i)
            except (OSError, IOError) as msg:
                logger.warn(
                    "Unable to remove file I/O error({0}): {1}".format(
                        msg.errno, msg.strerror
                    )
                )

        logger.progress("Sorting merged bam file", task_id=tasks_done, total=task_count)
        bam_handle.bam_sort(output_bam_list[0])
        tasks_done += 1
        logger.progress("Sorting merged bam file", task_id=tasks_done, total=task_count)

        logger.progress("Copying bam file into the output file",
                        task_id=tasks_done, total=task_count)
        bam_handle.bam_copy(output_bam_list[0], output_bam_file)
        tasks_done += 1
        logger.progress("Copying bam file into the output file",
                        task_id=tasks_done, total=task_count)

        compss_delete_file(output_bam_list[0])

        logger.info("BOWTIE2 ALIGNER: Alignments complete")

        barrier()
        try:
            shutil.rmtree(gz_data_path + "/tmp")
        except (OSError, IOError) as msg:
            logger.warn(
                "Already tidy I/O error({0}): {1}".format(
                    msg.errno, msg.strerror
                )
            )

        output_metadata = {
            "bam": Metadata(
                data_type=input_metadata['loc'].data_type,
                file_type="BAM",
                file_path=output_files["output"],
                sources=[input_metadata["genome"].file_path, input_metadata['loc'].file_path],
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "bowtie_aligner",
                    "parameters": command_params
                }
            )
        }

        return ({"bam": output_files["output"]}, output_metadata)

# ------------------------------------------------------------------------------
