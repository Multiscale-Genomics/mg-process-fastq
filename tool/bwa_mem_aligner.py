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

# ------------------------------------------------------------------------------


class bwaAlignerMEMTool(Tool):  # pylint: disable=invalid-name
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
        logger.info("BWA MEM Aligner")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(returns=bool, genome_file_name=IN, genome_idx=FILE_IN,
          amb_file=FILE_OUT, ann_file=FILE_OUT, bwt_file=FILE_OUT,
          pac_file=FILE_OUT, sa_file=FILE_OUT)
    def untar_index(  # pylint: disable=too-many-locals,too-many-arguments
            self, genome_file_name, genome_idx,
            amb_file, ann_file, bwt_file, pac_file, sa_file):
        """
        Extracts the BWA index files from the genome index tar file.

        Parameters
        ----------
        genome_file_name : str
            Location string of the genome fasta file
        genome_idx : str
            Location of the BWA index file
        amb_file : str
            Location of the amb index file
        ann_file : str
            Location of the ann index file
        bwt_file : str
            Location of the bwt index file
        pac_file : str
            Location of the pac index file
        sa_file : str
            Location of the sa index file

        Returns
        -------
        bool
            Boolean indicating if the task was successful
        """
        if "no-untar" in self.configuration and self.configuration["no-untar"] is True:
            return True

        gfl = genome_file_name.split("/")
        au_handle = alignerUtils()
        au_handle.bwa_untar_index(
            gfl[-1], genome_idx, amb_file, ann_file, bwt_file, pac_file, sa_file)

        return True

    @constraint(ComputingUnits="4")
    @task(returns=bool, genome_file_loc=FILE_IN, read_file_loc=FILE_IN,
          bam_loc=FILE_OUT, amb_file=FILE_IN, ann_file=FILE_IN, bwt_file=FILE_IN,
          pac_file=FILE_IN, sa_file=FILE_IN, mem_params=IN, isModifier=False)
    def bwa_aligner_single(  # pylint: disable=too-many-arguments, no-self-use
            self, genome_file_loc, read_file_loc, bam_loc,
            amb_file, ann_file, bwt_file, pac_file, sa_file,  # pylint: disable=unused-argument
            mem_params):
        """
        BWA MEM Aligner - Single Ended

        Parameters
        ----------
        genome_file_loc : str
            Location of the genomic fasta
        read_file_loc : str
            Location of the FASTQ file
        bam_loc : str
            Location of the output aligned bam file
        amb_file : str
            Location of the amb index file
        ann_file : str
            Location of the ann index file
        bwt_file : str
            Location of the bwt index file
        pac_file : str
            Location of the pac index file
        sa_file : str
            Location of the sa index file
        mem_params : dict
            Alignment parameters

        Returns
        -------
        bam_loc : str
            Location of the output file
        """
        if (
                os.path.isfile(read_file_loc) is False or
                os.path.getsize(read_file_loc) == 0):
            return False

        out_bam = read_file_loc + '.out.bam'

        au_handle = alignerUtils()
        logger.info(
            "BWA FINISHED: " + str(au_handle.bwa_mem_align_reads(
                genome_file_loc, out_bam, mem_params, read_file_loc))
        )

        try:
            with open(bam_loc, "wb") as f_out:
                with open(out_bam, "rb") as f_in:
                    f_out.write(f_in.read())
        except (OSError, IOError) as error:
            logger.fatal("SINGLE ALIGNER: I/O error({0}): {1}".format(error.errno, error.strerror))
            return False

        os.remove(out_bam)

        return True

    @constraint(ComputingUnits="4")
    @task(returns=bool, genome_file_loc=FILE_IN, read_file_loc1=FILE_IN,
          read_file_loc2=FILE_IN, bam_loc=FILE_OUT,
          amb_file=FILE_IN, ann_file=FILE_IN, bwt_file=FILE_IN,
          pac_file=FILE_IN, sa_file=FILE_IN, mem_params=IN, isModifier=False)
    def bwa_aligner_paired(  # pylint: disable=too-many-arguments, no-self-use, too-many-locals
            self, genome_file_loc, read_file_loc1, read_file_loc2, bam_loc,
            amb_file, ann_file, bwt_file, pac_file, sa_file, mem_params):  # pylint: disable=unused-argument
        """
        BWA MEM Aligner - Paired End

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
        amb_file : str
            Location of the amb index file
        ann_file : str
            Location of the ann index file
        bwt_file : str
            Location of the bwt index file
        pac_file : str
            Location of the pac index file
        sa_file : str
            Location of the sa index file
        mem_params : dict
            Alignment parameters

        Returns
        -------
        bam_loc : str
            Location of the output file
        """
        out_bam = read_file_loc1 + '.out.bam'
        au_handle = alignerUtils()
        logger.info(
            "BWA FINISHED: " + str(au_handle.bwa_mem_align_reads(
                genome_file_loc, out_bam, mem_params, read_file_loc1, read_file_loc2))
        )

        try:
            with open(bam_loc, "wb") as f_out:
                with open(out_bam, "rb") as f_in:
                    f_out.write(f_in.read())
        except (OSError, IOError) as error:
            logger.fatal("PARIED ALIGNER: I/O error({0}): {1}".format(error.errno, error.strerror))
            return False

        os.remove(out_bam)

        return True

    @staticmethod
    def get_mem_params(params):
        """
        Function to handle to extraction of commandline parameters and formatting
        them for use in the aligner for BWA MEM

        Parameters
        ----------
        params : dict

        Returns
        -------
        list
        """

        command_parameters = {
            "bwa_mem_min_seed_len_param": ["-k", True],
            "bwa_mem_band_width_param": ["-w", True],
            "bwa_mem_zdropoff_param": ["-d", True],
            "bwa_mem_reseeding_param": ["-r", True],
            "bwa_mem_insensitive_param": ["-c", True],
            "bwa_mem_paried_rescue_mode_param": ["-P", False],
            "bwa_mem_matching_score_param": ["-A", True],
            "bwa_mem_mismatch_penalty_param": ["-B", True],
            "bwa_mem_gap_open_penalty_param": ["-O", True],
            "bwa_mem_gap_ext_penalty_param": ["-E", True],
            "bwa_mem_clipping_penalty_param": ["-L", True],
            "bwa_mem_unpaired_penalty_param": ["-U", True],
            "bwa_mem_reads_interleaved_param": ["-p", False],
            "bwa_mem_complete_read_head_param": ["-R", True],
            "bwa_mem_alignment_threshold_param": ["-T", True],
            "bwa_mem_hard_clipping_param": ["-H", False],
            "bwa_mem_short_split_secondary_param": ["-M", False]
        }

        command_params = []
        for param in params:
            if param in command_parameters:
                if command_parameters[param][1] and params[param] != "":
                    command_params = command_params + [command_parameters[param][0], params[param]]
                else:
                    if command_parameters[param][0] and params[param] is not False:
                        command_params.append(command_parameters[param][0])

        return command_params

    def run(self, input_files, input_metadata, output_files):  # pylint: disable=too-many-branches,too-many-locals,too-many-statements
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

        tasks_done = 0
        task_count = 7

        untar_idx = True
        if "no-untar" in self.configuration and self.configuration["no-untar"] is True:
            untar_idx = False
            task_count = 5

        index_files = {
            "amb": input_files["genome"] + ".amb",
            "ann": input_files["genome"] + ".ann",
            "bwt": input_files["genome"] + ".bwt",
            "pac": input_files["genome"] + ".pac",
            "sa": input_files["genome"] + ".sa"
        }

        if untar_idx:
            logger.progress("Untar Index", task_id=tasks_done, total=task_count)
            self.untar_index(
                input_files["genome"],
                input_files["index"],
                index_files["amb"],
                index_files["ann"],
                index_files["bwt"],
                index_files["pac"],
                index_files["sa"]
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

        # compss_delete_file(fastq1)
        # if "fastq2" in input_files:
        #     compss_delete_file(fastq2)

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
            os.remove(fastq_file_gz)
            compss_delete_file(fastq_file_gz)
        except tarfile.TarError:
            logger.fatal("Split FASTQ files: Malformed tar file")
            return {}, {}

        tasks_done += 1
        logger.progress("FASTQ Splitter", task_id=tasks_done, total=task_count)

        # input and output share most metadata
        output_metadata = {}

        output_bam_file = output_files["output"]
        # output_bai_file = output_files["bai"]

        logger.info("BWA ALIGNER: Aligning sequence reads to the genome")
        logger.progress("ALIGNER - jobs = " + str(len(fastq_file_list)),
                        task_id=tasks_done, total=task_count)

        output_bam_list = []
        for fastq_file_pair in fastq_file_list:
            if "fastq2" in input_files:
                tmp_fq1 = os.path.join(gz_data_path, "tmp", fastq_file_pair[0])
                tmp_fq2 = os.path.join(gz_data_path, "tmp", fastq_file_pair[1])
                output_bam_file_tmp = tmp_fq1 + ".bam"
                output_bam_list.append(output_bam_file_tmp)

                logger.info("BWA MEM FILES: " + tmp_fq1 + " - " + tmp_fq2)
                self.bwa_aligner_paired(
                    str(input_files["genome"]), tmp_fq1, tmp_fq2, output_bam_file_tmp,
                    index_files["amb"],
                    index_files["ann"],
                    index_files["bwt"],
                    index_files["pac"],
                    index_files["sa"],
                    self.get_mem_params(self.configuration)
                )
            else:
                tmp_fq = os.path.join(gz_data_path, "tmp", fastq_file_pair[0])
                output_bam_file_tmp = tmp_fq + ".bam"
                output_bam_list.append(output_bam_file_tmp)

                logger.info("BWA MEM FILES: " + tmp_fq)
                self.bwa_aligner_single(
                    str(input_files["genome"]), tmp_fq, output_bam_file_tmp,
                    index_files["amb"],
                    index_files["ann"],
                    index_files["bwt"],
                    index_files["pac"],
                    index_files["sa"],
                    self.get_mem_params(self.configuration)
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

        logger.progress("Indexing output bam file",
                        task_id=tasks_done, total=task_count)
        bam_handle.bam_index(output_bam_file, output_files["bai"])
        tasks_done += 1
        logger.progress("Indexing output bam file",
                        task_id=tasks_done, total=task_count)

        compss_delete_file(output_bam_list[0])

        logger.info("BWA ALIGNER: Alignments complete")

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
                sources=sources,
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "bwa_aligner",
                    "parameters": self.get_mem_params(self.configuration),
                    "associated_files": [output_files["bai"]]
                }
            ),
            "bai": Metadata(
                data_type=input_metadata['loc'].data_type,
                file_type="BAI",
                file_path=output_files["bai"],
                sources=sources,
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "assembly": input_metadata["genome"].meta_data["assembly"],
                    "tool": "bs_seeker_aligner",
                    "associated_master": output_bam_file
                }
            )
        }

        return (
            {"bam": output_files["output"], "bai": output_files["bai"]},
            output_metadata
        )

# ------------------------------------------------------------------------------
