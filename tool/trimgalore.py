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

from utils import logger

try:
    if hasattr(sys, "_run_from_cmdl") is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import compss_wait_on  # pylint: disable=ungrouped-imports

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool


# ------------------------------------------------------------------------------

class trimgalore(Tool):
    """
    Tool for trimming FASTQ reads that are of low quality
    """

    def __init__(self, configuration=None):
        """
        Init function
        """
        logger.info("TrimGalore FASTQ read trimming")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(returns=bool, fastq_file_in=FILE_IN, fastq_file_out=FILE_OUT,
          params=IN, isModifier=False)
    def trimgalore_single(self, fastq_file_in, fastq_file_out, params):  # pylint: disable=no-self-use
        """
        Trims and removes low quality subsections and reads from a singed-ended
        FASTQ file

        Parameters
        ----------
        fastq_file_in : str
            Location of the input fastq file
        fastq_file_out : str
            Location of the output fastq file
        params : dict
            Parameters to use in TrimGalore

        Returns
        -------
        bool
            Indicator of the success of the function
        """

        command_line = "trim_galore " + " ".join(params) + " " + fastq_file_in
        logger.info("TRIM GALORE: command_line: " + command_line)

        try:
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()
        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}) - trim_galore: {1}\n{2}".format(
                msg.errno, msg.strerror, command_line))
            return False

        # Output file name used by TrimGalore
        tg_tmp_out = fastq_file_in.replace(".fastq", "_trimmed.fq")

        try:
            with open(fastq_file_out, "wb") as f_out:
                with open(tg_tmp_out, "rb") as f_in:
                    f_out.write(f_in.read())
        except IOError as error:
            logger.fatal(
                "I/O error({0}): {1}".format(error.errno, error.strerror))
            return False

        return True

    @task(returns=bool,
          fastq1_file_in=FILE_IN, fastq1_file_out=FILE_OUT,
          fastq2_file_in=FILE_IN, fastq2_file_out=FILE_OUT,
          params=IN, isModifier=False)
    def trimgalore_paired(  # pylint: disable=too-many-arguments
            self, fastq1_file_in, fastq1_file_out,
            fastq2_file_in, fastq2_file_out, params):  # pylint: disable=no-self-use
        """
        Trims and removes low quality subsections and reads from paired-end
        FASTQ files

        Parameters
        ----------
        fastq_file_in : str
            Location of the input fastq file
        fastq_file_out : str
            Location of the output fastq file
        params : dict
            Parameters to use in TrimGalore

        Returns
        -------
        bool
            Indicator of the success of the function
        """

        command_line = "trim_galore " + " ".join(params) + " " + fastq1_file_in + " " + fastq2_file_in
        logger.info("TRIM GALORE: command_line: " + command_line)

        try:
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()
        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}) - trim_galore: {1}\n{2}".format(
                msg.errno, msg.strerror, command_line))
            return False

        # Output file name used by TrimGalore
        tg_tmp_out_1 = fastq1_file_in.replace(".fastq", "_trimmed.fq")
        tg_tmp_out_2 = fastq2_file_in.replace(".fastq", "_trimmed.fq")

        try:
            with open(fastq1_file_out, "wb") as f_out:
                with open(tg_tmp_out_1, "rb") as f_in:
                    f_out.write(f_in.read())
        except IOError as error:
            logger.fatal("I/O error({0}): {1}".format(error.errno, error.strerror))
            return False

        try:
            with open(fastq2_file_out, "wb") as f_out:
                with open(tg_tmp_out_2, "rb") as f_in:
                    f_out.write(f_in.read())
        except IOError as error:
            logger.fatal("I/O error({0}): {1}".format(error.errno, error.strerror))
            return False

        return True

    @staticmethod
    def get_trimgalore_params(params):
        """
        Function to handle for extraction of commandline parameters

        Parameters
        ----------
        params : dict

        Returns
        -------
        list
        """
        command_params = []

        command_parameters = {
            # General options
            "tg_quality": ["--quality", True],
            "tg_fastqc": ["--fastqc", False],
            "tg_fastqc_args": ["--fastqc_args", True],
            "tg_adapter": ["--adapter", True],
            "tg_adapter2": ["--adapter2", True],
            "tg_illumina": ["--illumina", False],
            "tg_nextera": ["--nextera", False],
            "tg_small_rna": ["--small_rna", False],
            "tg_max_length": ["--max_length", True],
            "tg_stringency": ["--stringency", True],
            "tg_error_rate": ["-e", True],
            "tg_length": ["--length", True],
            "tg_max_n": ["--max_n", True],
            "tg_trim_n": ["--trim-n", False],
            "tg_output_dir": ["--output_dir", True],
            "tg_no_report_file": ["--no_report_file", False],
            "tg_clip_R1": ["--clip_R1", True],
            "tg_clip_R2": ["--clip_R2", True],
            "tg_three_prime_clip_R1": ["--three_prime_clip_R1", True],
            "tg_three_prime_clip_R2": ["--three_prime_clip_R2", True],
            # RRBS specific options
            "tg_rrbs": ["--rrbs", False],
            "tg_non_directional": ["--non_directional", False],
            "tg_keep": ["--keep", False],
            # Paired-end specific options
            "tg_paired": ["--paired", False],
            "tg_trim1": ["--trim1", False],
            "tg_retain_unpaired": ["--retain_unpaired", False],
            "tg_length_1": ["--length_1", True],
            "tg_length_2": ["--length_2", True],
        }

        for param in params:
            if param in command_parameters:
                if command_parameters[param][1]:
                    command_params = command_params + [command_parameters[param][0], params[param]]
                else:
                    if command_parameters[param][0]:
                        command_params.append(command_parameters[param][0])

        if "tg_phred33" in params and "tg_phred64" not in params:
            command_params.append(command_parameters["tg_phred33"][0])
        if "tg_phred64" in params and "tg_phred33" not in params:
            command_params.append(command_parameters["tg_phred64"][0])

        return command_params

    def run(self, input_files, input_metadata, output_files):
        """
        The main function to run TrimGalore to remove low quality and very short
        reads. TrimGalore uses CutAdapt and FASTQC for the analysis.

        Parameters
        ----------
        input_files : dict
            fastq1 : string
                Location of the FASTQ file
            fastq2 : string
                [OPTIONAL] Location of the paired end FASTQ file
        metadata : dict
            Matching metadata for the inpit FASTQ files

        Returns
        -------
        output_files : dict
            fastq1_trimmed : str
                Location of the trimmed FASTQ file
            fastq2_trimmed : str
                [OPTIONAL] Location of a trimmed paired end FASTQ file
        output_metadata : dict
            Matching metadata for the output files

        """
        command_params = self.get_trimgalore_params(self.configuration)

        if "fastq2" in input_files:
            pass
            # results = self.trimgalore_paired(
            #     input_files['input'], output_files['output'], command_params)
        else:
            results = self.trimgalore_single(
                input_files['fastq1'], output_files["fastq1_trimmed"], command_params)
        results = compss_wait_on(results)

        if results is False:
            logger.fatal("TrimGalore: run failed with error: {}", results)
            return ({}, {})

        output_files_created = {"fastq1_trimmed" : output_files["fastq1_trimmed"]}

        output_metadata = {
            "fastq1_out": Metadata(
                data_type=input_metadata["fastq1"].data_type,
                file_type="FASTQ",
                file_path=output_files["fastq1_trimmed"],
                sources=[input_metadata["fastq1"].file_path],
                taxon_id=input_metadata["fastq1"].taxon_id,
                meta_data={
                    "tool": "trim_galore"
                }
            )
        }

        if "fastq2" in input_files:
            output_files_created["fastq2_trimmed"] = output_files["fastq2_trimmed"]
            output_metadata = {
                "fastq2_out": Metadata(
                    data_type=input_metadata["fastq2"].data_type,
                    file_type="FASTQ",
                    file_path=output_files["fastq2_trimmed"],
                    sources=[input_metadata["fastq2"].file_path],
                    taxon_id=input_metadata["fastq2"].taxon_id,
                    meta_data={
                        "tool": "trim_galore"
                    }
                )
            }

        logger.info("TRIM GALORE: GENERATED FILES:", output_files)

        return (output_files_created, output_metadata)

# ------------------------------------------------------------------------------
