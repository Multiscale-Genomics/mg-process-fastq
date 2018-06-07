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
import itertools
import sys
import os
import tarfile

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    # from pycompss.api.api import compss_wait_on
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task  # pylint: disable=ungrouped-imports
    # from utils.dummy_pycompss import compss_wait_on  # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

# ------------------------------------------------------------------------------


class kallistoQuantificationTool(Tool):  # pylint: disable=invalid-name
    """
    Tool for quantifying RNA-seq alignments to calculate expression levels of
    genes within a genome.
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
        logger.info("Kallisto Quantification")
        Tool.__init__(self)

        if configuration is None:
            configuration = {
                "kallisto_bootstrap_param": 0
            }

        if "kallisto_bootstrap_param" not in configuration:
            configuration["kallisto_bootstrap_param"] = 0

        self.configuration.update(configuration)

    @task(
        cdna_idx_file=FILE_IN,
        fastq_file_loc=FILE_IN,
        kallisto_tar_file=FILE_OUT,
        bootstrap=IN)
    def kallisto_quant_single(
            self, cdna_idx_file, fastq_file_loc, kallisto_tar_file, bootstrap=0):
        """
        Kallisto quantifier for single end RNA-seq data

        Parameters
        ----------
        idx_loc : str
            Location of the output index file
        fastq_file_loc : str
            Location of the FASTQ sequence file

        Returns
        -------
        wig_file_loc : loc
            Location of the wig file containing the levels of expression
        """

        fq_stats = self.seq_read_stats(fastq_file_loc)

        output_dir = fastq_file_loc.split('/')

        std = fq_stats["std"]
        if std == 0.0:
            std = 1/fq_stats['mean']

        command_line = "kallisto quant -i " + cdna_idx_file + " "
        if bootstrap:
            command_line += "--bootstrap-samples=" + str(bootstrap) + " "
        command_line += "-t 4"
        command_line += " -o " + "/".join(output_dir[0:-1]) + "/"
        command_line += " --single -l " + str(fq_stats['mean']) + " "
        command_line += "-s " + str(std) + " " + fastq_file_loc

        logger.info("KALLISTO_QUANT COMMAND", command_line)

        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        self.compress_results(
            kallisto_tar_file,
            '/'.join(output_dir[0:-1]) + "/abundance.h5",
            '/'.join(output_dir[0:-1]) + "/abundance.tsv",
            '/'.join(output_dir[0:-1]) + "/run_info.json"
        )

        return True

    @task(
        fastq_file_loc_01=FILE_IN,
        fastq_file_loc_02=FILE_IN,
        cdna_idx_file=FILE_IN,
        kallisto_tar_file=FILE_OUT,
        bootstrap=IN)
    def kallisto_quant_paired(
            self, cdna_idx_file, fastq_file_loc_01, fastq_file_loc_02,
            kallisto_tar_file, bootstrap=0):
        """
        Kallisto quantifier for paired end RNA-seq data

        Parameters
        ----------
        idx_loc : str
            Location of the output index file
        fastq_file_loc_01 : str
            Location of the FASTQ sequence file
        fastq_file_loc_02 : str
            Location of the paired FASTQ sequence file

        Returns
        -------
        wig_file_loc : loc
            Location of the wig file containing the levels of expression
        """

        output_dir = fastq_file_loc_01.split('/')

        command_line = 'kallisto quant -i ' + cdna_idx_file + ' '
        if bootstrap:
            command_line += "--bootstrap-samples=" + str(bootstrap) + " "
        command_line += "-t 4"
        command_line += " -o " + "/".join(output_dir[0:-1]) + "/ "
        command_line += fastq_file_loc_01 + ' ' + fastq_file_loc_02

        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        logger.info("OUTPUT DIR: " + '/'.join(output_dir[0:-1]))

        self.compress_results(
            kallisto_tar_file,
            '/'.join(output_dir[0:-1]) + "/abundance.h5",
            '/'.join(output_dir[0:-1]) + "/abundance.tsv",
            '/'.join(output_dir[0:-1]) + "/run_info.json"
        )

        return True

    @staticmethod
    def seq_read_stats(file_in):
        """
        Calculate the mean and standard deviation of the reads in a fastq file

        Parameters
        ----------
        file_in : str
            Location of a FASTQ file

        Returns
        -------
        dict
            mean : Mean length of sequenced strands
            std  : Standard deviation of lengths of sequenced strands

        """

        from numpy import std
        from numpy import mean

        total_len = []
        with open(file_in, 'r') as file_handle:
            forthlines = itertools.islice(file_handle, 1, None, 4)
            for line in forthlines:
                line = line.rstrip()
                total_len.append(len(line))

        length_sd = std(total_len)
        length_mean = mean(total_len)

        return {'mean': length_mean, 'std': length_sd}

    @staticmethod
    def compress_results(kallisto_tar_file, abundance_h5_file, abundance_tsv_file, run_info_file):
        """
        Function to compress the Kallisto results into a tar file containing a
        single directory with the outputs from kallisto quant
        """
        output_file_pregz = kallisto_tar_file.replace('.tar.gz', '.tar')

        if os.path.isfile(kallisto_tar_file):
            os.remove(kallisto_tar_file)

        tar = tarfile.open(output_file_pregz, "w")
        tar.add(abundance_h5_file, arcname='kallisto/abundance.h5')
        tar.add(abundance_tsv_file, arcname='kallisto/abundance.tsv')
        tar.add(run_info_file, arcname='kallisto/run_info.json')
        tar.close()

        try:
            command_line = 'pigz ' + output_file_pregz
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()
        except OSError:
            logger.warn("OSERROR: pigz not installed, using gzip")
            command_line = 'gzip ' + output_file_pregz
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()

    def run(self, input_files, input_metadata, output_files):
        """
        Tool for calculating the level of expression

        Parameters
        ----------
        input_files : list
            Kallisto index file for the
            FASTQ file for the experiemtnal alignments
        input_metadata : list

        Returns
        -------
        array : list
            First element is a list of the index files. Second element is a
            list of the matching metadata
        """

        # input and output share most metadata
        output_metadata = {}

        # file_loc = input_files[1].split("/")
        # output_dir = "/".join(file_loc[0:-1])

        # abundance_h5_file = output_dir + "/abundance.h5"
        # abundance_tsv_file = output_dir + "/abundance.tsv"
        # run_info_file = output_dir + "/run_info.json"

        sources = [input_metadata["cdna"].file_path]

        if "fastq2" not in input_files:
            self.kallisto_quant_single(
                input_files["index"],
                input_files["fastq1"],
                output_files["kallisto_tar_file"],
                self.configuration["kallisto_bootstrap_param"]
            )
            sources.append(input_metadata["fastq1"].file_path)
            # results = compss_wait_on(results)
        elif "fastq2" in input_files:
            # handle error
            self.kallisto_quant_paired(
                input_files["index"],
                input_files["fastq1"],
                input_files["fastq2"],
                output_files["kallisto_tar_file"],
                self.configuration["kallisto_bootstrap_param"]
            )
            sources.append(input_metadata["fastq1"].file_path)
            sources.append(input_metadata["fastq2"].file_path)
            # results = compss_wait_on(results)
        else:
            return ({}, {})

        output_metadata = {
            "kallisto_tar_file": Metadata(
                data_type="data_rnaseq",
                file_type="TAR",
                file_path=output_files["kallisto_tar_file"],
                sources=sources,
                taxon_id=input_metadata["cdna"].taxon_id,
                meta_data={
                    "assembly": input_metadata["cdna"].meta_data["assembly"],
                    "tool": "kallisto_quant"
                }
            )
        }

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
