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

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    print ("[Warning] Cannot import \"pycompss\" API packages.")
    print ("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT, IN
    from dummy_pycompss import task
    from dummy_pycompss import compss_wait_on

from basic_modules.tool import Tool

# ------------------------------------------------------------------------------

class kallistoQuantificationTool(Tool):
    """
    Tool for quantifying RNA-seq alignments to calculate expression levels of
    genes within a genome.
    """

    def __init__(self):
        """
        Init function
        """
        print("Kallisto Quantification")
        Tool.__init__(self)

    @task(
        cdna_idx_file=FILE_IN,
        fastq_file_loc=FILE_IN,
        fq_stats=IN,
        abundance_h5_file=FILE_OUT,
        abundance_tsv_file=FILE_OUT,
        run_info_file=FILE_OUT)
    def kallisto_quant_single(
            self, cdna_idx_file, fastq_file_loc, fq_stats,
            abundance_h5_file, abundance_tsv_file, run_info_file):
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

        output_dir = fastq_file_loc.split('/')

        command_line = 'kallisto quant -i ' + cdna_idx_file + ' '
        command_line += ' -o ' + '/'.join(output_dir[0:-1]) + "/"
        command_line += ' --single -l ' + str(fq_stats['mean']) + ' '
        command_line += '-s ' + str(fq_stats['std']) + ' ' + fastq_file_loc

        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()

        return True

    @task(
        fastq_file_loc_01=FILE_IN,
        fastq_file_loc_02=FILE_IN,
        cdna_idx_file=FILE_IN,
        abundance_h5_file=FILE_OUT,
        abundance_tsv_file=FILE_OUT,
        run_info_file=FILE_OUT)
    def kallisto_quant_paired(
            self, cdna_idx_file, fastq_file_loc_01, fastq_file_loc_02,
            abundance_h5_file, abundance_tsv_file, run_info_file):
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
        command_line += '-o ' + '/'.join(output_dir[0:-1]) + "/ "
        command_line += fastq_file_loc_01 + ' ' + fastq_file_loc_02

        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        return True

    def seq_read_stats(self, file_in):
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

        return {'mean' : length_mean, 'std' : length_sd}

    def run(self, input_files, output_files, metadata=None):
        """
        Tool for calculating the level of expression

        Parameters
        ----------
        input_files : list
            Kallisto index file for the
            FASTQ file for the experiemtnal alignments
        metadata : list

        Returns
        -------
        array : list
            First element is a list of the index files. Second element is a
            list of the matching metadata
        """

        # input and output share most metadata
        output_metadata = {}

        file_loc = input_files[1].split("/")
        output_dir = "/".join(file_loc[0:-1])

        abundance_h5_file = output_dir + "/abundance.h5"
        abundance_tsv_file = output_dir + "/abundance.tsv"
        run_info_file = output_dir + "/run_info.json"

        if len(input_files) == 2:
            fq_stats = self.seq_read_stats(input_files[1])
            results = self.kallisto_quant_single(
                input_files[0], input_files[1], fq_stats, abundance_h5_file,
                abundance_tsv_file, run_info_file)
            results = compss_wait_on(results)
        elif len(input_files) == 3:
            # handle error
            results = self.kallisto_quant_paired(
                input_files[0], input_files[1], input_files[2],
                abundance_h5_file, abundance_tsv_file, run_info_file,)
            results = compss_wait_on(results)
        else:
            abundance_h5_file = None
            abundance_tsv_file = None
            run_info_file = None

        return (
            [abundance_h5_file, abundance_tsv_file, run_info_file],
            [output_metadata])

# ------------------------------------------------------------------------------
