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

import sys
import glob
import os
from subprocess import PIPE, Popen
import subprocess

from basic_modules.tool import Tool

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    # from pycompss.api.api import compss_wait_on
    # from pycompss.api.constraint import constraint
except ImportError:
    logger.info("[Warning] Cannot import \"pycompss\" API packages.")
    logger.info("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task  # pylint: disable=ungrouped-imports
    # from utils.dummy_pycompss import compss_wait_on  # pylint: disable=ungrouped-imports
    # from utils.dummy_pycompss import constraint


# ------------------------------------------------------------------------------

class tbNormalizeTool(Tool):  # pylint: disable=invalid-name
    """
    Tool for normalizing an adjacency matrix
    """

    def __init__(self):
        """
        Init function
        """
        logger.info("TADbit - Normalize")
        Tool.__init__(self)

    @task(bamin=FILE_IN, normalization=IN, resolution=IN, min_perc=IN,
          max_perc=IN, workdir=IN, biases=FILE_OUT, interactions_plot=FILE_OUT,
          filtered_bins_plot=FILE_OUT)
    def tb_normalize(self, bamin, normalization, resolution, min_perc, # pylint: disable=too-many-locals,too-many-statements,unused-argument,no-self-use,too-many-arguments
                     max_perc, workdir, ncpus="1", min_count=None, fasta=None,
                     mappability=None, rest_enzyme=None):
        """
        Function to normalize to a given resolution the Hi-C
        matrix

        Parameters
        ----------
        bamin : str
            Location of the tadbit bam paired reads
        normalization: str
            normalization(s) to apply. Order matters. Choices: [Vanilla, oneD]
        resolution : str
            Resolution of the Hi-C
        min_perc : str
            lower percentile from which consider bins as good.
        max_perc : str
            upper percentile until which consider bins as good.
        workdir : str
            Location of working directory
        ncpus : str
            Number of cpus to use
        min_count : str
            minimum number of reads mapped to a bin (recommended value
            could be 2500). If set this option overrides the perc_zero
        fasta: str
            Location of the fasta file with genome sequence, to compute GC content and
            number of restriction sites per bin. Required for oneD normalization
        mappability: str
            Location of the file with mappability, required for oneD normalization
        rest_enzyme: str
            For oneD normalization. Name of the restriction enzyme used to do the Hi-C experiment

        Returns
        -------
        hic_biases : str
            Location of HiC biases pickle file
        interactions : str
            Location of interaction decay vs genomic distance pdf
        filtered_bins : str
            Location of filtered_bins png

        """
        # chr_hic_data = read_matrix(matrix_file, resolution=int(resolution))

        logger.info("TB NORMALIZATION: {0} {1} {2} {3} {4} {5}".format(
            bamin, normalization, resolution, min_perc, max_perc, workdir))

        _cmd = [
            'tadbit', 'normalize',
            '--bam', bamin,
            '--normalization', normalization,
            '--workdir', workdir,
            '--resolution', resolution,
            '--cpus', str(ncpus)
            ]

        if min_perc:
            _cmd.append('--min_perc')
            _cmd.append(min_perc)
        if max_perc:
            _cmd.append('--max_perc')
            _cmd.append(max_perc)
        if min_count:
            _cmd.append('--min_count')
            _cmd.append(min_count)

        if normalization == 'oneD':
            _cmd.append('--fasta')
            _cmd.append(fasta)
            _cmd.append('--mappability')
            _cmd.append(mappability)
            _cmd.append('--renz')
            _cmd.append(rest_enzyme)

        output_metadata = {}
        output_files = []

        try:
            _ = subprocess.check_output(_cmd, stderr=subprocess.STDOUT,
                                        cwd=workdir)
        except subprocess.CalledProcessError as subp_err:
            logger.info(subp_err.output)
            if not min_count:
                logger.info("cis/trans ratio failed, trying with min_count. Disabling plot.")
                _cmd.append('--min_count')
                _cmd.append('10')
                _cmd.append('--normalize_only')
                try:
                    _ = subprocess.check_output(_cmd, stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError as subp_err:
                    logger.fatal(subp_err.output)

        os.chdir(workdir+"/04_normalization")
        for fl_file in glob.glob("biases_*.pickle"):
            output_files.append(os.path.abspath(fl_file))
            break
        for fl_file in glob.glob("interactions*.png"):
            output_files.append(os.path.abspath(fl_file))
            break
        for fl_file in glob.glob("filtered_bins_*.png"):
            output_files.append(os.path.abspath(fl_file))
            break

        return (output_files, output_metadata)

    def run(self, input_files, input_metadata, output_files):  # pylint: disable=too-many-locals
        """
        The main function for the normalization of the Hi-C matrix to a given resolution

        Parameters
        ----------
        input_files : list
            bamin : str
                Location of the tadbit bam paired reads
        metadata : dict
            normalization: str
                normalization(s) to apply. Order matters. Choices: [Vanilla, oneD]
            resolution : str
                Resolution of the Hi-C
            min_perc : str
                lower percentile from which consider bins as good.
            max_perc : str
                upper percentile until which consider bins as good.
            workdir : str
                Location of working directory
            ncpus : str
                Number of cpus to use
            min_count : str
                minimum number of reads mapped to a bin (recommended value
                could be 2500). If set this option overrides the perc_zero
            fasta: str
                Location of the fasta file with genome sequence, to compute GC content and
                number of restriction sites per bin. Required for oneD normalization
            mappability: str
                Location of the file with mappability, required for oneD normalization
            rest_enzyme: str
                For oneD normalization.
                Name of the restriction enzyme used to do the Hi-C experiment

        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects

        """

        bamin = input_files[0]

        if not os.path.isfile(bamin.replace('.bam', '.bam.bai')):
            logger.info('Creating bam index')
            _cmd = ['samtools', 'index', bamin]
            out, err = Popen(_cmd, stdout=PIPE, stderr=PIPE).communicate()
            logger.info(out)
            logger.info(err)

        resolution = '1000000'
        if 'resolution' in input_metadata:
            resolution = input_metadata['resolution']

        normalization = 'Vanilla'
        if 'normalization' in input_metadata:
            normalization = input_metadata['normalization']

        min_perc = max_perc = min_count = fasta = mappability = rest_enzyme = None
        ncpus = 1
        if 'ncpus' in input_metadata:
            ncpus = input_metadata['ncpus']
        if 'min_perc' in input_metadata:
            min_perc = input_metadata['min_perc']
        if 'max_perc' in input_metadata:
            max_perc = input_metadata['max_perc']
        if 'min_count' in input_metadata:
            min_count = input_metadata['min_count']
        if 'fasta' in input_metadata:
            fasta = input_metadata['fasta']
        if 'mappability' in input_metadata:
            mappability = input_metadata['mappability']
        if 'rest_enzyme' in input_metadata:
            rest_enzyme = input_metadata['rest_enzyme']

        root_name = os.path.dirname(os.path.abspath(bamin))
        if 'workdir' in input_metadata:
            root_name = input_metadata['workdir']

        # input and output share most metadata

        output_files, output_metadata = self.tb_normalize(bamin, normalization,
                                                          resolution, min_perc,
                                                          max_perc, root_name, ncpus,
                                                          min_count, fasta, mappability,
                                                          rest_enzyme)

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
