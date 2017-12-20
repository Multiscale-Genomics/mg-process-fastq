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

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, FILE_INOUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
    # from pycompss.api.constraint import constraint
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT, FILE_INOUT, IN
    from dummy_pycompss import task
    from dummy_pycompss import compss_wait_on
    #from dummy_pycompss import constraint

from basic_modules.tool import Tool

# ------------------------------------------------------------------------------

class tbNormalizeTool(Tool):
    """
    Tool for normalizing an adjacency matrix
    """

    def __init__(self):
        """
        Init function
        """
        print("TADbit - Normalize")
        Tool.__init__(self)

    @task(bamin=FILE_IN, resolution=IN, min_perc=IN, max_perc=IN, workdir=IN)
    # @constraint(ProcessorCoreCount=16)
    def tb_normalize(self, bamin, resolution, min_perc, max_perc, workdir, ncpus="1", min_count=None):
        """
        Function to normalize to a given resolution the Hi-C
        matrix

        Parameters
        ----------
        bamin : str
            Location of the tadbit bam paired reads
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

        Returns
        -------
        hic_biases : str
            Location of HiC biases pickle file
        interactions : str
            Location of interaction decay vs genomic distance pdf
        filtered_bins : str
            Location of filtered_bins png

        """
        #chr_hic_data = read_matrix(matrix_file, resolution=int(resolution))

        print("TB NORMALIZATION:", bamin, resolution, min_perc, max_perc, workdir)

        _cmd = [
            'tadbit', 'normalize',
            '--bam', bamin,
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

        output_metadata = {}
        output_files = []

        try:
            proc = subprocess.check_output(_cmd, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError:
            if not min_count:
                print("cis/trans ratio failed, trying with min_count. Disabling plot.")
                _cmd.append('--min_count')
                _cmd.append('10')
                _cmd.append('--normalize_only')
                proc = subprocess.check_output(_cmd, stderr=subprocess.STDOUT)

        print(proc)
        os.chdir(workdir+"/04_normalization")
        for fl in glob.glob("biases_*.pickle"):
            output_files.append(os.path.abspath(fl))
            break
        for fl in glob.glob("interactions*.png"):
            output_files.append(os.path.abspath(fl))
            break
        for fl in glob.glob("filtered_bins_*.png"):
            output_files.append(os.path.abspath(fl))
            break

        return (output_files, output_metadata)

    def run(self, input_files, output_files, metadata=None):
        """
        The main function for the normalization of the Hi-C matrix to a given resolution

        Parameters
        ----------
        input_files : list
            bamin : str
                Location of the tadbit bam paired reads
        metadata : dict
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



        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects

        """

        bamin = input_files[0]

        if not os.path.isfile(bamin.replace('.bam', '.bam.bai')):
            print('Creating bam index')
            _cmd = ['samtools', 'index', bamin]
            out, err = Popen(_cmd, stdout=PIPE, stderr=PIPE).communicate()
            print(out)
            print(err)

        resolution = '1000000'
        if 'resolution' in metadata:
            resolution = metadata['resolution']

        min_perc = max_perc = min_count = None
        ncpus = 1
        if 'ncpus' in metadata:
            ncpus = metadata['ncpus']
        if 'min_perc' in metadata:
            min_perc = metadata['min_perc']
        if 'max_perc' in metadata:
            max_perc = metadata['max_perc']
        if 'min_count' in metadata:
            min_count = metadata['min_count']

        root_name = os.path.dirname(os.path.abspath(bamin))
        if 'workdir' in metadata:
            root_name = metadata['workdir']

        # input and output share most metadata

        output_files, output_metadata = self.tb_normalize(bamin, resolution, min_perc, max_perc, root_name, ncpus, min_count)

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
