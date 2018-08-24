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
import os
import glob
import shutil
# from subprocess import CalledProcessError
from subprocess import PIPE
from subprocess import Popen

from basic_modules.tool import Tool
from utils import logger

from tool.common import format_utils

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
    # from utils.dummy_pycompss import compss_wait_on # pylint: disable=ungrouped-imports
    # from utils.dummy_pycompss import constraint # pylint: disable=ungrouped-imports


# ------------------------------------------------------------------------------

class tbSegmentTool(Tool):  # pylint: disable=invalid-name
    """
    Tool for finding tads and compartments in an adjacency matrix
    """

    def __init__(self):
        """
        Init function
        """
        logger.info("TADbit - Normalize")
        Tool.__init__(self)

    @task(bamin=FILE_IN, biases=FILE_IN, resolution=IN, workdir=IN,
          tad_dir=FILE_OUT, compartment_dir=FILE_OUT)
    def tb_segment(  # pylint: disable=too-many-locals,too-many-statements,unused-argument,no-self-use,too-many-arguments
            self, bamin, biases, resolution, callers, chromosomes,
            workdir, fasta=None, ncpus="1"):
        """
        Function to find tads and compartments in the Hi-C
        matrix

        Parameters
        ----------
        bamin : str
            Location of the tadbit bam paired reads
        biases : str
            Location of the pickle hic biases
        resolution : int
            Resolution of the Hi-C
        callers: str
            1 for ta calling, 2 for compartment calling
        workdir : str
            Location of working directory
        ncpus : int
            Number of cpus to use

        Returns
        -------
        compartments : str
            Location of tsv file with compartment definition
        tads : str
            Location of tsv file with tad definition
        filtered_bins : str
            Location of filtered_bins png

        """
        logger.info("TB SEGMENT: {0} {1} {2}".format(bamin, resolution, workdir))

        _cmd = [
            'tadbit', 'segment',
            '--nosql', '--mreads', bamin,
            '--workdir', workdir,
            '--resolution', resolution,
            '--cpu', str(ncpus),
            '--nosql'
            ]

        if '2' not in callers:
            _cmd.append('--only_tads')
        if '1' not in callers:
            _cmd.append('--only_compartments')

        if chromosomes:
            _cmd.append('--chromosomes')
            _cmd.append(chromosomes)
        if fasta:
            _cmd.append('--fasta')
            _cmd.append(fasta)

        if biases:
            _cmd.append('--biases')
            _cmd.append(biases)

        output_metadata = {}
        output_files = []

        out, err = Popen(_cmd, stdout=PIPE, stderr=PIPE).communicate()
        logger.info(out)
        logger.info(err)

        if '1' in callers:
            tad_dir = os.path.join(workdir, '06_segmentation',
                                   'tads_%s' % (format_utils.nice(int(resolution))))
            clean_headers(tad_dir)
            output_files.append(tad_dir)
        if '2' in callers:
            cmprt_dir = os.path.join(workdir, '06_segmentation',
                                     'compartments_%s' % (format_utils.nice(int(resolution))))
            clean_headers(cmprt_dir)
            output_files.append(cmprt_dir)

        return (output_files, output_metadata)

    def run(self, input_files, output_files, metadata=None):  # pylint: disable=too-many-locals
        """
        The main function to the predict TAD sites and compartments for a given resolution from
        the Hi-C matrix

        Parameters
        ----------
        input_files : list
            bamin : str
                Location of the tadbit bam paired reads
            biases : str
                Location of the pickle hic biases
        metadata : dict
            resolution : int
                Resolution of the Hi-C
            workdir : str
                Location of working directory
            ncpus : int
                Number of cpus to use



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
        if 'resolution' in metadata:
            resolution = metadata['resolution']

        ncpus = 1
        if 'ncpus' in metadata:
            ncpus = metadata['ncpus']
        biases = chromosomes = fasta = None
        if len(input_files) > 1:
            biases = input_files[1]
        if "chromosomes" in metadata:
            chromosomes = metadata['chromosomes']
        if "fasta" in metadata:
            fasta = metadata['fasta']
        callers = "1"
        if "callers" in metadata:
            callers = metadata['callers']

        root_name = os.path.dirname(os.path.abspath(bamin))
        if 'workdir' in metadata:
            root_name = metadata['workdir']

        # input and output share most metadata

        output_files, output_metadata = self.tb_segment(bamin, biases, resolution,
                                                        callers, chromosomes, root_name,
                                                        fasta, ncpus)

        return (output_files, output_metadata)


def clean_headers(fpath):
    """
        Replaces spaces by underscores in the headers of tsv files
    """
    os.chdir(fpath)

    for fl_files in glob.glob("*.tsv"):
        tsv_file = open(fl_files)
        line = tsv_file.readline()
        line = line.replace(' ', '_')
        dest_file = os.path.join(os.path.dirname(fl_files), 'vre_' + os.path.basename(fl_files))
        to_file = open(dest_file, mode="w")
        to_file.write(line)
        shutil.copyfileobj(tsv_file, to_file)
        os.unlink(fl_files)
