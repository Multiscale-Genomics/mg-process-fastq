"""
.. Copyright 2017 EMBL-European Bioinformatics Institute

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
import tarfile

import pysam

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT, IN
    from dummy_pycompss import task
    from dummy_pycompss import compss_wait_on

from basic_modules.tool import Tool


# ------------------------------------------------------------------------------

class bssMethylationCallerTool(Tool):
    """
    Script from BS-Seeker2 for methylation calling
    """

    def __init__(self):
        """
        Init function
        """
        print("BS-Seeker Methylation Caller")
        Tool.__init__(self)

    @task(
        bss_path=IN, bam_file=FILE_IN, genome_idx=FILE_IN,
        wig_file=FILE_OUT, cgmap_file=FILE_OUT, atcgmap_file=FILE_OUT)
    def bss_methylation_caller(
            self, bss_path, bam_file, genome_idx,
            wig_file, cgmap_file, atcgmap_file):
        """
        Takes the merged and sorted bam file and calls the methylation sites.
        Generates a wig file of the potential sites.
        This is performed by running the external program rather than
        reimplementing the code from the main function to make it easier when
        it comes to updating the changes in BS-Seeker2

        Parameters
        ----------
        bss_path : str
            Location of the Methylation caller script
        bam_file : str
            Location of the sorted bam alignment file
        genome_idx : str
            Location of the FASTA file
        wig_file : str
            Location of the wig results file
        cgmap_file : str
            Location of the CGmap results file
        atcgmap_file : str
            Location of the ATCGmap results file

        A full description of the

        Returns
        -------
        The wig, CTmap and ATCGmap files are returned to the matching locations.
        This is managed by pyCOMPS
        """

        g_dir = genome_idx.split("/")
        g_dir = "/".join(g_dir[:-1])

        tar = tarfile.open(genome_idx)
        tar.extractall(path=g_dir)
        tar.close()

        command_line = (
            "python " + bss_path + "/bs_seeker2-call_methylation.py "
            "--sorted --input " + str(bam_file) + " --wig " + str(wig_file) + " "
            "--CGmap " + str(cgmap_file) + " --ATCGmap " + str(atcgmap_file) + " "
            "--db " + genome_idx.replace('.tar.gz', '')).format()
        print ("command for methyl caller :", command_line)
        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        return True


    @task(bam_file=FILE_IN, returns=int)
    def check_header(self, bam_file):
        """
        Wrapper for the pysam SAMtools merge function

        Parameters
        ----------
        bam_file_1 : str
            Location of the bam file to merge into
        bam_file_2 : str
            Location of the bam file that is to get merged into bam_file_1
        """
        output = True

        bam_file_handle = pysam.AlignmentFile(bam_file, "rb")
        if ('SO' not in bam_file_handle.header['HD'] or
            bam_file_handle.header['HD']['SO'] == 'unsorted'):
            output = False
        bam_file_handle.close()

        return output


    def run(self, input_files, output_files, metadata=None):
        """
        Tool for methylation calling using BS-Seeker2.

        Parameters
        ----------
        input_files : list
            Sorted BAM file with the sequence alignments
        metadata : list

        Returns
        -------
        array : list
            Location of the output wig file
        """


        file_name = input_files[0]
        #gd = file_name.split("/")
        genome_idx = input_files[1]

        bss_path = metadata['bss_path']
        wig_file = file_name.replace('.bam', '.wig')
        cgmap_file = file_name.replace('.bam', '.cgmap.tsv')
        atcgmap_file = file_name.replace('.bam', '.atcgmap.tsv')

        # input and output share most metadata
        output_metadata = {}

        results = self.check_header(file_name)
        results = compss_wait_on(results)
        if results is False:
            output_metadata['error'] = "bss_methylation_caller: Could not process files {}, {}.".format(*input_files)
            wig_file = None
            cgmap_file = None
            atcgmap_file = None

            return ([wig_file, cgmap_file, atcgmap_file], [output_metadata])

        results = self.bss_methylation_caller(
            bss_path, file_name, genome_idx, wig_file, cgmap_file, atcgmap_file)
        results = compss_wait_on(results)

        return ([wig_file, cgmap_file, atcgmap_file], [output_metadata])

# ------------------------------------------------------------------------------
