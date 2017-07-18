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
import pysam

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT, IN
    from dummy_pycompss import task

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
        bss_path=IN, bam_file=FILE_IN, db_dir=IN,
        wig_file=FILE_OUT, cgmap_file=FILE_OUT, atcgmap_file=FILE_OUT)
    def bss_methylation_caller(self, bss_path, bam_file, db_dir, wig_file, cgmap_file, atcgmap_file):
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
        db_dir : str
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

        command_line = (
            "python " + bss_path + "/bs_seeker2-call_methylation.py "
            "--sorted --input " + str(bam_file) + " --wig " + str(wig_file) + " "
            "--CGmap " + str(cgmap_file) + " --ATCGmap " + str(atcgmap_file) + " "
            "--db " + db_dir).format()
        print ("command for methyl caller :", command_line)
        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        return True


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
        genome_dir = metadata['index_path']

        bss_path = metadata['bss_path']
        wig_file = file_name.replace('.bam', '.wig')
        cgmap_file = file_name.replace('.bam', '.cgmap.tsv')
        atcgmap_file = file_name.replace('.bam', '.atcgmap.tsv')

        # input and output share most metadata
        output_metadata = {}

        bam_file_handle = pysam.AlignmentFile(file_name, "rb")
        if ('SO' not in bam_file_handle.header['HD'] or
                bam_file_handle.header['HD']['SO'] == 'unsorted'):

            bam_file_handle.close()

            output_metadata['error'] = "bss_methylation_caller: Could not process files {}, {}.".format(*input_files)
            wig_file = None
            cgmap_file = None
            atcgmap_file = None

            return ([wig_file, cgmap_file, atcgmap_file], [output_metadata])

        self.bss_methylation_caller(
            bss_path, file_name, genome_dir, wig_file, cgmap_file, atcgmap_file)
        return ([wig_file, cgmap_file, atcgmap_file], [output_metadata])

# ------------------------------------------------------------------------------
