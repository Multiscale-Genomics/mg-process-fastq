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

try:
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT, IN
    from dummy_pycompss import task

from basic_modules.tool import Tool

# ------------------------------------------------------------------------------

class filterReadsTool(Tool):
    """
    Script from BS-Seeker2 for filtering FASTQ files to remove repeats
    """

    def __init__(self):
        """
        Init function
        """
        print("BS-Seeker FilterReads wrapper")
        Tool.__init__(self)

    @task(infile=FILE_IN, outfile=FILE_OUT, bss_path=IN)
    def bss_seeker_filter(self, infile, outfile, bss_path):
        """
        This is optional, but removes reads that can be problematic for the
        alignment of whole genome datasets.

        If performing RRBS then this step can be skipped

        This is a function that is installed as part of the BS-Seeker
        installation process.

        Parameters
        ----------
        infile : str
            Location of the FASTQ file

        Returns
        -------
        outfile : str
            Location of the filtered FASTQ file
        """
        command_line = (
            "python " + bss_path + "/FilterReads.py"
            " -i " + infile + ""
            " -o " + outfile + ""
        ).format()

        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        return True


    def run(self, input_files, output_files, metadata=None):
        """
        Tool for filtering duplicate entries from FASTQ files using BS-Seeker2

        Parameters
        ----------
        input_files : list
            FASTQ file
        metadata : list

        Returns
        -------
        array : list
            Location of the filtered FASTQ file
        """

        file_name = input_files[0]
        output_file = file_name + '.filtered.fastq'

        output_metadata = {}

        results = self.bss_seeker_filter(file_name, output_file, metadata["bss_path"])

        if results is False:
            pass

        return ([output_file], output_metadata)

# ------------------------------------------------------------------------------
