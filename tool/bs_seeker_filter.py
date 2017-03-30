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

import os

try:
    from pycompss.api.parameter import FILE_IN, FILE_OUT
    from pycompss.api.task import task
except ImportError :
    print "[Warning] Cannot import \"pycompss\" API packages."
    print "          Using mock decorators."
    
    from dummy_pycompss import *

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

from common import common

pwd = os.environ.get('PWD')
pwd_split = pwd.split('/')

if pwd_split[-1] != 'docs':
    on_rtd = os.environ.get('READTHEDOCS') == 'True'
    if on_rtd == False:
        from FilterReads import *

# ------------------------------------------------------------------------------

class filterReadsTool(Tool):
    """
    Script from BS-Seeker2 for filtering FASTQ files to remove repeats
    """
    
    def __init__(self):
        """
        Init function
        """
        print "BS-Seeker FilterReads wrapper"

    @task(infile = FILE_IN, outfile = FILE_OUT)
    def bss_seeker_filter(self, infile, outfile):
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
        FilterReads(infile, outfile, True)
        return True


    def run(self, input_files, metadata):
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
        
        # handle error
        if not self.bss_seeker_filter(file_name, output_file):
            output_metadata.set_exception(
                Exception(
                    "bs_seeker_filter: Could not process files {}, {}.".format(*input_files)))
            output_file = None
        return ([output_file], [output_metadata])

# ------------------------------------------------------------------------------