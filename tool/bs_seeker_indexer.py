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

import os, shutil, shlex, subprocess

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

#if pwd_split[-1] != 'docs':
#    on_rtd = os.environ.get('READTHEDOCS') == 'True'
#    if on_rtd == False:
#        from bs_index.wg_build import *

# ------------------------------------------------------------------------------

class bssIndexerTool(Tool):
    """
    Script from BS-Seeker2 for building the index for alignment. In this case
    it uses Bowtie2.
    """
    
    def __init__(self):
        """
        Init function
        """
        print "BS-Seeker Indexer wrapper"

    @task(fasta_file = FILE_IN, aligner = IN, aligner_path = IN, bss_path = IN, ref_path = IN, bam_out = FILE_INOUT)
    def bss_build_index(self, fasta_file, aligner, aligner_path, bss_path, ref_path, bam_out):
        """
        Function to submit the FASTA file for the reference sequence and build
        the required index file used by the aligner.
        
        Parameters
        ----------
        fasta_file : str
            Location of the genome FASTA file
        aligner : str
            Aligner to use by BS-Seeker2. Currently only bowtie2 is available in
            this build
        aligner_path : str
            Location of the aligners binary file
        ref_path : str
            Location of the indexes for the FASTA file
        bam_out : str
            Location of the output bam alignment file

        Returns
        -------
        bam_out : str
            Location of the output bam alignment file
        """
        command_line = ("python " + bss_path + "/bs_seeker2-build.py"
            " -f " + fasta_file + ""
            " --aligner " + aligner + " --path " + aligner_path + ""
            " --db " + ref_path).format()
        
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()

        return True


    def run(self, input_files, metadata):
        """
        Tool for indexing the genome assembly using BS-Seeker2. In this case it
        is using Bowtie2
        
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
        gd = file_name.split("/")
        genome_dir = '/' + '/'.join(gd[:-1])
        
        aligner      = metadata['aligner']
        aligner_path = metadata['aligner_path']
        bss_path     = metadata['bss_path']


        output_file = file_name + '.filtered.bam'
        
        # input and output share most metadata
        output_metadata = {}
        
        # handle error
        if not self.bss_build_index(file_name, aligner, aligner_path, bss_path, genome_dir, output_file):
            output_metadata.set_exception(
                Exception(
                    "bs_seeker_filter: Could not process files {}, {}.".format(*input_files)))
            output_file = None
        return ([output_file], [output_metadata])

# ------------------------------------------------------------------------------