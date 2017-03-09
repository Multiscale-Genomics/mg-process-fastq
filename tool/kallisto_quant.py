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
        print "Kallisto Quantification"
    
    @task(fastq_file_loc=FILE_IN, cdna_idx_file=FILE_IN, wig_file = FILE_OUT)
    def kallisto_quant_single(self, cdna_idx_file, fastq_file_loc):
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
        
        command_line = 'kallisto quant -i ' + cdna_idx_file + ' -o .'
        fq_stats = self.seq_read_stats(fastq[0])
        command_line += ' --single -l ' + str(fq_stats['mean']) + ' -s ' + str(fq_stats['std']) + ' ' + fastq_file_loc
        
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()
        
        return True
    
    
    @task(fastq_file_loc_01=FILE_IN, fastq_file_loc_02=FILE_IN, cdna_idx_file=FILE_IN, wig_file = FILE_OUT)
    def kallisto_quant_paired(self, cdna_idx_file, fastq_file_loc_01, fastq_file_loc_02):
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
        
        command_line = 'kallisto quant -i ' + cdna_idx_file + ' -o .'
        command_line += ' ' + fastq_file_loc_01 + ' ' + fastq_file_loc_02
        
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()
        
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
        with open(file_in, 'r') as f:
            forthlines = itertools.islice(f, 1, None, 4)
            for line in forthlines:
                line = line.rstrip()
                total_len.append(len(line))

        s = std(total_len)
        m = mean(total_len)

        return {'mean' : m, 'std' : s}
    
    
    def run(self, input_files, metadata):
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
        output_metadata = dict(
            data_type=metadata[0]["data_type"],
            file_type=metadata[0]["file_type"],
            meta_data=metadata[0]["meta_data"])
        
        results_wig_file = input_files[0].replace('.fastq', '.wig')
        
        if len(input_files) == 2:
            # handle error
            if not self.kallisto_quant_single(input_files[0], input_files[1], results_wig_file):
                output_metadata.set_exception(
                    Exception(
                        "kallisto_quant_single: Could not process files {}, {}.".format(*input_files)))
            results_wig_file = None
        elif len(input_files) == 3:
            # handle error
            if not self.kallisto_quant_paired(input_files[0], input_files[1], input_files[2], results_wig_file):
                output_metadata.set_exception(
                    Exception(
                        "kallisto_quant_single: Could not process files {}, {}.".format(*input_files)))
            results_wig_file = None
        else:
            output_metadata.set_exception(
                Exception(
                    "Error: Incorrect number of parameters {}, {}.".format(*input_files)))
            results_wig_file = None
        
        return ([results_wig_file], [output_metadata])

# ------------------------------------------------------------------------------
