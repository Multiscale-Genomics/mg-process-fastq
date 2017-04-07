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

try :
    from pycompss.api.parameter import FILE_IN, FILE_OUT, FILE_INOUT
    from pycompss.api.task import task
except ImportError :
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")
    
    from dummy_pycompss import *

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

import numpy as np
import h5py
import pytadbit

from pytadbit.parsers.genome_parser import parse_fasta
from pytadbit.parsers.map_parser import parse_map
from pytadbit.mapping import get_intersection

# ------------------------------------------------------------------------------

class tbParseMappingTool(Tool):
    """
    Tool for parsing the mapped reads and generating the list of paired ends
    that have a match at both ends.
    """
    
    def __init__(self):
        """
        Init function
        """
        print "TADbit parse mapping"
    
    @task(genome_seq = IN, enzyme_name = IN,
        window1_1 = FILE_IN, window1_2 = FILE_IN, window1_3 = FILE_IN, window1_4 = FILE_IN,
        window2_1 = FILE_IN, window2_2 = FILE_IN, window2_3 = FILE_IN, window2_4 = FILE_IN,
        reads = FILE_INOUT)
    @constraint(ProcessorCoreCount=16)
    def tb_parse_mapping_iter(self,
        genome_seq, enzyme_name,
        window1_1, window1_2, window1_3, window1_4,
        window2_1, window2_2, window2_3, window2_4,
        reads):
        """
        Function to map the aligned reads and return the matching pairs
        
        Parameters
        ----------
        genome_seq : dict
            Object containing the sequence of each of the chromosomes
        enzyme_name : str
            Name of the enzyme used to digest the genome
        window1_1 : str
            Location of the first window index file
        window1_2 : str
            Location of the second window index file
        window1_3 : str
            Location of the third window index file
        window1_4 : str
            Location of the fourth window index file
        window2_1 : str
            Location of the first window index file
        window2_2 : str
            Location of the second window index file
        window2_3 : str
            Location of the third window index file
        window2_4 : str
            Location of the fourth window index file
        reads : str
            Location of the reads thats that has a matching location at both
            ends of the paired reads
        
        
        Returns
        -------
        reads : str
            Location of the intersection of mapped reads that have matching
            reads in both pair end files
        
        """

        root_name = reads.split("/")
        
        reads1 = "/".join(root_name[0:-1]) + '/reads_1.tsv'
        reads2 = "/".join(root_name[0:-1]) + '/reads_2.tsv'
        
        parse_map(
            [window1_1, window1_2, window1_3, window1_4],
            [window2_1, window2_2, window2_3, window2_4],
            out_file1=reads1,
            out_file2=reads2,
            genome_seq=genome_seq,
            re_name=enzyme_name,
            verbose=True,
            ncpus=16
        )

        intersections = get_intersection(reads1, reads2, reads, verbose=True)
        
        return True


    @task(genome_seq = IN, enzyme_name = IN,
        window1_full = FILE_IN, window1_frag = FILE_IN,
        window2_full = FILE_IN, window2_frag = FILE_IN,
        reads = FILE_INOUT)
    @constraint(ProcessorCoreCount=16)
    def tb_parse_mapping_frag(self,
        genome_seq, enzyme_name,
        window1_full, window1_frag, window2_full, window2_frag,
        reads):
        """
        Function to map the aligned reads and return the matching pairs
        
        Parameters
        ----------
        genome_seq : dict
            Object containing the sequence of each of the chromosomes
        enzyme_name : str
            Name of the enzyme used to digest the genome
        window1_full : str
            Location of the first window index file
        window1_full : str
            Location of the first window index file
        window2_full : str
            Location of the first window index file
        window2_frag : str
            Location of the second window index file
        reads : str
            Location of the reads thats that has a matching location at both
            ends of the paired reads
        
        
        Returns
        -------
        reads : str
            Location of the intersection of mapped reads that have matching
            reads in both pair end files
        
        """

        root_name = reads.split("/")
        
        reads1 = "/".join(root_name[0:-1]) + '/reads_1.tsv'
        reads2 = "/".join(root_name[0:-1]) + '/reads_2.tsv'
        
        parse_map(
            [window1_full, window1_frag],
            [window2_full, window2_frag],
            out_file1=reads1,
            out_file2=reads2,
            genome_seq=genome_seq,
            re_name=enzyme_name,
            verbose=True,
            ncpus=16
        )

        intersections = get_intersection(reads1, reads2, reads, verbose=True)
        
        return True
    
    
    def run(self, input_files, metadata):
        """
        The main function to map the aligned reads and return the matching pairs
        
        Parameters
        ----------
        input_files : list
            genome_file : str
                Location of the genome FASTA file
            window1_1 : str
                Location of the first window index file
            window1_2 : str
                Location of the second window index file
            window1_3 : str
                Location of the third window index file
            window1_4 : str
                Location of the fourth window index file
            window2_1 : str
                Location of the first window index file
            window2_2 : str
                Location of the second window index file
            window2_3 : str
                Location of the third window index file
            window2_4 : str
                Location of the fourth window index file
        metadata : dict
            windows : list
                List of lists with the window sizes to be computed
            enzyme_namer : str
                Restricture enzyme name
            mapping : list
                The mapping function used. The options are iter or frag.
        
        
        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects
        
        """
        
        genome_file = input_files[0]
        
        enzyme_name = metadata['enzyme_name']
        mapping_list = metadata['mapping']

        root_name = fastq_file.split("/")
        reads  = "/".join(root_name[0:-1]) + '/both_map.tsv'
        
        genome_seq = parse_fasta(genome_file)
        
        # input and output share most metadata
        output_metadata = {}
        
        if mapping_list[0] == mapping_list[1]:
            if mapping_list[0] == 'iter':
                window1_1 = input_files[1]
                window1_2 = input_files[2]
                window1_3 = input_files[3]
                window1_4 = input_files[4]

                window2_1 = input_files[5]
                window2_2 = input_files[6]
                window2_3 = input_files[7]
                window2_4 = input_files[8]

                # handle error
                if not self.tb_parse_mapping_iter(genome_seq, enzyme_name,
                    window1_1, window1_2, window1_3, window1_4,
                    window2_1, window2_2, window2_3, window2_4,
                    reads):

                    output_metadata.set_exception(
                        Exception(
                            "tb_parse_mapping_iter: Could not process files {}, {}.".format(*input_files)))
                    reads = None
                return ([reads], output_metadata)
            elif mapping_list[0] == 'frag':
                window1_1 = input_files[1]
                window1_2 = input_files[2]
                
                window2_1 = input_files[3]
                window2_2 = input_files[4]
                
                # handle error
                if not self.tb_parse_mapping_iter(genome_seq, enzyme_name,
                    window1_full, window1_frag,
                    window2_full, window2_frag,
                    reads):

                    output_metadata.set_exception(
                        Exception(
                            "tb_parse_mapping_iter: Could not process files {}, {}.".format(*input_files)))
                    reads = None
                return ([reads], output_metadata)
            else:
                output_metadata.set_exception(
                    Exception(
                        "tb_parse_mapping_*: Unspecified mapping version {}, {}.".format(*metadata)))
                reads = None
                return ([reads], output_metadata)

# ------------------------------------------------------------------------------
