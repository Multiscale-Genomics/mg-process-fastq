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

from pytadbit.mapping.mapper import full_mapping

# ------------------------------------------------------------------------------

class tbFullMappingTool(Tool):
    """
    Tool for mapping fastq paired end files to the GEM index files
    """
    
    def __init__(self):
        """
        Init function
        """
        print "TADbit full_mapping"
    
    @task(gem_file = FILE_IN, fastq_file = FILE_IN, windows = IN, window1 = FILE_INOUT, window2 = FILE_INOUT, window3 = FILE_INOUT, window4 = FILE_INOUT)
    @constraint(ProcessorCoreCount=8)
    def tb_full_mapping(self, gem_file, fastq_file, windows, window1, window2, window3, window4):
        """
        Function to map the FASTQ files to the GEM file over different window
        sizes ready for alignment
        
        Parameters
        ----------
        gem_file : str
            Location of the genome GEM index file
        fastq_file_bgd : str
            Location of the FASTQ file
        windows : list
            List of lists with the window sizes to be computed
        window1 : str
            Location of the first window index file
        window2 : str
            Location of the second window index file
        window3 : str
            Location of the third window index file
        window4 : str
            Location of the fourth window index file
        
        
        Returns
        -------
        window1 : str
            Location of the first window index file
        window2 : str
            Location of the second window index file
        window3 : str
            Location of the third window index file
        window4 : str
            Location of the fourth window index file
        
        """
        od = fastq_file.split("/")
        output_dir = "/".join(od[0:-1])

        map_files = full_mapping(gem_file, fastq_file, output_dir, windows=windows, frag_map=False, nthreads=8, clean=True, temp_dir='/tmp/')
        
        return True
    
    
    def run(self, input_files, metadata):
        """
        The main function to map the FASTQ files to the GEM file over different
        window sizes ready for alignment
        
        Parameters
        ----------
        input_files : list
            gem_file : str
                Location of the genome GEM index file
            fastq_file_bgd : str
                Location of the FASTQ file
        metadata : dict
            windows : list
                List of lists with the window sizes to be computed
        
        
        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects
        
        """
        
        gem_file = input_files[0]
        fastq_file = input_files[1]
        windows = meta_data['windows']
        
        root_name = fastq_file.split("/")
        root_name[-1].replace('.fa', '')
        
        name = root_name[-1]
        
        window1 = '/'.join(root_name) + "_full_1-25.map"
        window2 = '/'.join(root_name) + "_full_1-50.map"
        window3 = '/'.join(root_name) + "_full_1-75.map"
        window4 = '/'.join(root_name) + "_full_1-100.map"
        
        # input and output share most metadata
        output_metadata = {}
        
        # handle error
        if not self.tb_full_mapping(gem_file, fastq_file, windows, window1, window2, window3, window4):
            output_metadata.set_exception(
                Exception(
                    "tb_full_mapping: Could not process files {}, {}.".format(*input_files)))
            window1 = None
            window2 = None
            window3 = None
            window4 = None
        return ([window1, window2, window3, window4], output_metadata)

# ------------------------------------------------------------------------------
