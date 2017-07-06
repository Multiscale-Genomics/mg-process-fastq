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
        Tool.__init__(self)
    
    @task(gem_file = FILE_IN, fastq_file = FILE_IN, windows = IN, window1 = FILE_INOUT, window2 = FILE_INOUT, window3 = FILE_INOUT, window4 = FILE_INOUT)
    @constraint(ProcessorCoreCount=32)
    def tb_full_mapping_iter(self, gem_file, fastq_file, windows, window1, window2, window3, window4):
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

        map_files = full_mapping(gem_file, fastq_file, output_dir, windows=windows, frag_map=False, nthreads=32, clean=True, temp_dir='/tmp/')
        
        return True
    

    @task(gem_file = FILE_IN, fastq_file = FILE_IN, enzyme_name = IN, windows = IN, full_file = FILE_INOUT, frag_file = FILE_INOUT)
    @constraint(ProcessorCoreCount=16)
    def tb_full_mapping_frag(self, gem_file, fastq_file, enzyme_name, windows, full_file, frag_file):
        """
        Function to map the FASTQ files to the GEM file based on fragments
        derived from the restriction enzyme that was used.
        
        Parameters
        ----------
        gem_file : str
            Location of the genome GEM index file
        fastq_file_bgd : str
            Location of the FASTQ file
        enzyme_name : str
            Restriction enzyme name (MboI)
        windows : list
            List of lists with the window sizes to be computed
        window_file : str
            Location of the first window index file
        
        
        Returns
        -------
        window_file : str
            Location of the window index file
        
        """
        od = fastq_file.split("/")
        output_dir = "/".join(od[0:-1])

        map_files = full_mapping(gem_file, fastq_file, output_dir, r_enz=enzyme_name, windows=windows, frag_map=True, nthreads=32, clean=True, temp_dir='/tmp/')
        
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
            enzyme_name : str
                Restriction enzyme used [OPTIONAL]
        
        
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
        root_name[-1] = root_name[-1].replace('.fa', '')
        
        name = root_name[-1]
        
        # input and output share most metadata
        output_metadata = {}
        
        if 'enzyme_name' in meta_data:
            full_file = '/'.join(root_name) + "_full_" + windows[0][0] + "-" + windows[0][1] + ".map"
            frag_file = '/'.join(root_name) + "_frag_" + windows[0][0] + "-" + windows[0][1] + ".map"

            # handle error
            if not self.tb_full_mapping_frag(gem_file, fastq_file, meta_data['enzyme_name'], windows, full_file, frag_file):
                output_metadata.set_exception(
                    Exception(
                        "tb_full_mapping_frag: Could not process files {}, {}.".format(*input_files)))
                full_file = None
                frag_file = None
            output_metadata['func'] = 'frag'
            return ([full_file, frag_file], output_metadata)
        else:
            window1 = '/'.join(root_name) + "_full_" + windows[0][0] + "-" + windows[0][1] + ".map"
            window2 = '/'.join(root_name) + "_full_" + windows[1][0] + "-" + windows[1][1] + ".map"
            window3 = '/'.join(root_name) + "_full_" + windows[2][0] + "-" + windows[2][1] + ".map"
            window4 = '/'.join(root_name) + "_full_" + windows[3][0] + "-" + windows[3][1] + ".map"

            # handle error
            if not self.tb_full_mapping_iter(gem_file, fastq_file, windows, window1, window2, window3, window4):
                output_metadata.set_exception(
                    Exception(
                        "tb_full_mapping_iter: Could not process files {}, {}.".format(*input_files)))
                window1 = None
                window2 = None
                window3 = None
                window4 = None
            output_metadata['func'] = 'iter'
            return ([window1, window2, window3, window4], output_metadata)

# ------------------------------------------------------------------------------
