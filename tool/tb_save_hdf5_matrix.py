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

from pytadbit.parsers.hic_parser import load_hic_data_from_reads

# ------------------------------------------------------------------------------

class tbSaveAdjacencyHDF5Tool(Tool):
    """
    Tool for filtering out experimetnal artifacts from the aligned data
    """
    
    def __init__(self):
        """
        Init function
        """
        print "TADbit save adjacency matrix"
    
    @task(adj_list = FILE_IN, adj_matrix = FILE_INOUT, resolution = IN, assembly = IN)
    def tb_matrix_hdf5(self, adj_list, adj_hdf5, resolution, assembly):
        """
        Function to the Hi-C matrix into an HDF5 file
        
        Parameters
        ----------
        adj_list : str
                Location of the adjacency list
        hdf5_file : str
            Location of the HDF5 output matrix file
        resolution : int
            Resolution to read teh Hi-C adjacency list at
        assembly : str
            Assembly used to align the experimental data to
        
        Returns
        -------
        hdf5_file : str
            Location of the HDF5 output matrix file
        
        """
        hic_data = load_hic_data_from_reads(adj_list, resolution=resolution)

        dSize = len(hic_data)
        d = np.zeros([dSize, dSize], dtype='int32')
        d += hic_data.get_matrix()
        
        f = h5py.File(adj_hdf5, "a")
        dset = f.create_dataset(str(resolution), (dSize, dSize), dtype='int32', chunks=True, compression="gzip")
        dset[0:dSize,0:dSize] += d
        f.close()
        
        return True
    
    
    def run(self, input_files, metadata):
        """
        The main function save the adjacency list from Hi-C into an HDF5 index
        file at the defined resolutions.
        
        Parameters
        ----------
        input_files : list
            adj_list : str
                Location of the adjacency list
            hdf5_file : str
                Location of the HDF5 output matrix file
        metadata : dict
            resolutions : list
                Levels of resolution for the adjacency list to be daved at
            assembly : str
                Assembly of the aligned sequences

        
        
        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects
        
        """
        
        adj_list = input_files[0]
        hdf5_file = input_files[1]
        resolutions = [1000000]
        if 'resolutions' in metadata:
            resolutions = metadata['resolutions']

        # input and output share most metadata
        output_metadata = {}
        
        for resolution in resolutions:
            # handle error
            if not self.tb_matrix_hdf5(adj_list, hdf5_file, resolution):
                output_metadata.set_exception(
                    Exception(
                        "tb_matrix_hdf5: Could not process files {}, {}.".format(*input_files)))
            
        return ([hdf5_file], output_metadata)

# ------------------------------------------------------------------------------
