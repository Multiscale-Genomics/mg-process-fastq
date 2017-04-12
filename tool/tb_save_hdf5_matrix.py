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
from pytadbit.parsers.hic_parser import read_matrix

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
    
    
    def tb_matrix_hdf5(self, adj_list, adj_hdf5, resolution, chromosomes):
        """
        Function to the Hi-C matrix into an HDF5 file

        This has to be run sequentially as it is not possible for multiple
        streams to write to the same HDF5 file. This is a run once and leave
        operatation. There also needs to be a check that no other process is
        writing to the HDF5 file at the same time. This should be done at the 
        stage and unstaging level to prevent to file getting written to by
        multiple processes and generating conflicts.
        
        This needs to include attributes for the chromosomes for each resolution
        - See the mg-rest-adjacency hdf5_reader for further details about the
          requirement. This prevents the need for secondary storage details
          outside of the HDF5 file.
        
        Parameters
        ----------
        adj_list : str
                Location of the adjacency list
        hdf5_file : str
            Location of the HDF5 output matrix file
        resolution : int
            Resolution to read teh Hi-C adjacency list at
        chromosomes : list
            List of listsd of the chromosome names and their size in the order that they are presented for indexing
        
        Returns
        -------
        hdf5_file : str
            Location of the HDF5 output matrix file
        
        """
        hic_data = read_matrix(adjlist_file, resolution=resolution)

        dSize = len(hic_data)
        d = np.zeros([dSize, dSize], dtype='int32')
        d += hic_data.get_matrix()
        
        f = h5py.File(hdf5_file, "a")
        dset = f.create_dataset(str(self.resolution), (dSize, dSize), dtype='int32', chunks=True, compression="gzip")
        dset.attrs['chromosomes'] = chromosomes
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
        
        adjlist_file = input_files[0]
        hdf5_file    = input_files[1]
        
        assembly = metadata['assembly']
        resolution = metadata['resolution']
        chromosomes = metadata['chromosomes_meta']

        root_name = adj_list.split("/")

        matrix_files = []
        tad_files = {}

        for resolution in resolutions:
            hic_data = load_hic_data_from_reads(adj_list, resolution=int(resolution))
            tad_files[resolution] = {}

            save_matrix_file = "/".join(root_name[0:-1]) + '/adjlist_map_' + str(resolution) + '.tsv'
            matrix_files.append(save_matrix_file)
            hic_data.write_matrix(save_matrix_file, normalized=normalized)

            # handle error
            if not self.tb_matrix_hdf5(adjlist_file, hdf5_file, resolution, chromosomes):
                output_metadata.set_exception(
                    Exception(
                        "tb_matrix_hdf5: Could not process files {}, {}.".format(*input_files)))
        
        return ([hdf5_file], [output_metadata])

# ------------------------------------------------------------------------------
