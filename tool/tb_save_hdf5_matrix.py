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

#from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

import numpy as np
import h5py
import pytadbit

from pytadbit import load_hic_data_from_reads
from pytadbit import read_matrix
from pytadbit.parsers.genome_parser import parse_fasta

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
        Tool.__init__(self)


    def tb_matrix_hdf5(self, hic_data, adj_hdf5, resolution, chromosomes):
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
        hic_data : hic_data
            Hi-C data object
        hdf5_file : str
            Location of the HDF5 output matrix file
        resolution : int
            Resolution to read teh Hi-C adjacency list at
        chromosomes : list
            List of listsd of the chromosome names and their size in the order
            that they are presented for indexing

        Returns
        -------
        hdf5_file : str
            Location of the HDF5 output matrix file

        """
        #hic_data = read_matrix(adj_list, resolution=resolution)

        d_size = len(hic_data)
        d = np.zeros([d_size, d_size], dtype='int32')
        d += hic_data.get_matrix()

        hdf5_handle = h5py.File(adj_hdf5, "a")
        dset = hdf5_handle.create_dataset(
            str(resolution),
            (d_size, d_size),
            dtype='int32',
            chunks=True,
            compression="gzip"
        )
        dset.attrs['chromosomes'] = chromosomes
        dset[0:d_size,0:d_size] += d
        hdf5_handle.close()

        return True


    def run(self, input_files, output_files, metadata=None):
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
            normalized : bool
                Whether the dataset should be normalised before saving


        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects

        """

        adjlist_file = input_files[0]
        genome_file = input_files[1]
        hdf5_file    = adjlist_file.replace(".tsv", ".hdf5")

        genome_seq = parse_fasta(genome_file)
        chromosomes = []
        for chr_id in genome_seq:
            chromosomes.append([chr_id, len(genome_seq[chr_id])])

        assembly = metadata['assembly']
        resolutions = metadata['resolutions']

        normalized = False
        if 'normalized' in metadata:
            normalized = metadata['normalized']

        #root_name = adjlist_file.split("/")

        #matrix_files = []
        tad_files = {}
        output_metadata = {}

        for resolution in resolutions:
            hic_data = load_hic_data_from_reads(adjlist_file, resolution=int(resolution))
            tad_files[resolution] = {}

            if normalized is False:
                hic_data.normalize_hic(iterations=9, max_dev=0.1)

            #save_matrix_file = "/".join(root_name[0:-1]) + '/adjlist_map_' + str(resolution) + '.tsv'
            #matrix_files.append(save_matrix_file)
            #hic_data.write_matrix(save_matrix_file, normalized=normalized)

            results = self.tb_matrix_hdf5(hic_data, hdf5_file, resolution, chromosomes)

        return ([hdf5_file], [output_metadata])

# ------------------------------------------------------------------------------
