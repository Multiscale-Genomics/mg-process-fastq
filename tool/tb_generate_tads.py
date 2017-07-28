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

import sys

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    from pycompss.api.constraint import constraint
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT, IN
    from dummy_pycompss import task
    from dummy_pycompss import constraint

#from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

#import numpy as np
#import h5py
#import pytadbit

from pytadbit import Chromosome
from pytadbit import read_matrix
from pytadbit import load_hic_data_from_reads

# ------------------------------------------------------------------------------

class tbGenerateTADsTool(Tool):
    """
    Tool for taking the adjacency lists and predicting TADs
    """

    def __init__(self):
        """
        Init function
        """
        print("TADbit - Generate TADs")
        Tool.__init__(self)

    @task(matrix_file = FILE_IN, resolution = IN, tad_file = FILE_OUT)
    @constraint(ProcessorCoreCount=16)
    def tb_generate_tads(self, expt_name, matrix_file, resolution, tad_file):
        """
        Function to the predict TAD sites for a given resolution from the Hi-C
        matrix

        Parameters
        ----------
        expt_name : str
                Location of the adjacency list
        matrix_file : str
            Location of the HDF5 output matrix file
        resolution : int
            Resolution to read the Hi-C adjacency list at
        tad_file : str
            Location of the output TAD file

        Returns
        -------
        tad_file : str
            Location of the output TAD file

        """
        chr_hic_data = read_matrix(matrix_file, resolution=int(resolution))

        my_chrom = Chromosome(name=expt_name, centromere_search=True)
        my_chrom.add_experiment(expt_name, hic_data=chr_hic_data, resolution=int(resolution))

        # Run core TADbit function to find TADs on each expt.
        # For the current dataset required 61GB of RAM
        my_chrom.find_tad(expt_name, n_cpus=15)

        exp = my_chrom.experiments[expt_name]
        exp.write_tad_borders(savedata=tad_file)

        return True


    def run(self, input_files, metadata):
        """
        The main function to the predict TAD sites for a given resolution from
        the Hi-C matrix

        Parameters
        ----------
        input_files : list
            adj_list : str
                Location of the adjacency list
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

        resolutions = [1000000]
        if 'resolutions' in metadata:
            resolutions = metadata['resolutions']

        normalized = False
        if 'normalized' in metadata:
            normalized = metadata['normalized']

        # input and output share most metadata
        output_metadata = {}

        root_name = adj_list.split("/")

        matrix_files = []
        tad_files = {}

        for resolution in resolutions:
            hic_data = load_hic_data_from_reads(adj_list, resolution=int(resolution))
            tad_files[resolution] = {}

            for chrom in hic_data.chromosomes.keys():
                save_matrix_file = "/".join(root_name[0:-1]) + '/' + metadata['expt_name'] + '_adjlist_map_' + chrom + '-' + chrom + '_' + str(resolution) + '.tsv'
                matrix_files.append(save_matrix_file)

                save_tad_file = "/".join(root_name[0:-1]) + '/' + metadata['expt_name'] + '_tad_' + chrom + '_' + str(resolution) + '.tsv'
                tad_files[resolution][chrom] = save_tad_file

                hic_data.write_matrix(save_matrix_file, (chrom, chrom), normalized=normalized)

                expt_name = metadata['expt_name'] + '_tad_' + chrom + '_' + str(resolution)

                results = self.tb_generate_tads(expt_name, save_matrix_file, resolution, save_tad_file)

        # Step to merge all the TAD files into a single bed file
        tad_bed_file = "/".join(root_name[0:-1]) + '/tads.tsv'

        fo = open(tad_bed_file, 'w')
        for resolution in tad_files.keys():
            for chrom in tad_files[resolution].keys():
                fi_tmp = open(tad_files[resolution][chrom], 'r')
                for line in fi_tmp:
                    line = line.split("\t")
                    line[-1] = line[-1].rstrip()
                    fo.write(str(chrom) + "\t" + line[1] + "\t" + line[2] + "\tTADs_" + str(resolution) + "\t" + line[3] + "\t.\n")
                fi_tmp.close()
        fo.close()

        return ([tad_bed_file], output_metadata)

# ------------------------------------------------------------------------------
