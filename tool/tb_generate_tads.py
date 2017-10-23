"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.

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
    from pycompss.api.parameter import FILE_IN, FILE_OUT, FILE_INOUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
    # from pycompss.api.constraint import constraint
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT, FILE_INOUT, IN
    from dummy_pycompss import task
    from dummy_pycompss import compss_wait_on
    #from dummy_pycompss import constraint

from basic_modules.tool import Tool

from pytadbit import Chromosome
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

    @task(adj_list=FILE_IN, resolution=IN, normalized=IN, returns=list)
    def tb_hic_chr(self, adj_list, resolution):
        """
        """
        print("TB LOADED HIC MATRIX")
        hic_data = load_hic_data_from_reads(adj_list, resolution=int(resolution))

        print("TB LOADED HIC MATRIX")

        return hic_data.chromosomes.keys()


    @task(expt_name=IN, adj_list=FILE_IN, chrom=IN, resolution=IN, normalized=IN, tad_file=FILE_OUT)
    # @constraint(ProcessorCoreCount=16)
    def tb_generate_tads(self, expt_name, adj_list, chrom, resolution, normalized, tad_file):
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
        #chr_hic_data = read_matrix(matrix_file, resolution=int(resolution))

        print("TB TAD GENERATOR:", expt_name, adj_list, chrom, resolution, normalized, tad_file)

        hic_data = load_hic_data_from_reads(adj_list, resolution=int(resolution))

        if normalized is False:
            hic_data.normalize_hic(iterations=9, max_dev=0.1)

        save_matrix_file = adj_list + "_" + str(chrom) + "_tmp.txt"
        hic_data.write_matrix(save_matrix_file, (chrom, chrom), normalized=True)

        chr_hic_data = hic_data.get_matrix((chrom, chrom))
        print("TB - chr_hic_data:", chr_hic_data)

        my_chrom = Chromosome(name=chrom, centromere_search=True)
        my_chrom.add_experiment(expt_name, hic_data=save_matrix_file, resolution=int(resolution))

        # Run core TADbit function to find TADs on each expt.
        my_chrom.find_tad(expt_name, n_cpus=15)

        exp = my_chrom.experiments[expt_name]
        exp.write_tad_borders(savedata=tad_file + ".tmp")

        with open(tad_file, "wb") as f_out:
            with open(tad_file + ".tmp", "rb") as f_in:
                f_out.write(f_in.read())

        return True

    @task(input_file=FILE_IN, chrom=IN, resolution=IN, output_file=FILE_INOUT)
    def tb_merge_tad_files(self, input_file, chrom, resolution, output_file):
        """
        """
        with open(output_file, 'a') as f_out:
            with open(input_file, 'r') as f_in:
                for line in f_in:
                    line = line.split("\t")
                    line[-1] = line[-1].rstrip()
                    f_out.write(str(chrom) + "\t" + line[1] + "\t" + line[2] + "\tTADs_" + str(resolution) + "\t" + line[3] + "\t.\n")

        return True

    def run(self, input_files, output_files, metadata=None):
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

        tad_files = {}

        results =[]
        for resolution in resolutions:
            print("TB LOADING Hi-C:", adj_list, resolution, normalized)
            hic_data_chr = self.tb_hic_chr(adj_list, resolution)
            hic_data_chr = compss_wait_on(hic_data_chr)
            print("TB LOADED Hi-C!")

            print("TB CHROMOSOMES", hic_data_chr)

            tad_files[resolution] = {}

            for chrom in hic_data_chr:
                save_tad_file = "/".join(root_name[0:-1]) + '/' + metadata['expt_name'] + '_tad_' + chrom + '_' + str(resolution) + '.tsv'
                tad_files[resolution][chrom] = save_tad_file

                expt_name = metadata['expt_name'] + '_tad_' + chrom + '_' + str(resolution)

                print("TB Generate TADS:", resolution, normalized)
                results.append(
                    self.tb_generate_tads(
                        expt_name, adj_list, chrom, resolution, normalized,
                        save_tad_file
                    )
                )

        results = compss_wait_on(results)

        # Step to merge all the TAD files into a single bed file
        tad_bed_file = "/".join(root_name[0:-1]) + '/' + metadata['expt_name'] + '_tads.tsv'

        print("TB tad_files:", tad_files)

        # Step to merge all the TAD files into a single bed file
        tad_bed_file = "/".join(root_name[0:-1]) + '/' + metadata['expt_name'] + '_tads.tsv'
        f_prepare = open(tad_bed_file, 'w')
        f_prepare.close()

        for resolution in tad_files:
            for chrom in tad_files[resolution]:
                results = self.tb_merge_tad_files(
                    tad_files[resolution][chrom],
                    chrom,
                    resolution,
                    tad_bed_file
                )
                results = compss_wait_on(results)

        return ([tad_bed_file], output_metadata)

# ------------------------------------------------------------------------------
