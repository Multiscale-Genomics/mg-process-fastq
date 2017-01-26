#!/usr/bin/python

"""
Copyright 2017 EMBL-European Bioinformatics Institute

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

# -*- coding: utf-8 -*-
'''process Hi-C paired end FastQ files'''
import argparse, time

try :
    from pycompss.api.parameter import *
    from pycompss.api.task import task
except ImportError :
    print "[Warning] Cannot import \"pycompss\" API packages."
    print "          Using mock decorators."
    
    from dummy_pycompss import *

from common import common

class process_hic:
    #@task(params = IN)
    def main(self, params):
        """
        Initial grouping to download, parse and filter the individual
        experiments.
        
        Returns: None
        
        Output: Raw counts for the experiment in a HiC adjacency matrix saved to
                the tmp_dir
        """
        from common import common
        from fastq2adjacency import fastq2adjacency
        
        genome      = params[0]
        dataset     = params[1]
        sra_id      = params[2]
        library     = params[3]
        enzyme_name = params[4]
        resolution  = params[5]
        tmp_dir     = params[6]
        data_dir    = params[7]
        expt        = params[8]
        same_fastq  = params[9]
        windows1    = params[10]
        windows2    = params[11]
        
        print "Got Params"
        
        print sra_id, library, resolution, time.time()
        
        f2a = fastq2adjacency()
        f2a.set_params(genome, dataset, sra_id, library, enzyme_name, resolution, tmp_dir, data_dir, expt, same_fastq, windows1, windows2)
        
        print "Set Params"
        cf = common()
        in_files = cf.getFastqFiles(sra_id, data_dir)
        
        map(f2a.mapWindows, [1, 2])

        f2a.parseGenomeSeq()

        f2a.parseMaps()

        f2a.mergeMaps()

        f2a.filterReads(conservative=True)

        # It is at this point that the resolution is used.
        f2a.load_hic_read_data()
        f2a.save_hic_split_data()
        chroms = f2a.get_chromosomes()
        
        for chrom in chroms:
            f2a.generate_tads(chrom)
        
        f2a.normalise_hic_data()
        f2a.save_hic_data()

    def merge_adjacency_data(self, adj_list):
        """
        Merged the HiC filtered data into a single dataset.
        
        Input:   list of all the params in a list of lists.
        
        Returns: None
        
        Output:  The merged adjacency file, split files for each chrA vs chrB 
                 combination and a file for each chromosome with the positions
                 of the predicted TADs
        """
        from pycompss.api.api import compss_wait_on
        from fastq2adjacency import fastq2adjacency
        
        f2a = fastq2adjacency()
        genome      = params[0][0]
        dataset     = params[0][1]
        sra_id      = params[0][2]
        library     = params[0][3]
        enzyme_name = params[0][4]
        resolution  = params[0][5]
        tmp_dir     = params[0][6]
        data_dir    = params[0][7]
        expt        = params[0][8]
        same_fastq  = params[0][9]
        windows1    = params[0][10]
        windows2    = params[0][11]
        f2a.set_params(genome, dataset, sra_id, library, enzyme_name, resolution, tmp_dir, data_dir, expt, same_fastq, windows1, windows2)
        
        f2a.load_hic_read_data()
        
        new_hic_data = f2a.hic_data
        
        for i in range(1,len(adj_list)):
            f2a = fastq2adjacency()
            genome      = params[i][0]
            dataset     = params[i][1]
            sra_id      = params[i][2]
            library     = params[i][3]
            enzyme_name = params[i][4]
            resolution  = params[i][5]
            tmp_dir     = params[i][6]
            data_dir    = params[i][7]
            expt        = params[i][8]
            same_fastq  = params[i][9]
            windows1    = params[i][10]
            windows2    = params[i][11]
            f2a.set_params(genome, dataset, sra_id, library, enzyme_name, resolution, tmp_dir, data_dir, expt, same_fastq, windows1, windows2)
            
            f2a.load_hic_read_data()
        
            new_hic_data += f2a.hic_data
        
        f2a = fastq2adjacency()
        genome      = params[0][0]
        dataset     = params[0][1]
        sra_id      = params[0][2] + "_all"
        library     = "MERGED"
        enzyme_name = "RANGE"
        resolution  = params[0][5]
        tmp_dir     = params[0][6]
        data_dir    = params[0][7]
        expt        = params[0][8]
        same_fastq  = params[0][9]
        windows1    = params[0][10]
        windows2    = params[0][11]
        f2a.set_params(genome, dataset, sra_id, library, enzyme_name, resolution, tmp_dir, data_dir, expt, same_fastq, windows1, windows2)
        
        f2a.hic_data = new_hic_data
        f2a.save_hic_data()
        f2a.save_hic_hdf5()
        f2a.save_hic_split_data()
        
        tad_done = []
        for chrom in f2a.get_chromosomes():
            tad_done.append(call_tads(genome, dataset, sra_id, library, enzyme_name, resolution, tmp_dir, data_dir, expt, same_fastq, windows1, windows2, chrom))
        tad_done.compss_wait_on(tad_done)
    
    #@task(genome = IN, dataset = IN, sra_id = IN, library = IN, enzyme_name = IN, resolution = IN, tmp_dir = IN, data_dir = IN, expt = IN, same_fastq = IN, windows1 = IN, windows2 = IN, chrom = IN, returns = int)
    def call_tads(self, genome, dataset, sra_id, library, enzyme_name, resolution, tmp_dir, data_dir, expt, same_fastq, windows1, windows2, chrom):
        """
        TAD calling for a given dataset and chromosome. This should be run from
        the merge step, but can be run individually. Relies on the split
        adjacency files having already been created.
        
        Input:   params for the f2a set up and the chromosome for analysis
        
        Returns: None
        
        Output:  File containing the predicted TADs.
        """
        from fastq2adjacency import fastq2adjacency
        
        f2a = fastq2adjacency()
        f2a.set_params(genome, dataset, sra_id, library, enzyme_name, resolution, tmp_dir, data_dir, expt, same_fastq, windows1, windows2)
        f2a.generate_tads(chrom)
        
        return 1

    def merge_hdf5_files(self, genome, dataset, resolutions, data_dir):
        """
        Merges the separate HDF5 files with each of the separate resolutions
        into a single 
        """
        
        #f = h5py.File(filename, "a")
        #dset = f.create_dataset(str(self.resolution), (dSize, dSize), dtype='int32', chunks=True, compression="gzip")
        #dset[0:dSize,0:dSize] += d
        #f.close()
        
        f_out = data_dir + genome + "_" + dataset + ".hdf5"
        final_h5 = h5py.File(f_out, "a")
        for resolution in resolutions:
            f = data_dir + genome + "_" + dataset + "_" + str(resolution) + ".hdf5"
            fin = h5py.File(f, "r")
            hdf5.h5o.copy(fin, str(resolution), final_h5, str(resolution))
            fin.close()
        final_h5.close()

if __name__ == "__main__":
    import sys
    import os
    
    start = time.time()
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="Load adjacency list into HDF5 file")
    #parser.add_argument("--genome", help="Genome name") #             default="GCA_000001405.22")
    parser.add_argument("--species", help="Species (homo_sapiens)")
    parser.add_argument("--assembly", help="Assembly (GRCh38)")
    parser.add_argument("--dataset", help="Name of the dataset") #    default="GSE63525")
    parser.add_argument("--expt_name")
    parser.add_argument("--expt_list", help="TSV detailing the SRA ID, library and restriction enzymeused that are to be treated as a single set")
    parser.add_argument("--tmp_dir", help="Temporary data dir")
    parser.add_argument("--data_dir", help="Data directory; location to download SRA FASTQ files and save results")

    # Get the matching parameters from the command line
    args = parser.parse_args()

    genome      = args.genome
    dataset     = args.dataset
    expt_name   = args.expt_name
    expt_list   = args.expt_list
    tmp_dir     = args.tmp_dir
    data_dir    = args.data_dir
    
    # A default value is only required for the first few steps to generate the
    # intial alignments and prepare the HiC data for loading. The resolutions
    # get passed later on when the pipeline splits and handles each on
    # individually via the process_block_size() function.
    #resolutions = [1000, 2500, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000, 10000000]
    resolutions = [1000000, 10000000]
    
    windows1 = ((1,25), (1,50), (1,75),(1,100))
    windows2 = ((1,25), (1,50), (1,75),(1,100))
    
    f = open(expt_list, "r")

    loading_list = []
    for line in f:
        line = line.rstrip()
        line = line.split("\t")
        
        #                                sra_id,  library, enzyme_name
        more_params = [[genome, dataset, line[0], line[1], line[2], resolution, tmp_dir, data_dir, expt_name, False, windows1, windows2] for resolution in resolutions]
        less_params = [genome, dataset, line[0], line[1], line[2], 1000, tmp_dir, data_dir, expt_name, False, windows1, windows2]
        more_loading_list += more_params
        less_loading_list += less_params

    print more_loading_list
    
    cf = common()
    
    # Get the assembly
    genome_fa = cf.getGenomeFile(data_dir, species, assembly)
    
    hic = process_hic()
    
    # Downloads the FastQ files and then maps then to the genome.
    map(hic.main, less_loading_list)
    
    # Generates the final adjacency matrix for a given resolutions
    hic.merge_adjacency_data(more_loading_list)
    
    # This merges the final set of HDF5 files into a single file ready for the
    # REST API
    hic.merge_hdf5_files(genome, dataset, resolutions, data_dir)
    
