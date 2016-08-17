#!/usr/bin/python
# -*- coding: utf-8 -*-
'''process Hi-C paired end FastQ files'''
import argparse, time
from pycompss.api.task import task
from pycompss.api.parameter import *

class process_hic:
    #@task(params = IN)
    def main(self, params):
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
        
        f2a.getFastqData()
        
        map(f2a.mapWindows, [1, 2])

        f2a.parseGenomeSeq()

        f2a.parseMaps()

        f2a.mergeMaps()

        f2a.filterReads(conservative=True)

        # It is at this point that the resolution is used.
        #f2a.load_hic_read_data()
        
        #f2a.save_hic_split_data()
        
        #chrom = f2a.get_chromosomes()

        #f2a.normalise_hic_data()

        #f2a.save_hic_data()

    #@task()
    #def MapReads2Genome(f2a, window):
    #    f2a.mapWindows(window)
    
    def merge_adjacency_data(self, adj_list):
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
        
        f2a.save_hic_split_data()
        
        map(f2a.generate_tads, f2a.get_chromosomes())


if __name__ == "__main__":
    import sys
    import os
    from pycompss.api.api import compss_wait_on
    
    start = time.time()
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="Load adjacency list into HDF5 file")
    parser.add_argument("--genome", help="Genome name") #             default="GCA_000001405.22")
    parser.add_argument("--dataset", help="Name of the dataset") #    default="GSE63525")
    parser.add_argument("--expt_name")
    parser.add_argument("--expt_list", help="TSV detailing the SRA ID, library and restriction enzymeused that are to be treated as a single set")
    parser.add_argument("--tmp_dir", help="Temporary data dir")
    parser.add_argument("--data_dir", help="Data directory; location to download SRA FASTQ files and save results")

    # Get the matching parameters from the command line
    args = parser.parse_args()

    #if len(sys.argv) < 9:
    #    parser.print_help()
    #    sys.exit(2)

    genome      = args.genome
    dataset     = args.dataset
    expt_name = args.expt_name
    expt_list = args.expt_list
    tmp_dir     = args.tmp_dir
    data_dir    = args.data_dir
    
    # A defauilt value is only required for the first few steps to generate the
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
        loading_list += more_params

    print loading_list
    
    hic = process_hic()
    map(hic.main, loading_list)
    
    
    
