#!/usr/bin/python
# -*- coding: utf-8 -*-
'''process Hi-C paired end FastQ files'''
import argparse, time
from pycompss.api.task import task
from pycompss.api.parameter import *

from fastq2adjacency import fastq2adjacency

@task()
def main(genome, dataset, sra_id, library, enzyme_name, resolution, tmp_dir, data_dir, expt):
    f2a = fastq2adjacency()
    f2a.set_params(genome, dataset, sra_id, library, enzyme_name, resolution, tmp_dir, data_dir, expt, same_fastq=False, windows1=windows1, windows2=windows2)
    
    f2a.getFastqData()
    
    map(MapReads2Genome, [1, 2])

    f2a.parseGenomeSeq()

    f2a.parseMaps()

    f2a.mergeMaps()

    f2a.filterReads(conservative=True)

    # It is at this point that the resolution is used.
    f2a.load_hic_read_data()

    f2a.normalise_hic_data()

    f2a.save_hic_data()

@task()
def MapReads2Genome(f2a, window):
    f2a.mapWindows(window)


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

    if len(sys.argv) < 9:
        parser.print_help()
        sys.exit(2)

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
        more_params = [[genome, dataset, line[0], line[1], line[2], resolution, tmp_dir, data_dir, expt, False, windows1, windows2] for expt in resolutions]
        loading_list += more_params


    map(main, loading_list)
