#!/usr/bin/env python

import argparse, os.path, sys

from pytadbit import Chromosome

from fastq2adjacency import fastq2adjacency

# Set up the command line parameters
parser = argparse.ArgumentParser(description="Load adjacency list into HDF5 file")
parser.add_argument("--genome", help="Genome name") #             default="GCA_000001405.22")
parser.add_argument("--dataset", help="Name of the dataset") #    default="GSE63525")
parser.add_argument("--sra_id", help="SRA ID") #                  default="SRR1658572")
parser.add_argument("--library", help="Library") #                default="HIC003")
parser.add_argument("--enzyme_name", help="Enzyme name (HboI)") # default="HboI")
parser.add_argument("--resolution", help="Resolution") #          default=1000000)
parser.add_argument("--tmp_dir", help="Temporary data dir")
parser.add_argument("--data_dir", help="Data directory; location of SRA FASTQ files")

# Get the matching parameters from the command line
args = parser.parse_args()

if len(sys.argv) < 9:
    parser.print_help()
    sys.exit(2)

genome      = args.genome
dataset     = args.dataset
sra_id      = args.sra_id
library     = args.library
enzyme_name = args.enzyme_name
resolution  = args.resolution
tmp_dir     = args.tmp_dir
data_dir    = args.data_dir

expt_name = "_".join([sra_id, library, resolution])

f2a = fastq2adjacency()

f2a.set_params(genome, dataset, sra_id, library, enzyme_name, resolution, tmp_dir, data_dir)
f2a.load_hic_read_data()

# initiate a chromosome object that will store all Hi-C data and analysis
my_chrom = Chromosome(name=expt_name, centromere_search=True)

# load Hi-C data
my_chrom.add_experiment(expt_name, hic_data=f2a.hic_data, resolution=resolution)

# Filter and normalise the matrix
# This should have already been done as part of the generate_adjacency_matrix.py
# script
#
# my_chrom.experiments['First Hi-C experiment'].filter_columns()
# my_chrom.experiments['First Hi-C experiment'].normalize_hic(iterations=0, max_dev=0.1)

# Run core TADbit function to find TADs on each expt.
# For the current dataset required 61GB of RAM
my_chrom.find_tad(expt_name, n_cpus=8)


exp = my_chrom.experiments[expt_name]
tad_file = data_root + '/' + dataset + '/' + library + '/' + expt_name + '_tads_' + str(resolution) + '.tsv'
tad_out = open(tad_file)
tad_out.write(exp.write_tad_borders())
tad_out.close()


