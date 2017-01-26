#!/usr/bin/env python

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

print "Creating intermediate save files ..."
f2a.save_hic_split_data()

print "Generating TADS:"
chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']
#for chrom in range(len(chroms)):
for chrom in chroms:
    print chrom
    #f2a.generate_tads(chroms[chrom])
    f2a.generate_tads(chrom)



