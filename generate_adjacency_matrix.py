#!/usr/bin/env python

"""
Copyright 2016 EMBL-European Bioinformatics Institute

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
parser.add_argument("--data_dir", help="Data directory; location to download SRA FASTQ files and save results")

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


f2a = fastq2adjacency()

windows1 = ((1,25), (1,50), (1,75),(1,100))
windows2 = ((1,25), (1,50), (1,75),(1,100))
#windows2 = ((101,125), (101,150), (101,175),(101,200))
f2a.set_params(genome, dataset, sra_id, library, enzyme_name, resolution, tmp_dir, data_dir, same_fastq=False, windows1=windows1, windows2=windows2)

f2a.getFastqData()

map(f2a.mapWindows, [1, 2])

f2a.parseGenomeSeq()

f2a.parseMaps()

f2a.mergeMaps()

f2a.filterReads(conservative=True)

# It is at this point that the resolution is used.
f2a.load_hic_read_data()

f2a.normalise_hic_data()

f2a.save_hic_data()

