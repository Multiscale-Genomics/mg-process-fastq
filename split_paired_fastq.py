#!/usr/bin/env python

import argparse, os.path, sys

from fastqreader import fastqreader
    
        
# Set up the command line parameters
parser = argparse.ArgumentParser(description="Load adjacency list into HDF5 file")
parser.add_argument("--input_1", help="File 1")
parser.add_argument("--input_2", help="File 2")
parser.add_argument("--output_tag", help="Inserted before the file descriptor and after the file name: e.g. 'matching' would convert file_id-1.fastq to file_id-1.matching.fastq")

args = parser.parse_args()
file1 = args.input_1
file2 = args.input_2
tag = args.output_tag

fqr = fastqreader()
fqr.openPairedFastQ(file1, file2)
fqr.createOutputFiles(tag)

r1 = fqr.next(1)
r2 = fqr.next(2)

count_r1 = 0
count_r2 = 0
count_r3 = 0

while fqr.eof(1) == False and fqr.eof(2) == False:
    r1_id = r1["id"].split(" ")
    r2_id = r2["id"].split(" ")
    
    if r1_id[0] == r2_id[0]:
        fqr.writeOutput(r1, 1)
        fqr.writeOutput(r2, 2)
        
        r1 = fqr.next(1)
        r2 = fqr.next(2)
        
        count_r1 += 1
        count_r2 += 1
        count_r3 += 1
    elif r1_id[0] < r2_id[0]:
        r1 = fqr.next(1)
        count_r1 += 1
    else:
        r2 = fqr.next(2)
        count_r2 += 1
    
    if count_r3 % 1000000 == 0:
        print count_r1, r1_id[0], count_r2, r2_id[0], count_r3
        fqr.incrementOutputFiles()

fqr.closePairedFastQ()
fqr.closeOutputFiles()

print count_r1, r1_id[0], count_r2, r2_id[0], count_r3
