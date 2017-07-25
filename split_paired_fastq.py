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

import argparse

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
fqr.openFastQ(file1, file2)
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

fqr.closeFastQ()
fqr.closeOutputFiles()

print count_r1, r1_id[0], count_r2, r2_id[0], count_r3
