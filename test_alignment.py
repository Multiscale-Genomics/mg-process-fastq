#!/usr/bin/python

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

import argparse, urllib2, gzip, shutil, shlex, subprocess, os.path, json

from pycompss.api.task import task
from pycompss.api.parameter import *

from common import common

if __name__ == "__main__":
    import sys
    import os
    
    from pycompss.api.api import compss_wait_on
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="Genome assembly and indexing")
    parser.add_argument("--species", help="Species (homo_sapiens)")
    parser.add_argument("--assembly", help="Assembly (GRCh38)")
    parser.add_argument("--project_id", help="Project ID of the dataset")
    parser.add_argument("--run_ids", help="JSON file with list of the experiment run IDs and background data (if available) of the dataset")
    parser.add_argument("--data_dir", help="Data directory; location to download ERR FASTQ files and save results")

    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    project     = args.project_id
    run_id_file = args.run_ids
    species     = args.species
    assembly    = args.assembly
    data_dir    = args.data_dir
    
    pcs = process_chipseq()
    cf = common()
    
    if data_dir[-1] != "/":
        data_dir += "/"
    
    try:
        os.makedirs(data_dir)
    except:
        pass

    try:
        os.makedirs(data_dir + project)
    except:
        pass
    
    try:
        os.makedirs(data_dir + species + "_" + assembly)
    except:
        pass
    
    # Get the assembly
    genome_fa = cf.getGenomeFile(data_dir, species, assembly)
