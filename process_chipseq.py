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

import shlex, subprocess, os.path

from pycompss.api.task import task
from pycompss.api.parameter import *

class process_chipseq:
    """
    Functions for downloading and processing Chip-seq FastQ files. Files are
    downloaded from the European Nucleotide Archive (ENA), then aligned,
    filtered and analysed for peak calling
    """
    
    def __init__ (self):
        """
        Initialise the module
        """
        self.ready = ""
    
if __name__ == "__main__":
    import sys
    import os
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="RNA-seq peak calling")
    parser.add_argument("--species", help="Species (Homo_sapiens)")
    parser.add_argument("--assembly", help="Assembly (GRCh38)")
    parser.add_argument("--project_id", help="Project ID of the dataset ()")
    parser.add_argument("--run_id", help="Experiment run ID of the dataset ()")
    parser.add_argument("--data_dir", help="Data directory; location to download ERR FASTQ files and save results")

    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    project    = args.project_id
    ena_err_id = args.run_id
    species    = args.species
    assembly   = args.assembly
    data_dir   = args.data_dir
    
    prs = process_rnaseq()
    
    try:
        os.makedirs(data_dir)
        
        data_dir += "/"
        os.makedirs(data_dir + project)
        os.makedirs(data_dir + species + "_" + assembly)
    except:
        pass
    
    # Optain the FastQ files
    in_files = pcs.getFastqFiles(ena_err_id, data_dir)
    in_file1 = in_files[0]
    in_file2 = in_files[1]
    
    # Get the assembly
    
    
    # Run BWA
    
    
    # Run BioBamBam to sort and highlight duplicates
    
    
    # Peak Calling with MACs
    
