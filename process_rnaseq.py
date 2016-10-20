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

import argparse, urllib2, gzip, shutil, shlex, subprocess, os.path

from common import common

from pycompss.api.task import task
from pycompss.api.parameter import *

class process_rnaseq:
    """
    Functions for downloading and processing RNA-seq FastQ files. Files are
    downloaded from the European Nucleotide Archive (ENA), then they are mapped
    to quantify the amount of cDNA
    """
    
    def __init__ (self):
        """
        Initialise the module
        """
        self.ready = ""
    
    
    def run_kallisto_indexing(self, data_dir, species, assembly):
        """
        Runs the Kallisto index program to generate a list of the indexes for each of the cDNAs
        """
        
        cdna_file = data_dir + species + '_' + assembly + '/' + species + '.' + assembly + '.cdna.all.fa.gz'
        cdna_idx_file = cdna_file + '.idx'
        
        command_line = 'kallisto index -i ' + cdna_idx_file + ' ' + cdna_file
        
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()
        
    
    
    def run_kallisto_quant(self, data_dir, species, assembly, project, fastq = [], single = False):
        """
        Kallisto function to map the paired end FastQ files to the cDNAs and generate the matching quatification files.
        """
        
        cdna_file = data_dir + species + '_' + assembly + '/' + species + '.' + assembly + '.cdna.all.fa.gz'
        cdna_idx_file = cdna_file + '.idx'
        
        command_line = 'kallisto quant -i ' + cdna_idx_file + ' -o ' + data_dir + project
        if single == True:
            fq_stats = self.seq_read_stats(fastq[0])
            command_line += ' --single -l ' + str(fq_stats['mean']) + ' -s ' + str(fq_stats['std']) + ' ' + fastq[0]
        else:
            command_line += ' ' + fastq[0] + ' ' + fastq[1]
        
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()
    
    
    def seq_read_stats(self, file_in):
        """
        Calculate the mean and standard deviation of the reads in a fastq file
        """
        
        from numpy import std
        from numpy import mean

        total_len = []
        with open(file_in, 'r') as f:
            forthlines = itertools.islice(f, 1, None, 4)
            for line in forthlines:
                line = line.rstrip()
                total_len.append(len(line))

        s = std(total_len)
        m = mean(total_len)

        return {'mean' : m, 'std' : s}
        

if __name__ == "__main__":
    import sys
    import os
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="Parse RNA-seq for expression analysis")
    parser.add_argument("--species", help="Species (Homo_sapiens)")
    parser.add_argument("--assembly", help="Assembly (GRCh38)")
    parser.add_argument("--project_id", help="Project ID of the dataset (PRJEB2445)")
    parser.add_argument("--run_id", help="Experiment run ID of the dataset (ERR030872)")
    parser.add_argument("--data_dir", help="Data directory; location to download ERR FASTQ files and save results")
    parser.add_argument("--local", help="Directory and data files already available", default=0)

    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    project    = args.project_id
    ena_err_id = args.run_id
    species    = args.species
    assembly   = args.assembly
    data_dir   = args.data_dir
    local      = args.local
    
    prs = process_rnaseq()
    cf = common()
    
    try:
        os.makedirs(data_dir)
    except:
        pass
    
    if data_dir[-1] != "/":
        data_dir += "/"
        #data_dir gzf.replace('.fastq.gz', '.fastq')
    
    try:
        os.makedirs(data_dir + project)
    except:
        pass
    
    try:
        os.makedirs(data_dir + species + "_" + assembly)
    except:
        pass
    
    
    
    # Optain the paired FastQ files
    if (local == 0):
        in_files = cf.getFastqFiles(ena_err_id, data_dir)
    else:
        in_files = [f for f in os.listdir(data_dir + project) if re.match(run_id, f)]
    
    # Get the cDNA
    cf.getcDNAFiles(data_dir, species, assembly)
    
    # Index the cDNA
    prs.run_kallisto_indexing(data_dir, species, assembly)
    
    # Quantify RNA-seq
    if len(in_files) == 1:
        prs.run_kallisto_quant(data_dir, species, assembly, project, in_files, True)
    else:
        prs.run_kallisto_quant(data_dir, species, assembly, project, in_files)
    
    
