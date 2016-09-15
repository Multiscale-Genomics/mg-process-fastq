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
    
    
    def getcDNAFiles(self, data_dir, species, assembly):
        """
        Function for downloading and extracting the CDNA files from the ensembl FTP
        """
        
        file_name = data_dir + species + '_' + assembly + '/' + species + '.' + assembly + '.cdna.all.fa.gz'
        
        if os.path.isfile(file_name) == False:
            cdna_file = urllib2.urlopen(
            'ftp://ftp.ensembl.org/pub/current_fasta/' + species.lower() + '/cdna/' + species + '.' + assembly + '.cdna.all.fa.gz')
            
            CHUNK = 16 * 1024
                    
            with open(file_name, 'wb') as fp:
                while True:
                    chunk = cdna_file.read(CHUNK)
                    if not chunk: break
                    fp.write(chunk)
    
    
    def getFastqFiles(self, ena_err_id, data_dir):
        """
        Function for downloading and extracting the FastQ files from the ENA
        """
        
        f_index = urllib2.urlopen(
        'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=' + str(ena_err_id) + '&result=read_run&fields=study_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp&download=txt')
        data = f_index.read()
        rows = data.split("\n")
        row_count = 0
        files = []
        gzfiles  = []
        for row in rows:
            if row_count == 0:
                row_count += 1
                continue
            
            row = row.rstrip()
            row = row.split("\t")
            
            if len(row) < 6:
                continue
            
            print row
            project = row[0]
            srr_id = row[1]
            fastq_files = row[6].split(';')
            row_count += 1
            
            for fastq_file in fastq_files:
                file_name = fastq_file.split("/")
                print fastq_file
                print data_dir + project + "/" + file_name[-1]
                print file_name[-1]
                
                req = urllib2.urlopen("ftp://" + fastq_file)
                CHUNK = 16 * 1024
                
                files.append(data_dir + project + "/" + file_name[-1].replace('.fastq.gz', '.fastq'))
                gzfiles.append(data_dir + project + "/" + file_name[-1])
                
                with open(data_dir + project + "/" + file_name[-1], 'wb') as fp:
                    while True:
                        chunk = req.read(CHUNK)
                        if not chunk: break
                        fp.write(chunk)
        
        for gzf in gzfiles:
            with gzip.open(gzf, 'rb') as f_in, open(gzf.replace('.fastq.gz', '.fastq'), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(gzf)
        
        return files
    
    
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
        
    
    
    def run_kallisto_quant(self, data_dir, species, assembly, project, fastq = []):
        """
        Kallisto function to map the paired end FastQ files to the cDNAs and generate the matching quatification files.
        """
        
        cdna_file = data_dir + species + '_' + assembly + '/' + species + '.' + assembly + '.cdna.all.fa.gz'
        cdna_idx_file = cdna_file + '.idx'
        
        command_line = 'kallisto quant -i ' + cdna_idx_file + ' -o ' + data_dir + project + ' ' + fastq[0] + ' ' + fastq[1]
        
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()
        

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
    except:
        pass
    
    data_dir += "/"
    
    try:
        os.makedirs(data_dir + project)
    except:
        pass
    
    try:
        os.makedirs(data_dir + species + "_" + assembly)
    except:
        pass
    
    # Optain the paired FastQ files
    in_files = prs.getFastqFiles(ena_err_id, data_dir)
    
    # Get the cDNA
    prs.getcDNAFiles(data_dir, species, assembly)
    
    # Index the cDNA
    prs.run_kallisto_indexing(data_dir, species, assembly)
    
    # Quantify RNA-seq
    prs.run_kallisto_quant(data_dir, species, assembly, project, in_files)
    
    
