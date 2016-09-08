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


try :
    import pysam
except ImportError :
    print "[Error] Cannot import \"pysam\" package. Have you installed it?"
    exit(-1)


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
    
    
    def getGenomeFile(self, data_dir, species, assembly):
        """
        Function for downloading and extracting the CDNA files from the ensembl FTP
        """
        
        file_name = data_dir + '/' + species + '_' + assembly + '/' + species + '.' + assembly + '.cdna.all.fa.gz'
        
        if os.path.isfile(file_name) == False:
            cdna_file = urllib2.urlopen(
            'ftp://ftp.ensembl.org/pub/current_fasta/' + species + '/cdna/' + species + '.' + assembly + '.cdna.all.fa.gz')
            
            CHUNK = 16 * 1024
                    
            with open(file_name, 'wb') as fp:
                while True:
                    chunk = cdna_file.read(CHUNK)
                    if not chunk: break
                    fp.write(chunk)
            
            self.bwa_index_genome(file_name)
    
    
    def bwa_index_genome(self, genome_file):
        """
        Create an index of the 
        """
        command_line = 'bwa index ' + genome_file
        
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()
        
        
    def bwa_align_reads(self, genome_file, data_dir, run_id):
        """
        Map the reads to the genome
        """
        
        reads_file = data_dir + '/' + run_id + '.fq'
        intermediate_file = data_dir + '/' + run_id + '.sai'
        intermediate_sam_file = data_dir + '/' + run_id + '.sam'
        output_bam_file = data_dir + '/' + run_id + '.bam'
        
        command_lines = [
            'bwa aln -q 5 -f ' + intermediate_file + ' ' + genome_file + ' ' + reads_file,
            'bwa samse -r "read group information" -f ' + intermediate_sam_file  + ' ' + genome_file + ' ' + intermediate_file + ' ' + reads_file,
            'samtools view -b -o ' + output_bam_file + ' ' + intermediate_sam_file
        ]
        
        for command_line in command_lines:
            args = shlex.split(command_line)
            p = subprocess.Popen(args)
            p.wait()
     
    
    def merge_bam(self, data_dir, final_id, run_ids=[]):
        """
        Merge together all the bams in a directory and sort to create the final
        bam ready to be filtered
        
        If run_ids is blank then the function looks for all bam files in the
        data_dir
        """
        out_bam_file = data_dir + '/' + final_id + '.bam'
        
        if len(run_ids) == 0:
            bam_files = [f for f in listdir(data_dir) if f.endswith(("sai"))]
        else:
            bam_files = [f + ".bam" for f in run_ids]
        
        for bam in bam_files:
            bam_sort_files.append(bam)
            bam_merge_files.append(["-o", bam + ".sorted.bam", "-T", bam + ".bam_sort", bam])
        
        map(pysam.sort, bam_sort_files)
        
        pysam.merge(out_bam_file, bam_merge_files)
    
        pysam.sort("-o", out_bam_file + '.sorted.bam', "-T", out_bam_file + ".bam_sort", out_bam_file)
    
        pysam.index(out_bam_file)
    
    
    def biobambam_filter_alignments(self, data_dir, run_id):
        """
        Sorts and filters the bam file.
        
        It is important that all duplicate alignments have been removed. This
        can be run as an intermediate step, but should always be run as the 
        """
        command_line = 'bamsormadup ' + data_dir + '/' + run_id + '.bam'
        
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()
    
    
    def macs2_peak_calling(self, data_dir, bam_file)
        
    
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
    
