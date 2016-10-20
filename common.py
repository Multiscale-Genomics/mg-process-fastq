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

from socket import error as SocketError
import errno

try :
    import pysam
except ImportError :
    print "[Error] Cannot import \"pysam\" package. Have you installed it?"
    exit(-1)

class common:
    """
    Functions for downloading and processing Mnase-seq FastQ files. Files are
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
        Function for downloading and extracting the DNA files from the ensembl FTP
        """
        
        file_name = data_dir + '/' + species + '_' + assembly + '/' + species + '.' + assembly + '.dna.toplevel.fa.gz'
        
        if os.path.isfile(file_name) == False:
            cdna_file = urllib2.urlopen(
            'ftp://ftp.ensembl.org/pub/current_fasta/' + species.lower() + '/dna/' + species[0].upper() + species[1:] + '.' + assembly + '.dna.toplevel.fa.gz')
            
            CHUNK = 16 * 1024
                    
            with open(file_name, 'wb') as fp:
                while True:
                    chunk = cdna_file.read(CHUNK)
                    if not chunk: break
                    fp.write(chunk)
            
            cdna_file.close()
            
            self.bwa_index_genome(file_name)
        
        return file_name
    
    
    def getcDNAFiles(self, data_dir, species, assembly):
        """
        Function for downloading and extracting the CDNA files from the ensembl FTP
        """
        
        file_name = data_dir + species + '_' + assembly + '/' + species + '.' + assembly + '.cdna.all.fa.gz'
        
        if os.path.isfile(file_name) == False:
            cdna_file = urllib2.urlopen(
            'ftp://ftp.ensembl.org/pub/current_fasta/' + species.lower() + '/cdna/' + species[0].upper() + species[1:] + '.' + assembly + '.cdna.all.fa.gz')
            
            CHUNK = 16 * 1024
                    
            with open(file_name, 'wb') as fp:
                while True:
                    chunk = cdna_file.read(CHUNK)
                    if not chunk: break
                    fp.write(chunk)
            
            cdna_file.close()
        
        return file_name
    
    
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
            
            project = row[0]
            srr_id = row[1]
            fastq_files = row[6].split(';')
            row_count += 1
            
            for fastq_file in fastq_files:
                file_name = fastq_file.split("/")
                print data_dir + '/' + project + "/" + file_name[-1]
                
                if os.path.isfile(data_dir + '/' + project + "/" + file_name[-1]) == False and os.path.isfile(data_dir + '/' + project + "/" + file_name[-1].replace('.fastq.gz', '.fastq')) == False:
                    req = urllib2.urlopen("ftp://" + fastq_file)
                    CHUNK = 16 * 1024
                    
                    #files.append(data_dir + '/' + project + "/" + file_name[-1].replace('.fastq.gz', '.fastq'))
                    gzfiles.append(data_dir + '/' + project + "/" + file_name[-1])
                    try:
                        with open(data_dir + '/' + project + "/" + file_name[-1], 'wb') as fp:
                            while True:
                                chunk = req.read(CHUNK)
                                if not chunk: break
                                fp.write(chunk)
                        
                        req.close()
                    except SocketError as e:
                        if e.errno != errno.ECONNRESET:
                            raise
                        print e
                        pass # Handle error here
                
                files.append(data_dir + '/' + project + "/" + file_name[-1].replace('.fastq.gz', '.fastq'))
        
        for gzf in gzfiles:
            with gzip.open(gzf, 'rb') as f_in, open(gzf.replace('.fastq.gz', '.fastq'), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(gzf)
        
        return files
    
    
    def bwa_index_genome(self, genome_file):
        """
        Create an index of the 
        """
        command_line = 'bwa index ' + genome_file
        
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()
        
        
    def bwa_align_reads(self, genome_file, data_dir, project_id, run_id):
        """
        Map the reads to the genome
        """
        
        reads_file = data_dir + '/' + project_id + '/' + run_id + '.fastq'
        intermediate_file = data_dir + '/' + project_id + '/' + run_id + '.sai'
        intermediate_sam_file = data_dir + '/' + project_id + '/' + run_id + '.sam'
        output_bam_file = data_dir + '/' + project_id + '/' + run_id + '.bam'
        
        command_lines = [
            'bwa aln -q 5 -f ' + intermediate_file + ' ' + genome_file + ' ' + reads_file,
            'bwa samse -f ' + intermediate_sam_file  + ' ' + genome_file + ' ' + intermediate_file + ' ' + reads_file,
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
        out_bam_file = data_dir + '/' + project_id + '/' + final_id + '.bam'
        
        if len(run_ids) == 0:
            bam_files = [f for f in listdir(data_dir + '/' + project_id) if f.endswith(("sai"))]
        else:
            bam_files = [f + ".bam" for f in run_ids]
        
        bam_sort_files = []
        bam_merge_files = []
        for bam in bam_files:
            bam_sort_files.append(bam)
            bam_merge_files.append(["-o", bam + ".sorted.bam", "-T", bam + ".bam_sort", bam])
        
        map(pysam.sort, bam_sort_files)
        
        pysam.merge(out_bam_file, bam_merge_files)
    
        pysam.sort("-o", out_bam_file + '.sorted.bam', "-T", out_bam_file + ".bam_sort", out_bam_file)
    
        pysam.index(out_bam_file)
