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

# -*- coding: utf-8 -*-
'''process whole genome bisulfate sequencing FastQ files'''
import argparse, time, urllib2, gzip, shutil
from pycompss.api.task import task
from pycompss.api.parameter import *
from pycompss.api.constraint import *

from common import common

from fastqreader import *
from FilterReads import *
from bs_index.wg_build import *

class process_wgbs:
    """
    Functions for downloading and processing whole genome bisulfate sequencings
    (WGBS) files. Files are downloaded from the European Nucleotide Archive
    (ENA), then filtered, aligned and analysed for points of methylation
    """

    def __init__ (self):
        """
        Initialise the module
        """
        self.ready=None


    @task(infile = FILE_IN, outfile = FILE_OUT, returns = int)
    def FilterFastQReads(self, infile, outfile):
        """
        This is optional, but removes reads that can be problematic for the
        alignment of whole genome datasets.
        
        If performing RRBS then this step can be skipped
        """
        FilterReads(infile, outfile, True)
        return 1


    @constraint(ProcessorCoreCount=8)
    @task(fasta_file = FILE_IN, aligner = IN, aligner_path = IN, ref_path = IN)
    def Builder(self, fasta_file, aligner, aligner_path, ref_path):
        """
        Function to submit the FASTA file for the reference sequence and build
        the required index file used by the aligner.
        
        This only needs to be done once, so there needs to be a check to ensure
        that if the index file have already been generated then they do no need
        to be analysed again
        """
        builder_exec = os.path.join(aligner_path,
                                    {BOWTIE   : 'bowtie-build',
                                     BOWTIE2  : 'bowtie2-build',
                                     SOAP     : '2bwt-builder',
                                     RMAP     : '' # do nothing
                                    }[aligner])

        build_command = builder_exec + { BOWTIE   : ' -f %(fname)s.fa %(fname)s',
                                         BOWTIE2  : ' -f %(fname)s.fa %(fname)s',
                                         SOAP     : ' %(fname)s.fa'
                                       }[aligner]
        
        wg_build(fasta_file, build_command, ref_path, aligner)
        

    def Splitter(self, in_file1, in_file2, tag):
        """
        Function to divide the FastQ files into separte sub files of 1000000
        sequences so that the aligner can get run in parallel.
        
        Returns: Returns a list of lists of the files that have been generated.
                 Each sub list containing the two paired end files for that
                 subset.
        """
        
        fqr = fastqreader()
        fqr.openPairedFastQ(file1, file2)
        fqr.createOutputFiles(tag)

        r1 = fqr.next(1)
        r2 = fqr.next(2)

        count_r1 = 0
        count_r2 = 0
        count_r3 = 0
        files_out = []

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
                fqr.incrementOutputFiles()
                f1 = self.fastq1.split("/")
                f1[-1] = f1[-1].replace(".fastq", "." + str(self.output_tag) + "_" + str(self.output_file_count) + ".fastq")
                f1.insert(-1, "tmp")
                
                f2 = self.fastq2.split("/")
                f2[-1] = f2[-1].replace(".fastq", "." + str(self.output_tag) + "_" + str(self.output_file_count) + ".fastq")
                f2.insert(-1, "tmp")
                
                files_out.append(["/".join(f1), "/".join(f2)])

        fqr.closePairedFastQ()
        fqr.closeOutputFiles()
        
        return files_out


    @constraint(ProcessorCoreCount=8)
    @task(input_fastq1 = FILE_IN, input_fastq2 = FILE_IN, aligner = IN, aligner_path = IN, genome_fasta = FILE_IN, returns = int)
    def Aligner(self, input_fastq1, input_fastq2, aligner, aligner_path, genome_fasta):
        """
        Alignment of the paired ends to the reference genome
        Generates bam files for the alignments
        This is performed by running the external program rather than reimplementing
        the code from the main function to make it easier when it comes to updating
        the changes in BS-Seeker2
        """
        import subprocess
        
        args = ["python", "bs_seeker2-align.py", "--input_1", input_fastq1, "--input_2", input_fastq2, "--aligner", aligner, "--path", aligner_path_path, "--genome", genome_fasta, "--bt2-p", "4"]
        p = subprocess.Popen(args)
        p.wait()
        
        return 1


    @constraint(ProcessorCoreCount=8)
    @task(bam_file = FILE_IN, output_prefix = IN, db_dir = IN, returns = int)
    def MethylationCaller(self, bam_file, output_prefix, db_dir):
        """
        Takes the merged and sorted bam file and calls the methylation sites.
        Generates a wig file of the potential sites.
        This is performed by running the external program rather than reimplementing
        the code from the main function to make it easier when it comes to updating
        the changes in BS-Seeker2
        """
        import subprocess
        
        args = ["python", "bs_seeker2-call_methylation.py", "--input", str(bam_file), "--output-prefix", str(output_prefix), "--db", db_dir]
        p = subprocess.Popen(args)
        p.wait()
        
        return 1

    
    def clean_up(self, data_dir):
        """
        Clears up the tmp folders
        """
        os.chdir(self.data_dir)
        
        try:
            shutil.rmtree("tmp")
        except:
            pass
        

if __name__ == "__main__":
    import sys
    import os
    from pycompss.api.api import compss_wait_on
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="Parse WGBS data")
    parser.add_argument("--genome", help="Genome name") #             default="GCA_000001405.22")
    parser.add_argument("--srr_id", help="SRR ID of the dataset") #   default="SRR1536575")
    parser.add_argument("--aligner", help="Aligner to use (eg bowtie2)") #   default="bowtie2")
    parser.add_argument("--tmp_dir", help="Temporary data dir")
    parser.add_argument("--data_dir", help="Data directory; location to download SRA FASTQ files and save results")
    parser.add_argument("--aligner_dir", help="Directory for the aligner program")
    parser.add_argument("--local", help="Directory and data files already available", default=0)

    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    srr_id   = args.srr_id
    genome   = args.genome
    aligner  = args.aligner
    aligner_dir = args.aligner_dir
    data_dir = args.data_dir + '/' + srr_id
    tmp_dir  = args.tmp_dir
    local = args.local
    
    genome_dir = args.data_dir + genome + "/chroms/"
    
    start = time.time()
    
    db_dir = ""
    
    pwgbs = process_wgbs()
    
    try:
        os.makedirs(data_dir)
        
        data_dir += "/"
    except:
        pass
    
    
    cf = common()
    
    # Optain the paired FastQ files
    if (local == 1):
        in_files = cf.getFastqFiles(srr_id, data_dir)
    else:
        in_files = [f for f in os.listdir(data_dir + '/' + project) if re.match(run_id, f)]
    
    in_file1 = in_files[0]
    in_file2 = in_files[1]
    out_file1 = in_file1.replace(".fastq", "_filtered.fastq")
    out_file2 = in_file2.replace(".fastq", "_filtered.fastq")
    
    
    
    # Run the FilterReads.py steps for the individual FastQ files
    x = []
    for l in [[in_file1, out_file1], [in_file2, out_file2]]:
        x.append(pwgbs.FilterFastQReads(l))
    x = compss_wait_on(x)
    
    # Run the bs_seeker2-builder.py steps
    pwgbs.Builder(genome_dir + genome + ".fa", "bowtie2", aligner_dir, genome_dir)
        
    # Split the paired fastq files
    tmp_fastq = pwgbs.Splitter(in_file1, in_file2, tag)
    
    bam_sort_files = []
    bam_merge_files = []
    fastq_for_alignment = []
    for bams in tmp_fastq:
        tmp = bams
        tmp.append(aligner, aligner_path, genome_fasta_files)
        fastq_for_alignment.append(tmp)
        bam_root = bams[0] + "_bspe.bam"
        bam_sort_files.append(bams[0] + "_bspe.bam")
        bam_merge_files.append(["-o", bam_root + ".sorted.bam", "-T", bam_root + ".bam_sort", bam_root])
    
    # Run the bs_seeker2-align.py steps on the split up fastq files
    x = []
    for ffa in fastq_for_alignment:
        x.append(pwgbs.Aligner(ffa))
    x = compss_wait_on(x)
    
    # Sort and merge the aligned bam files
    # Pre-sort the original input bam files
    x = map(pysam.sort, bam_sort_files)
    
    f_bam = in_file1.split("/")
    f_bam[-1] = f_bam[-1].replace(".fastq", ".sorted.bam")
    out_bam_file = "/".join(f_bam)
    
    pysam.merge(out_bam_file, bam_merge_files)
    
    pysam.sort("-o", out_bam_file + '.sorted.bam', "-T", out_bam_file + ".bam_sort", out_bam_file)
    
    pysam.index(out_bam_file)
    
    # Run the bs_seeker2-call_methylation.py steps
    x = pwgbs.MethylationCaller(out_bam_file, tag, db_dir)
    x = compss_wait_on(x)
    
    # Tidy up
    pwgbs.clean_up(data_dir)

