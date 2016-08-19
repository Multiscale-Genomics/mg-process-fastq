#!/usr/bin/python
# -*- coding: utf-8 -*-
'''process whole genome bisulfate sequencing FastQ files'''
import argparse, time
from pycompss.api.task import task
from pycompss.api.parameter import *


def getFastqFiles(self, sra_id, data_dir):
    f_index = urllib2.urlopen(
    'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=' + str(self.sra_id) + '&result=read_run&fields=study_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp&download=txt')
    data = f_index.read()
    rows = data.split("\n")
    row_count = 0
    files = []
    for row in rows:
        if row_count == 0:
            row_count += 1
            continue
        
        row = row.rstrip()
        row = row.split("\t")
        
        if len(row) == 0:
            continue
        
        project = row[0]
        srr_ir = row[1]
        fastq_files = row[6].split(';')
        row_count += 1
        
        for fastq_file in fastq_files:
            file_name = fastq_file.split("/")
            print fastq_file
            print data_dir + file_name[-1]
            print file_name[-1]
            files.append(data_dir + file_name[-1])
            dl_file = open(data_dir + file_name[-1], "w")
            dl = urllib2.urlopen("ftp://" + fastq_file)
            dl_file.write(dl.read())
            dl_file.close()
    
    return files


@task(infile = FILE_IN, outfile = FILE_OUT)
def FilterFastQReads(infile, outfile):
    """
    This is optional, but removes reads that can be problematic for the
    alignment of whole genome datasets.
    
    If performing RRBS then this step can be skipped
    """
    from FilterReads import *
    FilterReads(infile, outfile, True)


@constraint(ProcessorCoreCount=8)
@task()
def Builder(fasta_file, aligner, aligner_path, ref_path):
    """
    Function to submit the FASTA file for the reference sequence and build the
    required index file used by the aligner.
    """
    from bs_index.wg_build import *
    
    builder_exec = os.path.join(aligner_path or aligner_path[aligner],
                                {BOWTIE   : 'bowtie-build',
                                 BOWTIE2  : 'bowtie2-build',
                                 SOAP     : '2bwt-builder',
                                 RMAP     : '' # do nothing
                                }[options.aligner])

    build_command = builder_exec + { BOWTIE   : ' -f %(fname)s.fa %(fname)s',
                                     BOWTIE2  : ' -f %(fname)s.fa %(fname)s',
                                     SOAP     : ' %(fname)s.fa'
                                   }[aligner]
    
    wg_build(fasta_file, build_command, ref_path, aligner)
    

def Splitter(in_file1, in_file2, tag):
    """
    Function to divide the FastQ files into separte sub files of 1000000
    sequences so that the aligner can get run in parallel.
    
    Returns: Returns a list of lists of the files that have been generated.
             Each sub list containing the 2 paired end files for that subset.
    """
    from fastqreader import *
    
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
@task() # No. of CPU, RAM required
def Aligner(input_fastq1, input_fastq2, aligner, aligner_path, genome_fasta):
    """
    Alignment of the paired ends to the reference genome
    Generates bam files for the alignments
    This is performed by running the external program rather than reimplementing
    the code from teh main function to make it easier when it comes to updating
    the changes in BS-Seeker2
    """
    import subprocess
    
    args = ["python", "bs_seeker2-align.py", "--input_1", input_fastq1, "--input_2", input_fastq2, "--aligner", aligner, "--path", aligner_path_path, "--genome", genome_fasta, "--bt2-p", "4"]
    p = subprocess.Popen(args)
    p.wait()


@task() # No. of CPUs, RAM required
def MethylationCaller(bam_file, output_prefix, db_dir):
    """
    Takes the merged and sorted bam file and calls the methylation sites.
    Generates a wig file of the potential sites.
    This is performed by running the external program rather than reimplementing
    the code from teh main function to make it easier when it comes to updating
    the changes in BS-Seeker2
    """
    import subprocess
    
    args = ["python", "bs_seeker2-call_methylation.py", "--input", str(bam_file), "--output-prefix", str(output_prefix), "--db", db_dir]
    p = subprocess.Popen(args)
    p.wait()


if __name__ == "__main__":
    import sys
    import os
    from pycompss.api.api import compss_wait_on
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="Load adjacency list into HDF5 file")
    parser.add_argument("--genome", help="Genome name") #             default="GCA_000001405.22")
    parser.add_argument("--dataset", help="Name of the dataset") #    default="SRR1536575")
    parser.add_argument("--tmp_dir", help="Temporary data dir")
    parser.add_argument("--data_dir", help="Data directory; location to download SRA FASTQ files and save results")

    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    sra_id   = args.dataset
    genome   = args.genome
    data_dir = args.data_dir
    tmp_dir  = args.tmp_dir
    
    start = time.time()
    
    db_dir = ""
    
    # Optain the paired FastQ files
    in_files = getFastqFiles(sra_id, data_dir)
    in_file1 = in_files[0]
    in_file2 = in_files[1]
    out_file1 = in_file1.replace(".fastq", "_filtered.fastq")
    out_file2 = in_file2.replace(".fastq", "_filtered.fastq")
    
    # Run the FilterReads.py steps for the individual FastQ files
    x = map(FilterFastQReads, [[in_file1, out_file1], [in_file2, out_file2]])
    x = compss_wait_on(x)
    
    # Run the bs_seeker2-builder.py steps
    x = Builder(fasta_file, aligner, aligner_path, ref_path)
    x = compss_wait_on(x)
    
    # Split the paired fastq files
    tmp_fastq = Splitter(in_file1, in_file2, tag)
    
    bam_sort_files = []
    bam_merge_files = []
    fastq_for_alignment = []
    for bams in tmp_fastq:
        tmp = bams
        tmp.append(aligner, aligner_path, genome_fasta_files)
        fastq_for_alignment.append(tmp)
        bam_root = bams[0] + "_bspe.bam"
        bam_sort_files.append(bams[0] + "_bspe.bam")
        bam_sort_files.append(["-o", bam_root + ".sorted.bam", "-T", bam_root + ".bam_sort", bam_root])
    
    # Run the bs_seeker2-align.py steps on the split up fastq files
    x = map(Aligner, fastq_for_alignment)
    x = compss_wait_on(x)
    
    # Sort and merge the aligned bam files
    # Pre-sort the original input bam files
    x = map(pysam.sort, bam_sort_files)
    x = compss_wait_on(x)
    
    f_bam = in_file1.split("/")
    f_bam[-1] = f_bam[-1].replace(".fastq", ".sorted.bam")
    out_bam_file = "/".join(f_bam)
    
    pysam.merge(out_bam_file, bam_merge_files)
    
    pysam.sort("-o", out_bam_file + '.sorted.bam', "-T", out_bam_file + ".bam_sort", out_bam_file)
    
    pysam.index(out_bam_file)
    
    # Run the bs_seeker2-call_methylation.py steps
    x = MethylationCaller(out_bam_file, tag, db_dir)
    x = compss_wait_on(x)
    
    # Tidy up
    
