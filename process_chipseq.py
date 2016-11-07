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
    
    
    @task(bam_file_in = FILE_IN, bam_file_out = FILE_OUT, tmp_dir = IN)
    def biobambam_filter_alignments(self, bam_file_in, bam_file_out, tmp_dir):
        """
        Sorts and filters the bam file.
        
        It is important that all duplicate alignments have been removed. This
        can be run as an intermediate step, but should always be run as the 
        """
        command_line = 'bamsormadup --tmpfile=' +  + '.tmp'
        args = shlex.split(command_line)
        with open(bam_file_in, "r") as f_in:
            with open(bam_file_out, "w") as f_out:
                p = subprocess.Popen(args, stdin=f_in, stdout=f_out)
                p.wait()
    
    
    @task(data_dir = IN, project_id = IN, run_id = IN, background_id = IN)
    def macs2_peak_calling(self, data_dir, project_id, run_id, background_id):
        """
        Function to run MACS2 for peak calling
        
        background might need to be optional.
        """
        bam_file = data_dir + project_id + '/' + run_id + '.filtered.bam'
        bam_background_file = data_dir + project_id + '/' + background_id + '.filtered.bam'
        command_line = 'macs2 callpeak -t ' + bam_file + ' -n ' + run_id + ' -c ' + bam_background_file + ' --outdir ' + data_dir + project_id
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()
        
    
    def main(self, data_dir, tmp_dir, expt, genome_fa):
        """
        Main loop
        """
        
        cf = common()
        
        # Optain the FastQ files
        run_ids = []
        for run_id in expt["run_ids"]:
            run_ids.append(run_id)
            if (expt.has_key("local") == False):
                in_files = cf.getFastqFiles(expt["project_id"], data_dir, run_id)
        
        # Obtain background FastQ files
        bgd_ids = []
        for bgd_id in expt["bgd_ids"]:
            bgd_ids.append(run_id)
            if (expt.has_key("local") == False):
                in_files = cf.getFastqFiles(expt["project_id"], data_dir, bgd_id)
        
        # Run BWA
        for run_id in expt["run_ids"]:
            cf.bwa_align_reads(genome_fa["unzipped"], data_dir, expt["project_id"], run_id)
            out_bam = data_dir + expt["project_id"] + '/' + run_id + '.bam'
        
        for bgd_id in expt["bgd_ids"]:
            cf.bwa_align_reads(genome_fa["unzipped"], data_dir, expt["project_id"], bgd_id)
        
        final_run_id = expt["project_id"] + "_" + expt["group_name"] + "_run"
        final_bgd_id = expt["project_id"] + "_" + expt["group_name"] + "_bgd"
        
        # Merge Bam files
        cf.merge_bam(data_dir, expt["project_id"], final_run_id, expt["run_ids"])
        cf.merge_bam(data_dir, expt["project_id"], final_bgd_id, expt["bgd_ids"])
        
        # Filter the bams
        bam_file_in  = data_dir + project_id + '/' + final_run_id + '.bam'
        bam_file_out = data_dir + project_id + '/' + final_run_id + '.filtered.bam'
        self.biobambam_filter_alignments(bam_file_in, bam_file_out, tmp_dir + expt["project_id"] + '.run.tmp')
        bam_bgd_file_in  = data_dir + project_id + '/' + final_bgd_id + '.bam'
        bam_bgd_file_out = data_dir + project_id + '/' + final_bgd_id + '.filtered.bam'
        self.biobambam_filter_alignments(bam_bgd_file_in, bam_bgd_file_out, tmp_dir + expt["project_id"] + '.bgd.tmp')
        
        # MACS2 to call peaks
        self.macs2_peak_calling(data_dir, expt["project_id"], final_run_id, final_bgd_id)
        
    
if __name__ == "__main__":
    import sys
    import os
    
    from pycompss.api.api import compss_wait_on
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="ChIP-seq peak calling")
    parser.add_argument("--species", help="Species (homo_sapiens)")
    parser.add_argument("--assembly", help="Assembly (GRCh38)")
    parser.add_argument("--project_id", help="Project ID of the dataset")
    parser.add_argument("--run_ids", help="JSON file with list of the experiment run IDs and background data (if available) of the dataset")
    parser.add_argument("--data_dir", help="Data directory; location to download ERR FASTQ files and save results")
    parser.add_argument("--tmp_dir", help="Temporary Data directory")

    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    project     = args.project_id
    run_id_file = args.run_ids
    species     = args.species
    assembly    = args.assembly
    data_dir    = args.data_dir
    tmp_dir    = args.tmp_dir
    
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
        os.makedirs(tmp_dir + project)
    except:
        pass
    
    try:
        os.makedirs(data_dir + species + "_" + assembly)
    except:
        pass
    
    # Get the assembly
    genome_fa = cf.getGenomeFromENA(data_dir, species, assembly)
    
    with open(run_id_file) as data_file:    
        job_id_sets = json.load(data_file)
    
    x = []
    for expt in job_id_sets["expts"]:
        pcs.main(data_dir, tmp_dir, expt, genome_fa)
    
    
    
    
