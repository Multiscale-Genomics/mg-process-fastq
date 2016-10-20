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

import argparse, urllib2, gzip, shutil, shlex, subprocess, os, json

from common import common

from pycompss.api.task import task
from pycompss.api.parameter import *


try :
    import pysam
except ImportError :
    print "[Error] Cannot import \"pysam\" package. Have you installed it?"
    exit(-1)

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


class process_rnaseq:
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
    
    
    #@task(data_dir = IN, project_id = IN, run_ids = IN, returns int)
    def inps_peak_calling(self, data_dir, project_id, run_ids):
        """
        Convert Bam to Bed then make Nucleosome peak calls. These are saved as
        bed files That can then get displayed on genome browsers.
        """
        
        with cd('../../lib/iNPS'):
            for run_id in run_ids:
                bam_file = data_dir + project_id + '/' + run_id + '.bam'
                bed_file = data_dir + project_id + '/' + run_id + '.bed'
                bed_out_folder = data_dir + project_id + '/' + run_id + '.inp'
                
                command_line_1 = 'bedtools bamtobed -i ' + bam_file
                command_line_2 = 'python iNPS_V1.2.2.py -i ' + bed_file + ' -o ' + bed_out_folder
                
                
                args = shlex.split(command_line_1)
                with open(bed_file, "w") as f_out:
                    p = subprocess.Popen(args, stdout=f_out)
                    p.wait()
                
                args = shlex.split(command_line_2)
                p = subprocess.Popen(args)
                p.wait()
    
    @task(data_dir = IN, expt = IN, genome_fa = IN, returns = int)
    def main(self, data_dir, expt, genome_fa):
        """
        Main loop
        """
        cf = common()
        
        local_files = data_dir + expt["project_id"]
        
        # Optain the FastQ files
        run_ids = []
        run_fastq_files = []
        
        run_ids = []
        run_fastq_files = {}
        for run_id in expt["run_ids"]:
            run_ids.append(run_id)
            if (expt.has_key("local") == False):
                in_files = cf.getFastqFiles(run_id, data_dir)
            else:
                in_files = [f for f in os.listdir(local_files) if re.match(run_id, f)]
            run_fastq_files[run_id] = in_files
        
        # Run BWA
        paired = 0
        for run_id in expt["run_ids"]:
            if len(run_fastq_files[run_id]) > 1:
                paired = 1
                for i in range(1,len(run_fastq_files[run_id])+1):
                    cf.bwa_align_reads(genome_fa, data_dir, expt["project_id"], run_id + "_" + str(i))
            else:
                cf.bwa_align_reads(genome_fa, data_dir, expt["project_id"], run_id)
        
        self.inps_peak_calling(data_dir, expt["project_id"], expt["run_ids"])
        
        return 1
        
        """
        if expt.has_key("bgd_ids"):
            # Obtain background FastQ files
            bgd_ids = []
            bgd_fastq_files = []
            for bgd_id in expt["bgd_ids"]:
                bgd_ids.append(run_id)
                in_files = cf.getFastqFiles(bgd_id, data_dir)
                bgd_fastq_files[run_id] = in_files
            
            for bgd_id in expt["bgd_ids"]:
                if len(run_fastq_files[run_id]) > 1:
                    paired = 1
                    for i in range(1,len(run_fastq_files[run_id])+1):
                        cf.bwa_align_reads(genome_fa, data_dir, expt["project_id"], bgd_id + "_" + str(i))
                else:
                    cf.bwa_align_reads(genome_fa, data_dir, expt["project_id"], bgd_id)
            
            # DANPOS to call peaks
            self.danpos2_peak_calling(data_dir, expt["project_id"], expt["run_ids"], paired, expt["bgd_ids"][0])
        else:
            # DANPOS to call peaks
            self.danpos2_peak_calling(data_dir, expt["project_id"], expt["run_ids"], paired, None)
        """
    
if __name__ == "__main__":
    import sys
    import os
    
    from pycompss.api.api import compss_wait_on
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="Mnase-seq peak calling")
    parser.add_argument("--species", help="Species (Homo_sapiens)")
    parser.add_argument("--assembly", help="Assembly (GRCm38)")
    parser.add_argument("--project_id", help="Project ID of the dataset (PRJDA47577)")
    parser.add_argument("--run_ids", help="File with list of the experiment run IDs of the dataset")
    parser.add_argument("--data_dir", help="Data directory; location to download ERR FASTQ files and save results")

    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    project     = args.project_id
    run_id_file = args.run_ids
    species     = args.species
    assembly    = args.assembly
    data_dir    = args.data_dir
    
    cf = common()
    prs = process_rnaseq()
    
    try:
        os.makedirs(data_dir)
    except:
        pass
    
    if data_dir[-1] == "/":
        data_dir += "/"

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
    
    # Run main loop
    with open(run_id_file) as data_file:    
        job_id_sets = json.load(data_file)
    
    x = []
    for expt in job_id_sets["expts"]:
        x.append(prs.main(data_dir, expt, genome_fa))
    x = compss_wait_on(x)
    
