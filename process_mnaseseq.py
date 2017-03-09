#!/usr/bin/python

"""
.. Copyright 2017 EMBL-European Bioinformatics Institute
 
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

from basic_modules import Tool, Workflow, Metadata

from functools import wraps

from tool.common import common
from tool.common import cd

import tool

try :
    from pycompss.api.parameter import *
    from pycompss.api.task import task
except ImportError :
    print "[Warning] Cannot import \"pycompss\" API packages."
    print "          Using mock decorators."
    
    from dummy_pycompss import *


try :
    import pysam
except ImportError :
    print "[Error] Cannot import \"pysam\" package. Have you installed it?"
    exit(-1)

# ------------------------------------------------------------------------------

class process_mnaseseq(Workflow):
    """
    Functions for downloading and processing Mnase-seq FastQ files. Files are
    downloaded from the European Nucleotide Archive (ENA), then aligned,
    filtered and analysed for peak calling
    """
    
    def run(self, file_ids, metadata):
        """
        Main run function for processing MNase-Seq FastQ data. Pipeline aligns
        the FASTQ files to the genome using BWA. iNPS is then used for peak
        calling to identify nucleosome position sites within the genome.
                
        Parameters
        ----------
        files_ids : list
            List of file locations
        metadata : list
            Required meta data
        
        Returns
        -------
        outputfiles : list
            List of locations for the output bam, bed and tsv files
        """
        
        cwa_align = tool.bwaAlignerTool(self.configuration)
        out_bam, out_meta = cwa_align.run((files_ids[0], files_ids[1]), ())
        
        # Needs moving to its own tool
        inps = tool.inps(self.configuration)
        out_peak_bed, out_meta = inps.inps_peak_calling(out_bam, ())
        
        return (out_peak_bed[0])

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    import os
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="Mnase-seq peak calling")
    parser.add_argument("--species", help="Species (Mus_musculus)")
    parser.add_argument("--genome", help="Genome assembly FASTA file")
    parser.add_argument("--file", help="Location of FASTQ file")
    #parser.add_argument("--project_id", help="Project ID of the dataset (PRJDA47577)")
    #parser.add_argument("--run_ids", help="File with list of the experiment run IDs of the dataset")
    #parser.add_argument("--data_dir", help="Data directory; location to download ERR FASTQ files and save results")

    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    species    = args.species
    genome     = args.genome
    fastq_file = args.file
    
    #cf = common()
    ps = process_mnaseseq()
    
    #try:
    #    os.makedirs(data_dir)
    #except:
    #    pass
    #
    #if data_dir[-1] != "/":
    #    data_dir += "/"
    #
    #try:
    #    os.makedirs(data_dir + project)
    #except:
    #    pass
    #
    #try:
    #    os.makedirs(data_dir + species + "_" + genome)
    #except:
    #    pass
    
    # Get the assembly
    #genome_fa = cf.getGenomeFromENA(data_dir, species, genome, False)
    
    # Run main loop
    #with open(run_id_file) as data_file:    
    #    job_id_sets = json.load(data_file)
    #
    #for expt in job_id_sets["expts"]:
    #    local_files = data_dir + expt["project_id"]
    #    
    #    # Optain the FastQ files
    #    run_ids = []
    #    run_fastq_files = []
    #    
    #    run_ids = []
    #    run_fastq_files = {}
    #    for run_id in expt["run_ids"]:
    #        run_ids.append(run_id)
    #        if (expt.has_key("local") == False):
    #            in_files = cf.getFastqFiles(expt["project_id"], data_dir, run_id)
    #        else:
    #            in_files = [f for f in os.listdir(local_files) if re.match(run_id, f)]
    #        run_fastq_files[run_id] = in_files
    #    
    #    peak_files = []
    #    for run_id in run_fastq_files:
    #        peak_file = ps.run((genome_fa["unzipped"], run_fastq_files[run_id][0]), ())
    #        peak_files.append(peak_file)
    
    peak_file = ps.run((genome_fa, fastq_file), ())
    
    print peak_files

