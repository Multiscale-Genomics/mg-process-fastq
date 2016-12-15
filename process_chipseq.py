#!/usr/bin/python

"""
.. Copyright 2016 EMBL-European Bioinformatics Institute
 
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

from .. import Tool, Workflow, Metadata
from common import common
from dmp import dmp
import os

try
    from pycompss.api.parameter import FILE_IN, FILE_OUT
    from pycompss.api.task import task
    from pycompss.api.constraint import constraint
except ImportError :
    print "[Warning] Cannot import \"pycompss\" API packages."
    print "          Using mock decorators."
    
    from dummy_pycompss import *

try :
    import pysam
except ImportError :
    print "[Error] Cannot import \"pysam\" package. Have you installed it?"
    exit(-1)


class biobambamTool(Tool):
    """
    Tool to sort and filter bam files
    """
    
    @task(bam_file_in = FILE_IN, bam_file_out = FILE_OUT, tmp_dir = IN)
    def biobambam_filter_alignments(self, bam_file_in, bam_file_out, tmp_dir):
        """
        Sorts and filters the bam file.
        
        It is important that all duplicate alignments have been removed. This
        can be run as an intermediate step, but should always be run as the 
        """
        command_line = 'bamsormadup --tmpfile=' + tmp_dir
        args = shlex.split(command_line)
        with open(bam_file_in, "r") as f_in:
            with open(bam_file_out, "w") as f_out:
                p = subprocess.Popen(args, stdin=f_in, stdout=f_out)
                p.wait()
        
        return True
    
    
    def run(self, input_files, metadata):
        """
        Standard function to call a task
        """
        output_file = input_files[0].replace('.bam', '.filtered.bam')
        
        # handle error
        if not self.biobambam_filter_alignments(input_files[0], output_file):
            output_metadata.set_exception(
                Exception(
                    "biobambamTool: Could not process files {}, {}.".format(*input_files)))
output_file = None
        return ([output_file], [output_metadata])
    

class macs2Tool(Tool):
    """
    Tool for peak calling for ChIP-seq data
    """
    
    @task(name = IN, bam_file = FILE_IN, bam_file_bg = FILE_IN, peak_bed = FILE_OUT, summit_bed = FILE_OUT, narrowPeak = FILE_OUT, broadPeak = FILE_OUT, gappedPeak = FILE_OUT)
    def macs2_peak_calling(self, name, bam_file, bam_file_bg):
        """
        Function to run MACS2 for peak calling
        
        background might need to be optional.
        """
        command_line = 'macs2 callpeak -t ' + bam_file + ' -n ' + name + ' -c ' + bam_file_bg + ' --outdir ' + data_dir + project_id
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()
        
        peak_bed    = name + "_peaks.bed"
        summits_bed = name + "_summits.bed"
        narrowPeak  = name + "_narrowPeak"
        broadPeak   = name + "_broadPeak"
        gappedPeak  = name + "_gappedPeak"
        
        return True
    
    
    def run(self, input_files, metadata):
        """
        Standard function to call a task
        """
        
        # handle error
        if not self.macs2_peak_calling("TODO_Name", input_files[0], input_files[1]):
            output_metadata.set_exception(
                Exception(
                    "macs2_peak_calling: Could not process files {}, {}.".format(*input_files)))
output_file = None
        return ([output_file], [output_metadata])


class process_chipseq(Workflow):
    """
    Functions for downloading and processing Chip-seq FastQ files. Files are
    downloaded from the European Nucleotide Archive (ENA), then aligned,
    filtered and analysed for peak calling
    """
    
    def run(self, file_ids, metadata):
        """
        Main run function
        """
        
        # TODO - Handle multiple file and background files
        genome_fa = file_ids[0]
        file_loc = file_ids[1]
        file_bgd_loc = file_ids[2]
        
        cf = common()
        
        cf.bwa_align_reads(genome_fa["unzipped"], file_loc)
        out_bam = file_loc.replace('.fastq', '.bam')
        
        cf.bwa_align_reads(genome_fa["unzipped"], file_bgd_loc)
        out_bgd_bam = file_bgd_loc.replace('.fastq', '.bam')
        
        # TODO - Multiple files need merging into a single bam file
        
        # Filter the bams
        b3f = biobambamTool(self.configuration)
        b3f_file_out = b3f.run((out_bam), ())
        b3f_file_bgd_out = b3f.run((out_bgd_bam), ())
        
        # MACS2 to call peaks
        macs2 = macs2Tool(self.configuration)
        peak_bed, summits_bed, narrowPeak, broadPeak, gappedPeak = macs2.run((b3f_file_out,  b3f_file_bgd_out), ())
        
        return (b3f_file_out, b3f_file_bgd_out, peak_bed, summits_bed, narrowPeak, broadPeak, gappedPeak)
        
    
if __name__ == "__main__":
    import sys
    import os
    
    from pycompss.api.api import compss_wait_on
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="ChIP-seq peak calling")
    parser.add_argument("--species", help="Species (homo_sapiens)")
    parser.add_argument("--genome", help="Genome FASTA file")
    parser.add_argument("--file", help="Project ID of the dataset")
    parser.add_argument("--bgd_file", help="Project ID of the dataset")
    
    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    species     = args.species
    genome_fa   = args.genome
    file_loc    = args.data_dir
    file_bg_loc = args.tmp_dir
    
    pcs = process_chipseq()
    cf = common()
    
    #
    # MuG Tool Steps
    # --------------
    # 
    # 1. Create data files
    
    # Get the assembly
    genome_fa = cf.getGenomeFromENA(data_dir, species, assembly, False)
    
    #2. Register the data with the DMP
    from dmp import dmp
    
    da = dmp()
    
    print da.get_files_by_user("test")
    
    genome_file = da.set_file("test", genome_fa, "fasta", "Assembly", 9606, None)
    file_in = da.set_file("test", file_loc, "fasta", "ChIP-seq", 9606, None)
    file_bg_in = da.set_file("test", file_bg_loc, "fasta", "ChIP-seq", 9606, None)
    
    print da.get_files_by_user("test")
    
    # 3. Instantiate and launch the App
    from nnn import WorkflowApp
    app = WoekflowApp()
    results = app.launch(process_chipseq, [genome_file, file_in, file_bg_in], {})
    
    print da.get_files_by_user("test")
    
