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

import argparse, urllib2, gzip, shutil, shlex, subprocess, os.path, json

from basic_modules import Tool, Workflow, Metadata

from functools import wraps

from dmp import dmp

from tool import bwa_aligner
from tool import biobambam_filter
from tool import macs2

import os

try:
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

# ------------------------------------------------------------------------------

class process_chipseq(Workflow):
    """
    Functions for processing Chip-Seq FastQ files. Files are the aligned,
    filtered and analysed for peak calling
    """
    
    def run(self, file_ids, metadata):
        """
        Main run function for processing ChIP-seq FastQ data. Pipeline aligns
        the FASTQ files to the genome using BWA. MACS 2 is then used for peak
        calling to identify transcription factor binding sites within the
        genome.
                
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
        
        # TODO - Handle multiple file and background files
        genome_fa = file_ids[0]
        bwa_amb = file_ids[1]
        bwa_ann = file_ids[2]
        bwa_bwt = file_ids[3]
        bwa_pac = file_ids[4]
        bwa_sa  = file_ids[5]
        file_loc = file_ids[6]
        file_bgd_loc = file_ids[7]
        

        
        out_bam = file_loc.replace('.fastq', '.bam')
        
        bwa = bwa_aligner.bwaAlignerTool()
        out_bam = file_loc.replace('.fastq', '.bam')
        out_file_bam, out_bam_meta = bwa.run(
            [genome_fa, file_loc, out_bam, bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_sa],
            {}
        )
        
        if file_bgd_loc != None:
            out_bgd_bam = file_bgd_loc.replace('.fastq', '.bam')
            out_bgd_bam_res, out_bgd_bam_meta = bwa.run(
                [genome_fa, file_bgd_loc, out_bgd_bam, bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_sa],
                {}
            )
        
        # TODO - Multiple files need merging into a single bam file
        
        # Filter the bams
        b3f = biobambam_filter.biobambam()
        b3f_file_out = file_loc.replace('.fastq', '.filtered.bam')
        b3f_out = b3f.run([out_bam, b3f_file_out], {})
        
        if file_bgd_loc != None:
            b3f_bgd_file_out = file_bgd_loc.replace('.fastq', '.filtered.bam')
            b3f_bgd_file_out = b3f.run([out_bgd_bam, b3f_bgd_file_out], {})
        else:
            b3f_bgd_file_out = [None]
        
        # MACS2 to call peaks
        m = macs2.macs2()
        if file_bgd_loc != None:
            peak_bed, summits_bed, narrowPeak, broadPeak, gappedPeak = m.run([b3f_file_out[0],  b3f_bgd_file_out[0]], {})
        else:
            peak_bed, summits_bed, narrowPeak, broadPeak, gappedPeak = m.run([b3f_file_out[0]], {})
        
        return (b3f_file_out, b3f_file_bgd_out, peak_bed, summits_bed, narrowPeak, broadPeak, gappedPeak)

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    import os
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="ChIP-seq peak calling")
    parser.add_argument("--taxon_id", help="Taxon_ID (9606)")
    parser.add_argument("--genome", help="Genome FASTA file")
    parser.add_argument("--assembly", help="Genome assembly ID (GCA_000001405.25)")
    parser.add_argument("--file", help="Location of FASTQ input file")
    parser.add_argument("--bgd_file", help="Location of FASTQ background file", default=None)
    
    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    taxon_id    = args.taxon_id
    genome_fa   = args.genome
    assembly    = args.assembly
    file_loc    = args.file
    file_bg_loc = args.bgd_file
    
    #
    # MuG Tool Steps
    # --------------
    # 
    # 1. Create data files
    
    # Get the assembly
    
    #2. Register the data with the DMP
    da = dmp(test=True)
    
    print da.get_files_by_user("test")
    
    genome_file = da.set_file("test", genome_fa, "fasta", "Assembly", taxon_id, meta_data={'assembly' : assembly})
    genome_file_idx1 = da.set_file("test", genome_fa + ".amb", "fasta", "Assembly", taxon_id, meta_data={'assembly' : assembly})
    genome_file_idx2 = da.set_file("test", genome_fa + ".ann", "fasta", "Assembly", taxon_id, meta_data={'assembly' : assembly})
    genome_file_idx3 = da.set_file("test", genome_fa + ".bwt", "fasta", "Assembly", taxon_id, meta_data={'assembly' : assembly})
    genome_file_idx4 = da.set_file("test", genome_fa + ".pac", "fasta", "Assembly", taxon_id, meta_data={'assembly' : assembly})
    genome_file_idx5 = da.set_file("test", genome_fa + ".sa", "fasta", "Assembly", taxon_id, meta_data={'assembly' : assembly})
    file_in = da.set_file("test", file_loc, "fasta", "ChIP-seq", taxon_id, meta_data={'assembly' : assembly})
    file_bg_in = da.set_file("test", file_bg_loc, "fasta", "ChIP-seq", taxon_id, meta_data={'assembly' : assembly})
    
    print da.get_files_by_user("test")

    files = [
        genome_fa,
        genome_fa + ".amb",
        genome_fa + ".ann",
        genome_fa + ".bwt",
        genome_fa + ".pac",
        genome_fa + ".sa",
        file_loc,
        file_bg_loc
    ]
    
    # 3. Instantiate and launch the App
    #from basic_modules import WorkflowApp
    #app = WorkflowApp()
    #results = app.launch(process_chipseq, [genome_file, file_in, file_bg_in], {'user_id' : 'test'})

    pc = process_chipseq()
    results = pc.run(files, {'user_id' : 'test'})
    
    print da.get_files_by_user("test")
    
