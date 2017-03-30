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

import tool

try :
    from pycompss.api.parameter import *
    from pycompss.api.task import task
except ImportError :
    print "[Warning] Cannot import \"pycompss\" API packages."
    print "          Using mock decorators."
    
    from dummy_pycompss import *

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

        genome_fa = file_ids[0]
        bwa_amb = file_ids[1]
        bwa_ann = file_ids[2]
        bwa_bwt = file_ids[3]
        bwa_pac = file_ids[4]
        bwa_sa  = file_ids[5]
        file_loc = file_ids[6]
        file_bgd_loc = file_ids[7]
        
        cwa_align = tool.bwaAlignerTool(self.configuration)
        out_bam, out_meta = cwa_align.run(
            [genome_fa, file_bgd_loc, out_bgd_bam, bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_sa],
            {}
        )
        
        # Needs moving to its own tool
        inps = tool.inps(self.configuration)
        out_peak_bed, out_meta = inps.inps_peak_calling(out_bam, {})
        
        return (out_peak_bed[0])

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    import os
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="Mnase-seq peak calling")
    parser.add_argument("--assembly", help="Genome assembly ID (GCA_000001635.2)")
    parser.add_argument("--taxon_id", help="Taxon ID (10090)")
    parser.add_argument("--genome", help="Genome assembly FASTA file")
    parser.add_argument("--file", help="Location of FASTQ file")
    
    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    genome_fa  = args.genome
    taxon_id   = args.taxon_id
    assembly   = args.assembly
    fastq_file = args.file
    
    da = dmp(test=True)
    
    print da.get_files_by_user("test")
    
    genome_file = da.set_file("test", genome_fa, "fasta", "Assembly", taxon_id, meta_data={'assembly' : assembly})
    genome_file_idx1 = da.set_file("test", genome_fa + ".amb", "fasta", "Assembly", taxon_id, meta_data={'assembly' : assembly})
    genome_file_idx2 = da.set_file("test", genome_fa + ".ann", "fasta", "Assembly", taxon_id, meta_data={'assembly' : assembly})
    genome_file_idx3 = da.set_file("test", genome_fa + ".bwt", "fasta", "Assembly", taxon_id, meta_data={'assembly' : assembly})
    genome_file_idx4 = da.set_file("test", genome_fa + ".pac", "fasta", "Assembly", taxon_id, meta_data={'assembly' : assembly})
    genome_file_idx5 = da.set_file("test", genome_fa + ".sa", "fasta", "Assembly", taxon_id, meta_data={'assembly' : assembly})
    file_in = da.set_file("test", file_loc, "fastq", "Mnase-seq", taxon_id, meta_data={'assembly' : assembly})
    
    print da.get_files_by_user("test")

    files = [
        genome_fa,
        genome_fa + ".amb",
        genome_fa + ".ann",
        genome_fa + ".bwt",
        genome_fa + ".pac",
        genome_fa + ".sa",
        file_loc
    ]
    
    # 3. Instantiate and launch the App
    #from basic_modules import WorkflowApp
    #app = WorkflowApp()
    #results = app.launch(process_mnaseseq, [genome_fa, fastq_file], 'user_id' : 'test'})
    
    pm = process_mnaseseq()
    peak_file = pm.run(files, {'user_id' : 'test'})
    
    print peak_files

