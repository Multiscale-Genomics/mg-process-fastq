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
    parser.add_argument("--assembly", help="Genome assembly ID (GCA_000001635.2)")
    parser.add_argument("--taxon_id", help="Taxon ID (10090)")
    parser.add_argument("--species", help="Species (Mus_musculus)")
    parser.add_argument("--genome", help="Genome assembly FASTA file")
    parser.add_argument("--file", help="Location of FASTQ file")
    
    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    species    = args.species
    genome     = args.genome
    taxon_id   = args.taxon_id
    assembly   = argsd.assembly
    fastq_file = args.file
    
    ps = process_mnaseseq()
    
    da = dmp()
    
    print da.get_files_by_user("test")
    
    genome_file = da.set_file("test", genome_fa, "fasta", "Assembly", taxon_id, {'assembly' : assembly})
    file_in = da.set_file("test", file_loc, "fastq", "ChIP-seq", taxon_id, {'assembly' : assembly})
    
    print da.get_files_by_user("test")
    
    # 3. Instantiate and launch the App
    from basic_modules import WorkflowApp
    app = WorkflowApp()
    results = app.launch(process_mnaseseq, [genome_fa, fastq_file], {})
    #peak_file = ps.run((genome_fa, fastq_file), ())
    
    print peak_files

