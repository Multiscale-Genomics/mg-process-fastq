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

import argparse, urllib2, gzip, shutil, shlex, subprocess, os.path


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

class process_rnaseq(Workflow):
    """
    Functions for downloading and processing RNA-seq FastQ files. Files are
    downloaded from the European Nucleotide Archive (ENA), then they are mapped
    to quantify the amount of cDNA
    """
    
    def run(self, file_ids, metadata):
        """
        Main run function for processing RNA-Seq FastQ data. Pipeline aligns
        the FASTQ files to the genome using Kallisto. Kallisto is then also
        used for peak calling to identify levels of expression.
                
        Parameters
        ----------
        files_ids : list
            List of file locations (genome FASTA, FASTQ_01, FASTQ_02 (for
            paired ends))
        metadata : list
            Required meta data
        
        Returns
        -------
        outputfiles : list
            List of locations for the output bam, bed and tsv files
        """
        
        genome_fa = file_ids[0]
        
        # Index the cDNA
        # This could get moved to the general tools section
        ki = tool.kallistoIndexerTool(self.configuration)
        genome_idx_loc = file_loc.replace('.fa', '.idx')
        gi_out = ki.run((genome_fa, genome_id_loc), ())
        
        # Quantification
        kq = tool.kallistoQuantificationTool(self.configuration)
        
        if len(file_ids) == 2:
            results = kq.run((file_ids[0], file_ids[0]), ())
        elif len(file_ids) == 3:
            results = kq.run((file_ids[0], file_ids[1], file_ids[2]), ())
        
        return results[0]


# ------------------------------------------------------------------------------


if __name__ == "__main__":
    import sys
    import os
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="Parse RNA-seq for expression analysis")
    parser.add_argument("--assembly", help="Genome assembly ID (GCA_000001405.25)")
    parser.add_argument("--taxon_id", help="Taxon_ID (9606)")
    parser.add_argument("--genome", help="Location of the genome cDNA FASTA file")
    parser.add_argument("--file", help="Location of the FASTQ file")
    parser.add_argument("--file2", help="[OPTIONAL] Location of the paired end FASTQ file", default=None)

    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    assembly   = args.assembly
    taxon_id    = args.taxon_id
    genome_fa   = args.genome
    file_loc    = args.file
    paired_file = args.file2
    
    da = dmp(test=True)
    
    print da.get_files_by_user("test")
    
    genome_file = da.set_file("test", genome_fa, "fasta", "cDNA", taxon_id, {'assembly' : assembly})
    file_in = da.set_file("test", file_loc, "fasta", "RNA-seq", taxon_id, {'assembly' : assembly})
    
    print da.get_files_by_user("test")
    
    # 3. Instantiate and launch the App
    #from basic_modules import WorkflowApp
    #app = WorkflowApp()
    #results = app.launch(process_rnaseq, [genome_file, file_in], {})

    files = [
        genome_fa,
        file_loc
    ]

    pr = process_rnaseq()
    if paired_file == None:
        resutls = pr.run(files, {'user_id' : 'test'})
    else:
        files.append(paired_file)
        file2_in = da.set_file("test", paired_file, "fasta", "RNA-seq", taxon_id, {'assembly' : assembly})
        resutls = pr.run(files, {'user_id' : 'test'})
    
    print da.get_files_by_user("test")
    
    
