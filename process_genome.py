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

from dmp import dmp

from tool import bowtie_indexer
from tool import bwa_indexer
from tool import gem_indexer

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

class process_genome(Workflow):
    """
    Workflow to download and pre-index a given genome
    """
    
    def __init__(self):
        """
        Initialise the class
        """
    
    def run(self, file_ids, metadata):
        """
        The downloading can be done using the current common.py functions. These
        should be prevented from running the indexing step as this will be done
        as part of this workflow.
        
        
        Parameters
        ----------
        file_ids : list
            List of file locations
        metadata : list
            Required meta data
        
        Returns
        -------
        outputfiles : list
            List of locations for the output bam, bed and tsv files
        """
        
        genome_fa = file_ids[0]
        
        # Bowtie2 Indexer
        bt = bowtie_indexer.bowtieIndexerTool()
        bti, btm = bt.run([genome_fa], metadata)
        
        # BWA Indexer
        bwa = bwa_indexer.bwaIndexerTool()
        bwai, bwam = bwa.run([genome_fa], metadata)
        
        # GEM Indexer
        gem = gem_indexer.gemIndexerTool()
        gemi, gemm = gem.run(([genome_fa], metadata)

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    import os
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="Index the genome file")
    parser.add_argument("--species", help="Species (9606)")
    parser.add_argument("--genome", help="Genome FASTA file")
    parser.add_argument("--assembly", help="Assembly ID")

    
    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    species     = args.species
    
    #
    # MuG Tool Steps
    # --------------
    # 
    # 1. Create data files
    #    This should have already been done by the VRE - Potentially. If these
    #    Are ones that are present in the ENA then I would need to download them
    
    # Get the assembly
    genome_fa   = args.genome
    assembly    = args.assembly
    
    #2. Register the data with the DMP
    da = dmp(test=True)
    
    print da.get_files_by_user("test")
    
    genome_file = da.set_file("test", genome_fa, "fasta", "Assembly", species, meta_data={'assembly' : assembly})
    
    print da.get_files_by_user("test")
    
    # 3. Instantiate and launch the App
    #from basic_modules import WorkflowApp
    #app = WorkflowApp()
    #results = app.launch(process_genome, [genome_file], {})

    pg = process_genome()
    results = pg.run([genome_fa], {'user_id' : 'test', 'meta_data' : {'assembly' : assembly}})

    
    print da.get_files_by_user("test")
    
