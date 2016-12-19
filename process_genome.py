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


class processs_genome(Workflow):
    """
    Workflow to download and pre-index a given genome
    
    The downloading can be done using the current common.py functions. These
    should be prevented from running the indexing step as this will be done as
    part of this workflow.
    """
    
    def run(self, file_ids, metadata):
        """
        Main run function
        """
        
        genome_fa = file_ids[0]
        
        # Bowtie2 Indexer
        bt = bowtieIndexerTool(self.configuration)
        bti = bt.run((genome_fa), ())
        
        # BWA Indexer
        bwa = bwaIndexerTool(self.configuration)
        bwai = bwa.run((genome_fa), ())
        
        # GEM Indexer
        gem = gemIndexerTool(self.configuration)
        gemi = gem.run((genome_fa), ())


if __name__ == "__main__":
    import sys
    import os
    
    from pycompss.api.api import compss_wait_on
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="Index the genome file")
    parser.add_argument("--species", help="Species (9606)")
    parser.add_argument("--genome", help="Genome FASTA file")
    
    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    species     = args.species
    file_loc    = args.data_dir
    file_bg_loc = args.tmp_dir
    
    pcs = process_genome()
    
    #
    # MuG Tool Steps
    # --------------
    # 
    # 1. Create data files
    
    # Get the assembly
    genome_fa   = args.genome
    #genome_fa = cf.getGenomeFromENA(data_dir, species, assembly, False)
    
    #2. Register the data with the DMP
    from dmp import dmp
    
    da = dmp()
    
    print da.get_files_by_user("test")
    
    genome_file = da.set_file("test", genome_fa, "fasta", "Assembly", species, None)
    
    print da.get_files_by_user("test")
    
    # 3. Instantiate and launch the App
    from nnn import WorkflowApp
    app = WorkflowApp()
    results = app.launch(process_genome, [genome_file], {})
    
    print da.get_files_by_user("test")
    
