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

class indexerTool(Tool):
    """
    Tool for running indexers over a genome FASTA file
    """
    
    @task()
    def bowtie2_indexer():
        """
        Bowtie2 Indexer
        """
    
    
    @task()
    def bwa_indexer():
        """
        BWA Indexer
        """
    
    
    @task()
    def gem_indexer():
        """
        GEM Indexer
        """


class processs_genome(Workflow):
    """
    Workflow to download and pre-index a given genome
    
    The downloading can be done using the current common.py functions. These
    should be prevented from running the indexing step as this will be done as
    part of this workflow.
    """
    
