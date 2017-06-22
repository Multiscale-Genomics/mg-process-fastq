#!usr/bin/env python3

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

import pytest
import random
import os
import sys
import pysam

from tool import bs_seeker_methylation_caller

def test_bs_seeker_methylation_caller():    
    
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genomefa_file = resource_path + "bsSeeker.Mouse.GRCm38.fasta_bowtie2"    
    bam_file = resource_path + "bsSeeker.Mouse.GRCm38.bam"

    pysam.sort("-o", str(bam_file), str(bam_file))

    home = os.path.expanduser('~')
    
    bsmc = bs_seeker_methylation_caller.bssMethylationCallerTool()
    bsmc.run(
        [genomefa_file, bam_file ],
        {'bss_path' : home + "/bin"}
    )
   