#!usr/bin/python

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

from tool import bs_seeker_aligner

def test_bs_seeker_aligner():    
    
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genomefa_file = resource_path + "bsSeeker.Mouse.GRCm38.fasta"
    fastq1_file = resource_path + "bsSeeker.Mouse.GRCm38_1.fastq"
    fastq2_file = resource_path + "bsSeeker.Mouse.GRCm38_2.fastq"
    bt2_1file = resource_path + "bsSeeker.Mouse.GRCm38_2.fasta_bowtie2/C_C2T.1.bt2"
    bt2_2file = resource_path + "bsSeeker.Mouse.GRCm38_2.fasta_bowtie2/C_C2T.2.bt2"
    bt2_3file = resource_path + "bsSeeker.Mouse.GRCm38_2.fasta_bowtie2/C_C2T.3.bt2"
    bt2_4file = resource_path + "bsSeeker.Mouse.GRCm38_2.fasta_bowtie2/C_C2T.4.bt2"
    bt2_rev_1file = resource_path + "bsSeeker.Mouse.GRCm38_2.fasta_bowtie2/C_C2T.rev.1.bt2"
    bt2_rev_2file = resource_path + "bsSeeker.Mouse.GRCm38_2.fasta_bowtie2/C_C2T.rev.2.bt2"
    out_file = resource_path + "bsSeeker.Mouse.GRCm38.bam"
    
    bsa = bs_seeker_aligner.bssAlignerTool()
    bsa.run(
        [genomefa_file, 
         fastq1_file, 
         fastq2_file, 
         out_file, 
         bt2_1file, 
         bt2_2file, 
         bt2_3file, 
         bt2_4file, 
         bt2_rev_1file, 
         bt2_rev_2file ],
            {"aligner":"bowtie2", "aligner_path":"/Users/reham/lib/bowtie2-2.3.2","bss_path":"/Users/reham/lib/BSseeker2-2.1.2Beta"})
    
    
test_bs_seeker_aligner()


   