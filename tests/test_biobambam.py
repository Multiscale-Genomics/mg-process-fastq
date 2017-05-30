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

import pytest
import random
import os.path

from tool import biobambam_filter

def test_biobambam():
    bbb = biobambam_filter.biobambam()
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    
    bbb.biobambam_filter_alignments(
        resource_path + "fastQForSelRegion.bam",
        resource_path +"_output_test.bam",resource_path
    )