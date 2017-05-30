#!usr/bin/python

import pytest
import random
import os.path

from tool import biobambam_filter

def test_biobambam():
    bbb = biobambam_filter.biobambam()
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    
    bbb.biobambam_filter_alignments(resource_path+"fastQForSelRegion.bam",resource_path+"_output_test.bam",resource_path)
    
    
    
test_biobambam()