#!usr/bin/python

import pytest
import random

from tool import biobambam_filter

def test_biobambam():
    bbb = biobambam_filter.biobambam()
    
    bbb.biobambam_filter_alignments("/Users/reham/Documents/fastQForSelRegion.bam","/Users/reham/Documents/_output_test.bam","/Users/reham/Documents/temp")
    
    
    
test_biobambam()