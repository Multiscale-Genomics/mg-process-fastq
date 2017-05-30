#!usr/bin/python

import pytest
import random
import os.path

from tool import macs2

def test_macs2():
    m = macs2.macs2()
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    m.macs2_peak_calling("selectGenomeRegion_peakFiles",resource_path+"fastQForSelRegion.bam")
    
    
test_macs2()