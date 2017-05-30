#!usr/bin/python

import pytest
import random

from tool import macs2

def test_macs2():
    m = macs2.macs2()
    m.macs2_peak_calling("I_Dont_know","/Users/reham/Documents/fastQForSelRegion.bam")
    
    
test_macs2()