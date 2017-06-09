#!usr/bin/python

import pytest
import random
import os.path

from tool import bwa_indexer

def test_bwa_indexer ():
    
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genome_fa = resource_path + "macs2.Human.GCA_000001405.22.fasta"
    fastqFile = resource_path + "macs2.Human.DRR000150.22.fastq"
    
    files = [
        genome_fa,
        genome_fa + ".amb",
        genome_fa + ".ann",
        genome_fa + ".bwt",
        genome_fa + ".pac",
        genome_fa + ".sa"
    ]
    
    bwa_amb = files[1]
    bwa_ann = files[2]
    bwa_bwt = files[3]
    bwa_pac = files[4]
    bwa_sa  = files[5]
    
    print (bwa_amb)
    print (bwa_ann)
    print (bwa_bwt)
    print (bwa_pac)
    print (bwa_sa)
    
    bwaT = bwa_indexer.bwaIndexerTool()
    bwaT.run([genome_fa], {'assembly' : 'test'}) 
    
def test_bwa_indexer_02 ():
    
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genome_fa = resource_path + "inps.Mouse.GRCm38.fasta"
    fastqFile = resource_path + "fastQForMouseRegion.fastq"
    
    files = [
        genome_fa,
        genome_fa + ".amb",
        genome_fa + ".ann",
        genome_fa + ".bwt",
        genome_fa + ".pac",
        genome_fa + ".sa"
    ]
    
    bwa_amb = files[1]
    bwa_ann = files[2]
    bwa_bwt = files[3]
    bwa_pac = files[4]
    bwa_sa  = files[5]
    
    print (bwa_amb)
    print (bwa_ann)
    print (bwa_bwt)
    print (bwa_pac)
    print (bwa_sa)
    
    bwaT = bwa_indexer.bwaIndexerTool()
    bwaT.run([genome_fa], {'assembly' : 'test'}) 
    
