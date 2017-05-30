#!usr/bin/python

import pytest
import random
import os.path

from tool import bwa_aligner

def test_bwa_aligner():
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genome_fa = resource_path+ "selectGenomeRegion.fasta"
    fastqFile = resource_path+ "fastQForSelRegion.fastq"
    out_bam = fastqFile.replace('.fastq', '.bam')
    
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
    
    bwaT = bwa_aligner.bwaAlignerTool()
    bwaT.bwa_aligner(genome_fa, fastqFile, out_bam, bwa_amb, bwa_ann, bwa_bwt,bwa_pac, bwa_sa) 
    
    print __file__

test_bwa_aligner()

