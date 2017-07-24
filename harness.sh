#!/bin/bash

# pytest -s tests/test_bwa_indexer.py
# pytest -s tests/test_bwa_aligner.py 

# echo "before bams, after tests"
# ls -l /home/travis/build/Multiscale-Genomics/mg-process-fastq/tests/data

# bgzip --help
# htsfile --help
# tabix --help

echo $PATH
set -x
bamsormadup --tmpfile="/home/travis/build/Multiscale-Genomics/mg-process-fastq/tests/data" </home/travis/build/Multiscale-Genomics/mg-process-fastq/tests/data/macs2.Human.DRR000150.22.bam >testFile.f.bam 
set +x
echo "after bams"
ls /home/travis/build/Multiscale-Genomics/mg-process-fastq/tests/data

pytest -s tests/test_biobambam.py
echo "after bams test"
ls /home/travis/build/Multiscale-Genomics/mg-process-fastq/tests/data
