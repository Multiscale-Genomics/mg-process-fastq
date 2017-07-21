#!/bin/bash

pytest -s tests/test_bwa_indexer.py
pytest -s tests/test_bwa_aligner.py 
exec ls /home/travis/build/Multiscale-Genomics/mg-process-fastq/tests/data
exec bamsormadup --tmpfile="/home/travis/build/Multiscale-Genomics/mg-process-fastq/tests/data" </home/travis/build/Multiscale-Genomics/mg-process-fastq/tests/data/macs2.Human.DRR000150.22.bam >testFile.f.bam
pytest -s tests/test_biobambam.py
