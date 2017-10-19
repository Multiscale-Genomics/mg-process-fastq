#!/bin/bash

python tests/test_toolchains.py --pipeline genome --verbose 1
tc1=$?
ls tests/data/
./tidy_data.sh
ls tests/data/

python tests/test_toolchains.py --pipeline chipseq --verbose 1
tc2=$?
ls tests/data/
./tidy_data.sh
ls tests/data/

python tests/test_toolchains.py --pipeline rnaseq
tc3=$?
ls tests/data/
./tidy_data.sh
ls tests/data/

python tests/test_toolchains.py --pipeline wgbs
tc4=$?
ls tests/data/
./tidy_data.sh
ls tests/data/

rc=$(($tc1+$tc2+$tc3+$tc4))

if [[ $rc != 0 ]]; then exit $rc; fi

# pytest  tests/test_bwa_indexer.py
# pytest  tests/test_bwa_aligner.py
# pytest  tests/test_biobambam.py
# pytest  tests/test_bs_seeker_indexer.py
# pytest  tests/test_bs_seeker_aligner.py
# pytest  tests/test_bs_seeker_filter.py
# samtools
# ls /home/travis/build/Multiscale-Genomics/mg-process-fastq/tests/data
# pytest  tests/test_bs_seeker_methylation_caller.py
# pytest -s tests/test_bowtie_indexer.py
# pytest -s tests/test_kallisto_indexer.py
# pytest -s tests/test_kallisto_quant.py
# pytest -s tests/test_macs2.py
# pytest -s tests/test_paired_splitter.py
# pytest -s tests/test_single_splitter.py

# python_version=$(python --version 2>&1)

# echo $python_version
# if [[ $python_version == *"3."* ]]; then
#     pytest -s tests/test_inps.py
# fi
