#!/bin/bash

pytest -s tests/test_bwa_indexer.py
pytest -s tests/test_bwa_aligner.py 
pytest -s tests/test_biobambam.py
pytest -s tests/test_bs_seeker_indexer.py
pytest -s tests/test_bs_seeker_aligner.py
pytest -s tests/test_bs_seeker_filter.py
pytest -s tests/test_bs_seeker_methylation_caller.py
pytest -s tests/test_bowtie_indexer.py
pytest -s tests/test_kallisto_indexer.py
pytest -s tests/test_kallisto_quant.py
pytest -s tests/test_macs2.py
pytest -s tests/test_paired_splitter.py
pytest -s tests/test_single_splitter.py
