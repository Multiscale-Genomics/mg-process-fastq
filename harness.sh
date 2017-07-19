#!/bin/bash

pytest -s tests/test_bwa_indexer.py
pytest -s tests/test_bwa_aligner.py 
pytest -s tests/test_biobambam.py
