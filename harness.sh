#!/bin/bash

pytest tests/test_bwa_indexer.py
pytest tests/test_bwa_aligner.py
ls tests/data
pytest tests/test_biobambam.py
