"""
.. Copyright 2017 EMBL-European Bioinformatics Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

from __future__ import print_function

import argparse
import pytest # pylint: disable=unused-import

def chipseq_toolchain():
    """
    Runs the tests for all of the tools from the ChIP-seq pipeline
    """
    pytest.main(
        [
            '-m chipseq',
            'tests/test_bwa_indexer.py',
            'tests/test_bwa_aligner.py',
            'tests/test_biobambam.py',
            'tests/test_macs2.py',
        ]
    )

def hic_toolchain():
    """
    Runs the tests for all of the tools from the Hi-C pipeline
    """
    pytest.main(
        [
            '-m hic',
            'tests/test_gem_indexer.py',
            'tests/test_tb_full_mapping.py',
            'tests/test_tb_parse_mapping.py',
            'tests/test_tb_filter.py',
            'tests/test_tb_generate_tads.py',
            'tests/test_tb_save_hdf5_matrix.py',
        ]
    )

def mnaseseq_toolchain():
    """
    Runs the tests for all of the tools from the MNase-seq pipeline
    """
    pytest.main(
        [
            '-m mnaseseq',
            'tests/test_bwa_indexer.py',
            'tests/test_bwa_aligner.py',
            'tests/test_inps.py',
        ]
    )

def rnaseq_toolchain():
    """
    Runs the tests for all of the tools from the RNA-seq pipeline
    """
    pytest.main(
        [
            '-m rnaseq',
            'tests/test_kallisto_indexer.py',
            'tests/test_kallisto_quant.py',
        ]
    )

def wgbs_toolchain():
    """
    Runs the tests for all of the tools from the WGBS pipeline
    """
    pytest.main(
        [
            '-m wgbs',
            'tests/test_bs_seeker_filter.py',
            'tests/test_bs_seeker_indexer.py',
            'tests/test_bs_seeker_aligner.py',
            'tests/test_bs_seeker_methylation_caller.py',
        ]
    )

if __name__ == '__main__':
    import sys
    sys._run_from_cmdl = True

    PARSER = argparse.ArgumentParser(description="Test runner for tool chains")
    PARSER.add_argument("--pipeline", default=None)

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    PIPELINES = ARGS.pipeline
    PIPELINES = PIPELINES.split(",")
    print(PIPELINES)

    if 'chipseq' in PIPELINES or 'all' in PIPELINES:
        print('CHIPSEQ')
        chipseq_toolchain()

    if 'hic' in PIPELINES or 'all' in PIPELINES:
        print('HIC')
        hic_toolchain()

    if 'mnaseseq' in PIPELINES or 'all' in PIPELINES:
        print('MNASESEQ')
        mnaseseq_toolchain()

    if 'rnaseq' in PIPELINES or 'all' in PIPELINES:
        print('RNASEQ')
        rnaseq_toolchain()

    if 'wgbs' in PIPELINES or 'all' in PIPELINES:
        print('WGBS')
        wgbs_toolchain()
