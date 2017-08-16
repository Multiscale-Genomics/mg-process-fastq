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

def genome_toolchain(verbose=False):
    """
    Runs the tests for all of the tools from the Genome indexing pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m genome tests/test_bowtie_indexer.py
       pytest -m genome tests/test_bwa_indexer.py
       pytest -m genome tests/test_gem_indexer.py
    """

    v_str = ''
    if verbose == 1:
        v_str = '-s'

    return pytest.main(
        [
            '-m genome',
            v_str,
            'tests/test_bowtie_indexer.py',
            'tests/test_bwa_indexer.py',
            'tests/test_gem_indexer.py',
        ]
    )

def chipseq_toolchain(verbose=False):
    """
    Runs the tests for all of the tools from the ChIP-seq pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m chipseq tests/test_bwa_indexer.py
       pytest -m chipseq tests/test_bwa_aligner.py
       pytest -m chipseq tests/test_biobambam.py
       pytest -m chipseq tests/test_macs2.py
    """

    v_str = ''
    if verbose == 1:
        v_str = '-s'

    return pytest.main(
        [
            '-m chipseq',
            v_str,
            'tests/test_bwa_indexer.py',
            'tests/test_bwa_aligner.py',
            'tests/test_biobambam.py',
            'tests/test_macs2.py',
        ]
    )

def hic_toolchain(verbose=False):
    """
    Runs the tests for all of the tools from the Hi-C pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m hic tests/test_gem_indexer.py
       pytest -m hic tests/test_tb_full_mapping.py
       pytest -m hic tests/test_tb_parse_mapping.py
       pytest -m hic tests/test_tb_filter.py
       pytest -m hic tests/test_tb_generate_tads.py
       pytest -m hic tests/test_tb_save_hdf5_matrix.py
    """

    v_str = ''
    if verbose == 1:
        v_str = '-s'

    return pytest.main(
        [
            '-m hic',
            v_str,
            'tests/test_gem_indexer.py',
            'tests/test_tb_full_mapping.py',
            'tests/test_tb_parse_mapping.py',
            'tests/test_tb_filter.py',
            'tests/test_tb_generate_tads.py',
            'tests/test_tb_save_hdf5_matrix.py',
        ]
    )

def mnaseseq_toolchain(verbose=False):
    """
    Runs the tests for all of the tools from the MNase-seq pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m mnaseseq tests/test_bwa_indexer.py
       pytest -m mnaseseq tests/test_bwa_aligner.py
       pytest -m mnaseseq tests/test_inps.py
    """

    v_str = ''
    if verbose == 1:
        v_str = '-s'

    return pytest.main(
        [
            '-m mnaseseq',
            v_str,
            'tests/test_bwa_indexer.py',
            'tests/test_bwa_aligner.py',
            'tests/test_inps.py',
        ]
    )

def rnaseq_toolchain(verbose=False):
    """
    Runs the tests for all of the tools from the RNA-seq pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m rnaseq tests/test_kallisto_indexer.py
       pytest -m rnaseq tests/test_kallisto_quant.py
    """

    v_str = ''
    if verbose == 1:
        v_str = '-s'

    return pytest.main(
        [
            '-m rnaseq',
            v_str,
            'tests/test_kallisto_indexer.py',
            'tests/test_kallisto_quant.py',
        ]
    )

def wgbs_toolchain(verbose=0):
    """
    Runs the tests for all of the tools from the WGBS pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m wgbs tests/test_bs_seeker_filter.py
       pytest -m wgbs tests/test_bs_seeker_indexer.py
       pytest -m wgbs tests/test_fastq_splitter.py
       pytest -m wgbs tests/test_bs_seeker_aligner.py
       pytest -m wgbs tests/test_bs_seeker_methylation_caller.py
    """

    v_str = ''
    if verbose == 1:
        v_str = '-s'

    return pytest.main(
        [
            '-m wgbs',
            v_str,
            'tests/test_bs_seeker_filter.py',
            'tests/test_bs_seeker_indexer.py',
            'tests/test_fastq_splitter.py',
            'tests/test_bs_seeker_aligner.py',
            'tests/test_bs_seeker_methylation_caller.py',
        ]
    )

if __name__ == '__main__':
    import sys
    sys._run_from_cmdl = True

    PARSER = argparse.ArgumentParser(description="Test runner for tool chains")
    PARSER.add_argument(
        "--pipeline",
        required=True,
        type=str,
        choices=['genome', 'chipseq', 'hic', 'mnaseseq', 'rnaseq', 'wgbs', 'all'],
        help=""
    )
    PARSER.add_argument("--verbose", type=int, default=0)

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()
    #print(ARGS)

    PIPELINES = ARGS.pipeline

    if PIPELINES is None:
        PARSER.print_help()
        sys.exit(1)

    PIPELINES = PIPELINES.split(",")
    print(PIPELINES)

    VERBOSE = ARGS.verbose

    if 'genome' in PIPELINES or 'all' in PIPELINES:
        print('GENOME')
        if genome_toolchain(VERBOSE) > 0:
            sys.exit(1)

    if 'chipseq' in PIPELINES or 'all' in PIPELINES:
        print('CHIPSEQ')
        if chipseq_toolchain(VERBOSE) > 0:
            sys.exit(1)

    if 'hic' in PIPELINES or 'all' in PIPELINES:
        print('HIC')
        if hic_toolchain(VERBOSE) > 0:
            sys.exit(1)

    if 'mnaseseq' in PIPELINES or 'all' in PIPELINES:
        print('MNASESEQ')
        if mnaseseq_toolchain(VERBOSE) > 0:
            sys.exit(1)

    if 'rnaseq' in PIPELINES or 'all' in PIPELINES:
        print('RNASEQ')
        if rnaseq_toolchain(VERBOSE) > 0:
            sys.exit(1)

    if 'wgbs' in PIPELINES or 'all' in PIPELINES:
        print('WGBS')
        if wgbs_toolchain(VERBOSE) > 0:
            sys.exit(1)
