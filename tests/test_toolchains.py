"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.

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

    params = ['-m genome']

    if verbose is True:
        params.append('-s')

    params.append('tests/test_bowtie_indexer.py')
    params.append('tests/test_bwa_indexer.py')
    params.append('tests/test_gem_indexer.py')

    return pytest.main(params)

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

    params = ['-m chipseq']

    if verbose is True:
        params.append('-s')

    params.append('tests/test_bwa_indexer.py')
    params.append('tests/test_bwa_aligner.py')
    params.append('tests/test_biobambam.py')
    params.append('tests/test_macs2.py')

    return pytest.main(params)

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

    params = ['-m hic']

    if verbose is True:
        params.append('-s')

    params.append('tests/test_gem_indexer.py')
    params.append('tests/test_tb_full_mapping.py')
    params.append('tests/test_tb_parse_mapping.py')
    params.append('tests/test_tb_filter.py')
    params.append('tests/test_tb_generate_tads.py')
    params.append('tests/test_tb_save_hdf5_matrix.py')

    return pytest.main(params)

def mnaseseq_toolchain(verbose=False):
    """
    Runs the tests for all of the tools from the MNase-seq pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m mnaseseq tests/test_bwa_indexer.py
       pytest -m mnaseseq tests/test_bwa_aligner.py
       pytest -m mnaseseq tests/test_inps.py
    """

    params = ['-m mnaseseq']

    if verbose is True:
        params.append('-s')

    params.append('tests/test_bwa_indexer.py')
    params.append('tests/test_bwa_aligner.py')
    params.append('tests/test_inps.py')

    return pytest.main(params)

def rnaseq_toolchain(verbose=False):
    """
    Runs the tests for all of the tools from the RNA-seq pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m rnaseq tests/test_kallisto_indexer.py
       pytest -m rnaseq tests/test_kallisto_quant.py
    """

    params = ['-m rnaseq']

    if verbose is True:
        params.append('-s')

    params.append('tests/test_kallisto_indexer.py')
    params.append('tests/test_kallisto_quant.py')

    return pytest.main(params)

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

    params = ['-m wgbs']

    if verbose is True:
        params.append('-s')

    params.append('tests/test_bs_seeker_filter.py')
    params.append('tests/test_bs_seeker_indexer.py')
    params.append('tests/test_bs_seeker_aligner.py')
    params.append('tests/test_bs_seeker_methylation_caller.py')

    return pytest.main(params)

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
    PARSER.add_argument("--verbose", action="store_const", const=True, default=False)


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

    # if 'hic' in PIPELINES or 'all' in PIPELINES:
    #     print('HIC')
    #     if hic_toolchain(VERBOSE) > 0:
    #         sys.exit(1)

    # if 'mnaseseq' in PIPELINES or 'all' in PIPELINES:
    #     print('MNASESEQ')
    #     if mnaseseq_toolchain(VERBOSE) > 0:
    #         sys.exit(1)

    if 'rnaseq' in PIPELINES or 'all' in PIPELINES:
        print('RNASEQ')
        if rnaseq_toolchain(VERBOSE) > 0:
            sys.exit(1)

    if 'wgbs' in PIPELINES or 'all' in PIPELINES:
        print('WGBS')
        if wgbs_toolchain(VERBOSE) > 0:
            sys.exit(1)
