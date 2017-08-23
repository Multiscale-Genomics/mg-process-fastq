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

def genome_pipeline(verbose=False):
    """
    Runs the tests for the Genome indexing pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m genome tests/test_pipeline_genome.py
    """

    v_str = ''
    if verbose == 1:
        v_str = '-s'

    return pytest.main(
        [
            '-m genome',
            v_str,
            'tests/test_pipeline_genome.py',
        ]
    )

def chipseq_pipeline(verbose=False):
    """
    Runs the tests for the ChIP-seq pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m chipseq tests/test_pipeline_chipseq.py
    """

    v_str = ''
    if verbose == 1:
        v_str = '-s'

    return pytest.main(
        [
            '-m chipseq',
            v_str,
            'tests/test_chipseq.py',
        ]
    )

def hic_pipeline(verbose=False):
    """
    Runs the tests for the Hi-C pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m hic tests/test_pipeline_tb.py
    """

    v_str = ''
    if verbose == 1:
        v_str = '-s'

    return pytest.main(
        [
            '-m hic',
            v_str,
            'tests/test_pipeline_tb.py',
        ]
    )

def mnaseseq_pipeline(verbose=False):
    """
    Runs the tests for the MNase-seq pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m mnaseseq tests/test_pipeline_mnaseseq.py
    """

    v_str = ''
    if verbose == 1:
        v_str = '-s'

    return pytest.main(
        [
            '-m mnaseseq',
            v_str,
            'tests/test_pipeline_mnaseseq.py',
        ]
    )

def rnaseq_pipeline(verbose=False):
    """
    Runs the tests for the RNA-seq pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m rnaseq tests/test_pipeline_rnaseq.py
    """

    v_str = ''
    if verbose == 1:
        v_str = '-s'

    return pytest.main(
        [
            '-m rnaseq',
            v_str,
            'tests/test_pipeline_rnaseq.py',
        ]
    )

def wgbs_pipeline(verbose=0):
    """
    Runs the tests for the WGBS pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m wgbs tests/test_pipeline_wgbs.py
    """

    v_str = ''
    if verbose == 1:
        v_str = '-s'

    return pytest.main(
        [
            '-m wgbs',
            v_str,
            'tests/test_pipeline_wgbs.py',
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
        if genome_pipeline(VERBOSE) > 0:
            sys.exit(1)

    if 'chipseq' in PIPELINES or 'all' in PIPELINES:
        print('CHIPSEQ')
        if chipseq_pipeline(VERBOSE) > 0:
            sys.exit(1)

    if 'hic' in PIPELINES or 'all' in PIPELINES:
        print('HIC')
        if hic_pipeline(VERBOSE) > 0:
            sys.exit(1)

    if 'mnaseseq' in PIPELINES or 'all' in PIPELINES:
        print('MNASESEQ')
        if mnaseseq_pipeline(VERBOSE) > 0:
            sys.exit(1)

    if 'rnaseq' in PIPELINES or 'all' in PIPELINES:
        print('RNASEQ')
        if rnaseq_pipeline(VERBOSE) > 0:
            sys.exit(1)

    if 'wgbs' in PIPELINES or 'all' in PIPELINES:
        print('WGBS')
        if wgbs_pipeline(VERBOSE) > 0:
            sys.exit(1)
