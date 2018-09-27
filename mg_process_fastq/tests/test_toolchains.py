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
import subprocess
import pytest


def all_toolchain(verbose=False):
    """
    Runs the tests for all of the tools

    This set is only required for determining code coverage.
    """

    params = []

    if verbose is True:
        params.append('-s')

    params.append('mg_process_fastq/tests/test_bowtie_indexer.py')
    params.append('mg_process_fastq/tests/test_bwa_indexer.py')
    params.append('mg_process_fastq/tests/test_gem_indexer.py')

    params.append('mg_process_fastq/tests/test_fastqc_validation.py')

    params.append('mg_process_fastq/tests/test_bowtie2_aligner.py')
    params.append('mg_process_fastq/tests/test_bwa_aligner.py')

    params.append('mg_process_fastq/tests/test_biobambam.py')

    params.append('mg_process_fastq/tests/test_macs2.py')

    # params.append('mg_process_fastq/tests/test_bsgenome.py')
    # params.append('mg_process_fastq/tests/test_idear.py')

    params.append('mg_process_fastq/tests/test_inps.py')

    params.append('mg_process_fastq/tests/test_kallisto_indexer.py')
    params.append('mg_process_fastq/tests/test_kallisto_quant.py')

    params.append('mg_process_fastq/tests/test_bs_seeker_filter.py')
    params.append('mg_process_fastq/tests/test_bs_seeker_indexer.py')
    params.append('mg_process_fastq/tests/test_bs_seeker_aligner.py')
    params.append('mg_process_fastq/tests/test_bs_seeker_methylation_caller.py')

    return pytest.main(params)


def genome_toolchain(verbose=False):
    """
    Runs the tests for all of the tools from the Genome indexing pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m genome mg_process_fastq/tests/test_bowtie_indexer.py
       pytest -m genome mg_process_fastq/tests/test_bwa_indexer.py
       pytest -m genome mg_process_fastq/tests/test_gem_indexer.py
    """

    params = ['-m genome']

    if verbose is True:
        params.append('-s')

    params.append('mg_process_fastq/tests/test_bowtie_indexer.py')
    params.append('mg_process_fastq/tests/test_bwa_indexer.py')
    params.append('mg_process_fastq/tests/test_gem_indexer.py')

    return pytest.main(params)


def biobambam_toolchain(verbose=False):
    """
    Runs the tests for all of the tools from the BWA pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m chipseq mg_process_fastq/tests/test_fastqc_validation.py
       pytest -m chipseq mg_process_fastq/tests/test_bwa_indexer.py
       pytest -m chipseq mg_process_fastq/tests/test_bwa_aligner.py
       pytest -m chipseq mg_process_fastq/tests/test_biobambam.py
    """

    params = ['-m chipseq']

    if verbose is True:
        params.append('-s')

    params.append('mg_process_fastq/tests/test_fastqc_validation.py')
    params.append('mg_process_fastq/tests/test_bwa_indexer.py')
    params.append('mg_process_fastq/tests/test_bwa_aligner.py')
    params.append('mg_process_fastq/tests/test_biobambam.py')

    return pytest.main(params)


def bowtie2_toolchain(verbose=False):
    """
    Runs the tests for all of the tools from the BWA pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m bowtie2 mg_process_fastq/tests/test_fastqc_validation.py
       pytest -m bowtie2 mg_process_fastq/tests/test_bowtie_indexer.py
       pytest -m bowtie2 mg_process_fastq/tests/test_bowtie2_aligner.py
    """

    params = ['-m bowtie2']

    if verbose is True:
        params.append('-s')

    params.append('mg_process_fastq/tests/test_fastqc_validation.py')
    params.append('mg_process_fastq/tests/test_bowtie_indexer.py')
    params.append('mg_process_fastq/tests/test_bowtie2_aligner.py')

    return pytest.main(params)


def bwa_toolchain(verbose=False):
    """
    Runs the tests for all of the tools from the BWA pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m bwa mg_process_fastq/tests/test_fastqc_validation.py
       pytest -m bwa mg_process_fastq/tests/test_bwa_indexer.py
       pytest -m bwa mg_process_fastq/tests/test_bwa_aligner.py
    """

    params = ['-m bwa']

    if verbose is True:
        params.append('-s')

    params.append('mg_process_fastq/tests/test_fastqc_validation.py')
    params.append('mg_process_fastq/tests/test_bwa_indexer.py')
    params.append('mg_process_fastq/tests/test_bwa_aligner.py')

    return pytest.main(params)


def chipseq_toolchain(verbose=False):
    """
    Runs the tests for all of the tools from the ChIP-seq pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m chipseq mg_process_fastq/tests/test_fastqc_validation.py
       pytest -m chipseq mg_process_fastq/tests/test_bwa_indexer.py
       pytest -m chipseq mg_process_fastq/tests/test_bwa_aligner.py
       pytest -m chipseq mg_process_fastq/tests/test_biobambam.py
       pytest -m chipseq mg_process_fastq/tests/test_macs2.py
    """

    params = ['-m chipseq']

    if verbose is True:
        params.append('-s')

    params.append('mg_process_fastq/tests/test_fastqc_validation.py')
    params.append('mg_process_fastq/tests/test_bwa_indexer.py')
    params.append('mg_process_fastq/tests/test_bwa_aligner.py')
    params.append('mg_process_fastq/tests/test_biobambam.py')
    params.append('mg_process_fastq/tests/test_macs2.py')

    return pytest.main(params)


def idamidseq_toolchain(verbose=False):
    """
    Runs the tests for all of the tools from the ChIP-seq pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m idamidseq mg_process_fastq/tests/test_fastqc_validation.py
       pytest -m idamidseq mg_process_fastq/tests/test_bwa_indexer.py
       pytest -m idamidseq mg_process_fastq/tests/test_bwa_aligner.py
       pytest -m idamidseq mg_process_fastq/tests/test_biobambam.py
       pytest -m idamidseq mg_process_fastq/tests/test_bsgenome.py
       pytest -m idamidseq mg_process_fastq/tests/test_idear.py
    """

    params = ['-m idamidseq']

    if verbose is True:
        params.append('-s')

    params.append('mg_process_fastq/tests/test_fastqc_validation.py')
    params.append('mg_process_fastq/tests/test_bwa_indexer.py')
    params.append('mg_process_fastq/tests/test_bwa_aligner.py')
    params.append('mg_process_fastq/tests/test_biobambam.py')
    params.append('mg_process_fastq/tests/test_bsgenome.py')
    params.append('mg_process_fastq/tests/test_idear.py')

    return pytest.main(params)


def hic_toolchain(verbose=False):
    """
    Runs the tests for all of the tools from the Hi-C pipeline

    Runs the following tests:

    .. code-block:: none

<<<<<<< HEAD:mg_process_fastq/tests/test_toolchains.py
       pytest -m hic mg_process_fastq/tests/test_fastqc_validation.py
       pytest -m hic mg_process_fastq/tests/test_gem_indexer.py
       pytest -m hic mg_process_fastq/tests/test_tb_full_mapping.py
       pytest -m hic mg_process_fastq/tests/test_tb_parse_mapping.py
       pytest -m hic mg_process_fastq/tests/test_tb_filter.py
       pytest -m hic mg_process_fastq/tests/test_tb_generate_tads.py
       pytest -m hic mg_process_fastq/tests/test_tb_save_hdf5_matrix.py
=======
       pytest -m hic tests/test_fastqc_validation.py
       pytest -m hic tests/test_gem_indexer.py
       pytest -m hic tests/test_tb_full_mapping.py
       pytest -m hic tests/test_tb_parse_mapping.py
       pytest -m hic tests/test_tb_filter.py
       pytest -m hic tests/test_tb_normalize.py
       pytest -m hic tests/test_tb_segment.py
       pytest -m hic tests/test_tb_generate_tads.py
       pytest -m hic tests/test_tb_bin.py
       pytest -m hic tests/test_tb_save_hdf5_matrix.py
>>>>>>> master:tests/test_toolchains.py
    """

    params = ['-m hic']

    if verbose is True:
        params.append('-s')

    params.append('tests/test_fastqc_validation.py')
    params.append('tests/test_gem_indexer.py')
    params.append('tests/test_tb_full_mapping.py')
    params.append('tests/test_tb_parse_mapping.py')
    params.append('tests/test_tb_filter.py')
    params.append('tests/test_tb_normalize.py')
    params.append('tests/test_tb_segment.py')
    params.append('tests/test_tb_generate_tads.py')
    params.append('tests/test_tb_bin.py')
    params.append('tests/test_tb_model.py')
    params.append('tests/test_tb_save_hdf5_matrix.py')

    return pytest.main(params)


def mnaseseq_toolchain(verbose=False):
    """
    Runs the tests for all of the tools from the MNase-seq pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m mnaseseq mg_process_fastq/tests/test_fastqc_validation.py
       pytest -m mnaseseq mg_process_fastq/tests/test_bwa_indexer.py
       pytest -m mnaseseq mg_process_fastq/tests/test_bwa_aligner.py
       pytest -m mnaseseq mg_process_fastq/tests/test_inps.py
    """

    params = ['-m mnaseseq']

    if verbose is True:
        params.append('-s')

    params.append('mg_process_fastq/tests/test_fastqc_validation.py')
    params.append('mg_process_fastq/tests/test_bwa_indexer.py')
    params.append('mg_process_fastq/tests/test_bwa_aligner.py')
    params.append('mg_process_fastq/tests/test_inps.py')

    return pytest.main(params)


def rnaseq_toolchain(verbose=False):
    """
    Runs the tests for all of the tools from the RNA-seq pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m rnaseq mg_process_fastq/tests/test_fastqc_validation.py
       pytest -m rnaseq mg_process_fastq/tests/test_kallisto_indexer.py
       pytest -m rnaseq mg_process_fastq/tests/test_kallisto_quant.py
    """

    params = ['-m rnaseq']

    if verbose is True:
        params.append('-s')

    params.append('mg_process_fastq/tests/test_fastqc_validation.py')
    params.append('mg_process_fastq/tests/test_kallisto_indexer.py')
    params.append('mg_process_fastq/tests/test_kallisto_quant.py')

    return pytest.main(params)


def wgbs_toolchain(verbose=0):
    """
    Runs the tests for all of the tools from the WGBS pipeline

    Runs the following tests:

    .. code-block:: none

       pytest -m wgbs mg_process_fastq/tests/test_fastqc_validation.py
       pytest -m wgbs mg_process_fastq/tests/test_bs_seeker_filter.py
       pytest -m wgbs mg_process_fastq/tests/test_bs_seeker_indexer.py
       pytest -m wgbs mg_process_fastq/tests/test_bs_seeker_aligner.py
       pytest -m wgbs mg_process_fastq/tests/test_bs_seeker_methylation_caller.py
    """

    params = ['-m wgbs']

    if verbose is True:
        params.append('-s')

    params.append('mg_process_fastq/tests/test_fastqc_validation.py')
    params.append('mg_process_fastq/tests/test_bs_seeker_filter.py')
    params.append('mg_process_fastq/tests/test_bs_seeker_indexer.py')
    params.append('mg_process_fastq/tests/test_bs_seeker_aligner.py')
    params.append('mg_process_fastq/tests/test_bs_seeker_methylation_caller.py')

    return pytest.main(params)


def tidy_data():
    """
    Runs the tidy_data.sh script
    """
    print("TIDY DATA")
    try:
        command_line = './tidy_data.sh'
        process = subprocess.Popen(command_line, shell=True)
        process.wait()
    except (IOError, OSError) as msg:
        print("I/O error({0}): {1}\n{2}".format(
            msg.errno, msg.strerror, command_line))


if __name__ == '__main__':
    import sys
    sys._run_from_cmdl = True  # pylint: disable=protected-access

    PARSER = argparse.ArgumentParser(description="Test runner for tool chains")
    PARSER.add_argument(
        "--pipeline",
        required=True,
        type=str,
        choices=[
            'genome', 'biobambam', 'bowtie2', 'bwa', 'chipseq', 'hic', 'idamidseq', 'mnaseseq',
            'rnaseq', 'wgbs', 'all'
        ],
        help=""
    )
    PARSER.add_argument("--verbose", action="store_const", const=True, default=False)
    PARSER.add_argument("--tidy", action="store_const", const=True, default=False)

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    PIPELINES = ARGS.pipeline

    if PIPELINES is None:
        PARSER.print_help()
        sys.exit(1)

    PIPELINES = PIPELINES.split(",")
    print(PIPELINES)

    VERBOSE = ARGS.verbose
    TIDY = ARGS.tidy

    if 'all' in PIPELINES:
        print('ALL')

        if all_toolchain(VERBOSE) > 0:
            sys.exit(1)

        if TIDY:
            tidy_data()

    if 'genome' in PIPELINES:
        print('GENOME')

        if genome_toolchain(VERBOSE) > 0:
            sys.exit(1)

        if TIDY:
            tidy_data()

    if 'biobambam' in PIPELINES:
        print('BIOBAMBAM')

        if biobambam_toolchain(VERBOSE) > 0:
            sys.exit(1)

        if TIDY:
            tidy_data()

    if 'bowtie2' in PIPELINES:
        print('BOWTIE2')

        if bowtie2_toolchain(VERBOSE) > 0:
            sys.exit(1)

        if TIDY:
            tidy_data()

    if 'bwa' in PIPELINES:
        print('BWA')

        if bwa_toolchain(VERBOSE) > 0:
            sys.exit(1)

        if TIDY:
            tidy_data()

    if 'chipseq' in PIPELINES:
        print('CHIPSEQ')
        if chipseq_toolchain(VERBOSE) > 0:
            sys.exit(1)

        if TIDY:
            tidy_data()

    if 'hic' in PIPELINES:
        print('HIC')
        if hic_toolchain(VERBOSE) > 0:
            sys.exit(1)

        if TIDY:
            tidy_data()

    if 'idamidseq' in PIPELINES:
        print('IDAMIDSEQ')
        if idamidseq_toolchain(VERBOSE) > 0:
            sys.exit(1)

        if TIDY:
            tidy_data()

    if 'mnaseseq' in PIPELINES:
        print('MNASESEQ')
        if mnaseseq_toolchain(VERBOSE) > 0:
            sys.exit(1)

        if TIDY:
            tidy_data()

    if 'rnaseq' in PIPELINES:
        print('RNASEQ')
        if rnaseq_toolchain(VERBOSE) > 0:
            sys.exit(1)

        if TIDY:
            tidy_data()

    if 'wgbs' in PIPELINES:
        print('WGBS')
        if wgbs_toolchain(VERBOSE) > 0:
            sys.exit(1)

        if TIDY:
            tidy_data()
