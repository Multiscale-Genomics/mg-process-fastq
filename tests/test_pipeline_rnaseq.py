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

import os.path
import pytest

from basic_modules.metadata import Metadata
from process_rnaseq import process_rnaseq


@pytest.mark.rnaseq
@pytest.mark.pipeline
def test_rnaseq_pipeline():
    """
    Test case to ensure that the RNA-seq pipeline code works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    files = {
        'cdna': resource_path + 'kallisto.Human.GRCh38.fasta',
        'fastq1': resource_path + 'kallisto.Human.ERR030872_1.fastq',
        'fastq2': resource_path + 'kallisto.Human.ERR030872_2.fastq'
    }

    metadata = {
        "cdna": Metadata(
            "Assembly", "fasta", files['cdna'], None,
            {'assembly': 'GCA_000001405.22'}),
        "fastq1": Metadata(
            "data_rna_seq", "fastq", files['fastq1'], None,
            {'assembly': 'GCA_000001405.22'}
        ),
        "fastq2": Metadata(
            "data_rna_seq", "fastq", files['fastq2'], None,
            {'assembly': 'GCA_000001405.22'}
        ),
    }

    files_out = {
        "index": 'tests/data/kallisto.idx',
        "kallisto_tar": 'tests/data/kallisto.tar.gz'
    }

    rs_handle = process_rnaseq()
    rs_files, rs_meta = rs_handle.run(files, metadata, files_out)  # pylint: disable=unused-variable

    # Checks that the returned files matches the expected set of results
    assert len(rs_files) == 4

    # Add tests for all files created
    for f_out in rs_files:
        print("RNA SEQ RESULTS FILE:", f_out)
        assert rs_files[f_out] == files_out[f_out]
        assert os.path.isfile(rs_files[f_out]) is True
        assert os.path.getsize(rs_files[f_out]) > 0


@pytest.mark.sleuth
@pytest.mark.pipeline
def test_rnaseq_pipeline_sleuth_00():
    """
    Test case to ensure that the Sleuth pipeline code works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    files = {
        'cdna': resource_path + 'sleuth.Human.GRCh38.cdna.fasta',
        'fastq1': resource_path + 'sleuth.Human.ERR030872_1.fastq',
        'fastq2': resource_path + 'sleuth.Human.ERR030872_2.fastq'
    }

    metadata = {
        "cdna": Metadata(
            "Assembly", "fasta", files['cdna'], None,
            {'assembly': 'GRCh38'}),
        "fastq1": Metadata(
            "data_rna_seq", "fastq", files['fastq1'], None,
            {'assembly': 'GRCh38'}
        ),
        "fastq2": Metadata(
            "data_rna_seq", "fastq", files['fastq2'], None,
            {'assembly': 'GRCh38'}
        ),
    }

    files_out = {
        "index": 'tests/data/kallisto.idx',
        "kallisto_tar": 'tests/data/ERR030872.tar.gz'
    }

    rs_handle = process_rnaseq()
    rs_files, rs_meta = rs_handle.run(files, metadata, files_out)  # pylint: disable=unused-variable

    # Checks that the returned files matches the expected set of results
    assert len(rs_files) == 4

    # Add tests for all files created
    for f_out in rs_files:
        print("RNA SEQ RESULTS FILE:", f_out)
        assert rs_files[f_out] == files_out[f_out]
        assert os.path.isfile(rs_files[f_out]) is True
        assert os.path.getsize(rs_files[f_out]) > 0


@pytest.mark.sleuth
@pytest.mark.pipeline
def test_rnaseq_pipeline_sleuth_01():
    """
    Test case to ensure that the Sleuth pipeline code works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    files = {
        'cdna': resource_path + 'sleuth.Human.GRCh38.cdna.fasta',
        'fastq1': resource_path + 'sleuth.Human.ERR030856.fastq',
    }

    metadata = {
        "cdna": Metadata(
            "Assembly", "fasta", files['cdna'], None,
            {'assembly': 'GRCh38'}),
        "fastq1": Metadata(
            "data_rna_seq", "fastq", files['fastq1'], None,
            {'assembly': 'GRCh38'}
        ),
    }

    files_out = {
        "index": 'tests/data/kallisto.idx',
        "kallisto_tar": 'tests/data/ERR030856.tar.gz'
    }

    rs_handle = process_rnaseq()
    rs_files, rs_meta = rs_handle.run(files, metadata, files_out)  # pylint: disable=unused-variable

    # Checks that the returned files matches the expected set of results
    assert len(rs_files) == 4

    # Add tests for all files created
    for f_out in rs_files:
        print("RNA SEQ RESULTS FILE:", f_out)
        assert rs_files[f_out] == files_out[f_out]
        assert os.path.isfile(rs_files[f_out]) is True
        assert os.path.getsize(rs_files[f_out]) > 0


@pytest.mark.sleuth
@pytest.mark.pipeline
def test_rnaseq_pipeline_sleuth_01():
    """
    Test case to ensure that the Sleuth pipeline code works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    files = {
        'cdna': resource_path + 'sleuth.Human.GRCh38.cdna.fasta',
        'fastq1': resource_path + 'sleuth.Human.ERR030857.fastq',
    }

    metadata = {
        "cdna": Metadata(
            "Assembly", "fasta", files['cdna'], None,
            {'assembly': 'GRCh38'}),
        "fastq1": Metadata(
            "data_rna_seq", "fastq", files['fastq1'], None,
            {'assembly': 'GRCh38'}
        ),
    }

    files_out = {
        "index": 'tests/data/kallisto.idx',
        "kallisto_tar": 'tests/data/ERR030857.tar.gz'
    }

    rs_handle = process_rnaseq()
    rs_files, rs_meta = rs_handle.run(files, metadata, files_out)  # pylint: disable=unused-variable

    # Checks that the returned files matches the expected set of results
    assert len(rs_files) == 4

    # Add tests for all files created
    for f_out in rs_files:
        print("RNA SEQ RESULTS FILE:", f_out)
        assert rs_files[f_out] == files_out[f_out]
        assert os.path.isfile(rs_files[f_out]) is True
        assert os.path.getsize(rs_files[f_out]) > 0


@pytest.mark.sleuth
@pytest.mark.pipeline
def test_rnaseq_pipeline_sleuth_01():
    """
    Test case to ensure that the Sleuth pipeline code works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    files = {
        'cdna': resource_path + 'sleuth.Human.GRCh38.cdna.fasta',
        'fastq1': resource_path + 'sleuth.Human.ERR030858.fastq',
    }

    metadata = {
        "cdna": Metadata(
            "Assembly", "fasta", files['cdna'], None,
            {'assembly': 'GRCh38'}),
        "fastq1": Metadata(
            "data_rna_seq", "fastq", files['fastq1'], None,
            {'assembly': 'GRCh38'}
        ),
    }

    files_out = {
        "index": 'tests/data/kallisto.idx',
        "kallisto_tar": 'tests/data/ERR030858.tar.gz'
    }

    rs_handle = process_rnaseq()
    rs_files, rs_meta = rs_handle.run(files, metadata, files_out)  # pylint: disable=unused-variable

    # Checks that the returned files matches the expected set of results
    assert len(rs_files) == 4

    # Add tests for all files created
    for f_out in rs_files:
        print("RNA SEQ RESULTS FILE:", f_out)
        assert rs_files[f_out] == files_out[f_out]
        assert os.path.isfile(rs_files[f_out]) is True
        assert os.path.getsize(rs_files[f_out]) > 0


@pytest.mark.sleuth
@pytest.mark.pipeline
def test_rnaseq_pipeline_sleuth_01():
    """
    Test case to ensure that the Sleuth pipeline code works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    files = {
        'cdna': resource_path + 'sleuth.Human.GRCh38.cdna.fasta',
        'fastq1': resource_path + 'sleuth.Human.ERR030903.fastq',
    }

    metadata = {
        "cdna": Metadata(
            "Assembly", "fasta", files['cdna'], None,
            {'assembly': 'GRCh38'}),
        "fastq1": Metadata(
            "data_rna_seq", "fastq", files['fastq1'], None,
            {'assembly': 'GRCh38'}
        ),
    }

    files_out = {
        "index": 'tests/data/kallisto.idx',
        "kallisto_tar": 'tests/data/ERR030903.tar.gz'
    }

    rs_handle = process_rnaseq()
    rs_files, rs_meta = rs_handle.run(files, metadata, files_out)  # pylint: disable=unused-variable

    # Checks that the returned files matches the expected set of results
    assert len(rs_files) == 4

    # Add tests for all files created
    for f_out in rs_files:
        print("RNA SEQ RESULTS FILE:", f_out)
        assert rs_files[f_out] == files_out[f_out]
        assert os.path.isfile(rs_files[f_out]) is True
        assert os.path.getsize(rs_files[f_out]) > 0
