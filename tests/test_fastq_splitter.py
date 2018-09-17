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

import os.path
import pytest
import gzip

from basic_modules.metadata import Metadata

from tool.fastq_splitter import fastq_splitter


@pytest.mark.wgbs
def test_splitter_00():
    """
    Extract the compressed FASTQ files
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    fastq_file_1 = resource_path + "bsSeeker.Mouse.SRR892982_1.fastq"
    fastq_file_2 = resource_path + "bsSeeker.Mouse.SRR892982_2.fastq"

    files = [fastq_file_1, fastq_file_2]

    for fastq_gz in files:
        with gzip.open(fastq_gz + '.gz', 'rb') as fgz_in:
            with open(fastq_gz, 'wb') as f_out:
                f_out.write(fgz_in.read())

    assert os.path.isfile(fastq_file_1) is True
    assert os.path.getsize(fastq_file_1) > 0
    assert os.path.isfile(fastq_file_2) is True
    assert os.path.getsize(fastq_file_2) > 0


@pytest.mark.wgbs
def test_splitter_01():
    """
    Function to test paired splitter
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    fastq_1file = resource_path + "bsSeeker.Mouse.SRR892982_1.fastq"
    fastq_2file = resource_path + "bsSeeker.Mouse.SRR892982_2.fastq"

    fqs_handle = fastq_splitter()
    results = fqs_handle.run(
        {
            "fastq1": fastq_1file,
            "fastq2": fastq_2file
        },
        {
            "fastq1": Metadata(
                "data_rnaseq", "fastq", [], None,
                {'assembly': 'test'}),
            "fastq2": Metadata(
                "data_rnaseq", "fastq", [], None,
                {'assembly': 'test'})
        },
        {"output": fastq_1file + ".tar.gz"}
    )

    print("WGBS - PAIRED RESULTS:", results)

    assert os.path.isfile(results[0]["output"]) is True
    assert os.path.getsize(results[0]["output"]) > 0

    os.remove(results[0]["output"])


@pytest.mark.wgbs
def test_splitter_02():
    """
    Function to test single splitter
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    fastq_file = resource_path + "bsSeeker.Mouse.SRR892982_1.fastq"

    fqs_handle = fastq_splitter()
    results = fqs_handle.run(
        {
            "fastq1": fastq_file,
        },
        {
            "fastq1": Metadata(
                "data_rnaseq", "fastq", [], None,
                {'assembly': 'test'})
        },
        {"output": fastq_file + ".tar.gz"}
    )

    print("WGBS - SINGLE RESULTS:", results)

    assert os.path.isfile(results[0]["output"]) is True
    assert os.path.getsize(results[0]["output"]) > 0

    os.remove(results[0]["output"])
