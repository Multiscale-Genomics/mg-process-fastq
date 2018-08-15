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

import os
import gzip
import pytest

from basic_modules.metadata import Metadata

from tool import bs_seeker_filter


@pytest.mark.wgbs
def test_bs_seeker_filter_00():
    """
    Extract the compressed FASTQ files
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    fastq_file_1 = resource_path + "bsSeeker.Mouse.SRR892982_1.fastq"
    fastq_file_2 = resource_path + "bsSeeker.Mouse.SRR892982_2.fastq"

    with gzip.open(fastq_file_1 + '.gz', 'rb') as fgz_in:
        with open(fastq_file_1, 'w') as f_out:
            f_out.write(fgz_in.read())

    with gzip.open(fastq_file_2 + '.gz', 'rb') as fgz_in:
        with open(fastq_file_2, 'w') as f_out:
            f_out.write(fgz_in.read())

    assert os.path.isfile(fastq_file_1) is True
    assert os.path.getsize(fastq_file_1) > 0
    assert os.path.isfile(fastq_file_2) is True
    assert os.path.getsize(fastq_file_2) > 0


@pytest.mark.wgbs
def test_bs_seeker_filter_01():
    """
    Test that it is possible to call the BSseeker filter
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    home = os.path.expanduser('~')

    input_files = {
        "fastq": resource_path + "bsSeeker.Mouse.SRR892982_1.fastq"
    }

    output_files = {
        "fastq_filtered": resource_path + "bsSeeker.Mouse.SRR892982_1.filtered.fastq"
    }

    metadata = {
        "fastq": Metadata(
            "data_wgbs", "fastq", input_files["fastq"], None,
            {'assembly': 'test'})
    }

    config_param = {
        "aligner": "bowtie2",
        "aligner_path": home + "/lib/bowtie2-2.3.4-linux-x86_64",
        "bss_path": home + "/lib/BSseeker2",
        "execution": resource_path
    }

    bsi = bs_seeker_filter.filterReadsTool(config_param)
    bsi.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["fastq_filtered"]) is True
    assert os.path.getsize(output_files["fastq_filtered"]) > 0


@pytest.mark.wgbs
def test_bs_seeker_filter_02():
    """
    Test that it is possible to call the BSseeker filter
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    home = os.path.expanduser('~')

    input_files = {
        "fastq": resource_path + "bsSeeker.Mouse.SRR892982_2.fastq"
    }

    output_files = {
        "fastq_filtered": resource_path + "bsSeeker.Mouse.SRR892982_2.filtered.fastq"
    }

    metadata = {
        "fastq": Metadata(
            "data_wgbs", "fastq", input_files["fastq"], None,
            {'assembly': 'test'})
    }

    config_param = {
        "aligner": "bowtie2",
        "aligner_path": home + "/lib/bowtie2-2.3.4-linux-x86_64",
        "bss_path": home + "/lib/BSseeker2",
        "execution": resource_path
    }

    bsi = bs_seeker_filter.filterReadsTool(config_param)
    bsi.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["fastq_filtered"]) is True
    assert os.path.getsize(output_files["fastq_filtered"]) > 0
