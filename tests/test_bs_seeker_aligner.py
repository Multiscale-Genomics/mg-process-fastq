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
import pytest # pylint: disable=unused-import

from basic_modules.metadata import Metadata

from tool import bs_seeker_aligner

@pytest.mark.wgbs
def test_bs_seeker_aligner():
    """
    Test to ensure bs-Seeker aligner works
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    home = os.path.expanduser('~')

    input_files = {
        "genome": resource_path + "bsSeeker.Mouse.GRCm38.fasta",
        "index": resource_path + "bsSeeker.Mouse.GRCm38.fasta.bt2.tar.gz",
        "fastq1": resource_path + "bsSeeker.Mouse.GRCm38_1_filtered.fastq",
        "fastq2": resource_path + "bsSeeker.Mouse.GRCm38_2_filtered.fastq",
    }

    output_files = {
        "bam": resource_path + "bsSeeker.Mouse.GRCm38_1_filtered.bam",
        "bai": resource_path + "bsSeeker.Mouse.GRCm38_1_filtered.bai"
    }

    metadata = {
        "genome": Metadata(
            "Assembly", "fasta", input_files["genome"], None,
            {'assembly' : 'test'}),
        "index": Metadata(
            "index_bowtie", "index", input_files["genome"], None,
            {'assembly' : 'test'}),
        "fastq1": Metadata(
            "data_wgbs", "fastq", input_files["fastq1"], None,
            {'assembly' : 'test'}),
        "fastq2": Metadata(
            "data_wgbs", "fastq", input_files["fastq2"], None,
            {'assembly' : 'test'}),
        "aligner" : "bowtie2",
        "aligner_path" : home + "/lib/bowtie2-2.3.2",
        "bss_path" : home + "/lib/BSseeker2"
    }

    bsa = bs_seeker_aligner.bssAlignerTool()
    bsa.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["bam"]) is True
    assert os.path.getsize(output_files["bam"]) > 0
    assert os.path.isfile(output_files["bai"]) is True
    assert os.path.getsize(output_files["bai"]) > 0
