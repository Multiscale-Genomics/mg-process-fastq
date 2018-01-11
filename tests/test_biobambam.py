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
import pytest # pylint: disable=unused-import

from basic_modules.metadata import Metadata

from tool import biobambam_filter

@pytest.mark.chipseq
def test_biobambam_chipseq():
    """
    Test case to ensure that BioBamBam works
    """

    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "input": resource_path + "macs2.Human.DRR000150.22.bam"
    }

    output_files = {
        "output": resource_path + "macs2.Human.DRR000150.22_filtered.bam"
    }

    metadata = {
        "input": Metadata(
            "data_chipseq", "fastq", [], None,
            {'assembly' : 'test'}),
    }

    bbb = biobambam_filter.biobambam()
    bbb.run(input_files, metadata, output_files)

    assert os.path.isfile(resource_path + "macs2.Human.DRR000150.22_filtered.bam") is True
    assert os.path.getsize(resource_path + "macs2.Human.DRR000150.22_filtered.bam") > 0

@pytest.mark.idamidseq
def test_biobambam_idamidseq():
    """
    Test case to ensure that BioBamBam works
    """

    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    bam_files = [
        resource_path + "idear.Human.SRR3714775.bam",
        resource_path + "idear.Human.SRR3714776.bam",
        resource_path + "idear.Human.SRR3714777.bam",
        resource_path + "idear.Human.SRR3714778.bam"
    ]

    for bam_file in bam_files:
        input_files = {
            "input": bam_file
        }

        output_files = {
            "output": bam_file.replace(".bam", "_filtered.bam")
        }

        metadata = {
            "input": Metadata(
                "data_damid_seq", "bam", [], None,
                {'assembly' : 'test'}),
        }

        bbb = biobambam_filter.biobambam()
        bbb.run(input_files, metadata, output_files)

        assert os.path.isfile(bam_file.replace(".bam", "_filtered.bam")) is True
        assert os.path.getsize(bam_file.replace(".bam", "_filtered.bam")) > 0
