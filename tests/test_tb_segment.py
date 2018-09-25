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

from tool.tb_segment import tbSegmentTool


@pytest.mark.hic
def test_tb_segment():
    """
    Test case to detect tads and compartments in a sample bam file
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    bamin = resource_path + "tb_paired_reads.bam"
    hic_biases = resource_path + "tb_biases_100kb.pickle"

    files = [bamin, hic_biases]

    metadata = {
        'resolution': "100000",
        'chromosome_names': 'chr21',
        "callers": [
            "1",
            "2"
        ],
        'ncpus': 4
    }

    ts_handle = tbSegmentTool()
    ts_files, ts_meta = ts_handle.run(files, [], metadata)

    assert len(ts_files) == 2
    assert len(ts_meta) == 2
