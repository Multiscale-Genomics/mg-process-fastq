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
from pysam import AlignmentFile  # pylint: disable=no-name-in-module
import pytest

from mg_process_fastq.tool.tb_normalize import tbNormalizeTool


@pytest.mark.hic
def test_tb_normalize():
    """
    Test case to normalize a sample bam file to 100K
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    bamin = resource_path + "tb_paired_reads.bam"

    metadata = {
        'resolution': "100000",
        'min_perc': '2',
        'max_perc': '99.8',
        'ncpus': '4'
    }

    bamfile = AlignmentFile(bamin, 'rb')
    if len(bamfile.references) == 1:
        metadata["min_count"] = "10"
    bamfile.close()

    tn_handle = tbNormalizeTool()
    tn_files, tn_meta = tn_handle.run([bamin], [], metadata)

    assert len(tn_files) > 1
    assert len(tn_meta) > 1
