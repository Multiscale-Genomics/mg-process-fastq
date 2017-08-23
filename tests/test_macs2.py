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

from tool.macs2 import macs2

@pytest.mark.chipseq
def test_macs2():
    """
    Function to test MACS2
    """
    macs_handle = macs2()
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    bam_file_stub = resource_path + "macs2.Human.DRR000150.22.filtered"
    bam_file = bam_file_stub + '.bam'

    macs_handle.run([bam_file], [])

    assert os.path.isfile(bam_file_stub + '_peaks.narrowPeak') is True
    assert os.path.getsize(bam_file_stub + '_peaks.narrowPeak') > 0
    assert os.path.isfile(bam_file_stub + '_summits.bed') is True
    assert os.path.getsize(bam_file_stub + '_summits.bed') > 0
