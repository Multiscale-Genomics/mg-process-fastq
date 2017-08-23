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

from tool import biobambam_filter

@pytest.mark.chipseq
def test_biobambam():
    """
    Test case to ensure that BioBamBam works
    """
    bbb = biobambam_filter.biobambam()
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    bbb.run(
        [resource_path + "macs2.Human.DRR000150.22.bam"],
        []
    )

    assert os.path.isfile(resource_path + "macs2.Human.DRR000150.22.filtered.bam") is True
    assert os.path.getsize(resource_path + "macs2.Human.DRR000150.22.filtered.bam") > 0
