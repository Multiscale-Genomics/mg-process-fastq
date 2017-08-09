"""
.. Copyright 2017 EMBL-European Bioinformatics Institute

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

from tool import bs_seeker_filter

@pytest.mark.wgbs
def test_bs_seeker_filter_01():
    """
    Test that it is possible to call the BSseeker filter
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genomefa_file = resource_path + "bsSeeker.Mouse.GRCm38_1.fastq"
    home = os.path.expanduser('~')

    bsi = bs_seeker_filter.filterReadsTool()
    bs_files, bs_meta = bsi.run(
        [genomefa_file],
        [],
        {
            "aligner" : "bowtie",
            "aligner_path" : home + "/bin",
            "bss_path" : home + "/lib/BSseeker2"
        }
    )

    assert os.path.isfile(bs_files[0]) is True
    assert os.path.getsize(bs_files[0]) > 0

@pytest.mark.wgbs
def test_bs_seeker_filter_02():
    """
    Test that it is possible to call the BSseeker filter
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genomefa_file = resource_path + "bsSeeker.Mouse.GRCm38_2.fastq"
    home = os.path.expanduser('~')

    bsi = bs_seeker_filter.filterReadsTool()
    bs_files, bs_meta = bsi.run(
        [genomefa_file],
        [],
        {
            "aligner" : "bowtie",
            "aligner_path" : home + "/bin",
            "bss_path" : home + "/lib/BSseeker2"
        }
    )

    assert os.path.isfile(bs_files[0]) is True
    assert os.path.getsize(bs_files[0]) > 0
