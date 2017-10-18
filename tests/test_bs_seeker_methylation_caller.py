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

from tool import bs_seeker_methylation_caller

@pytest.mark.wgbs
def test_bs_seeker_methylation_caller():
    """
    Test that it is possible to call the methylation called by BS seeker
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genome_fa_file = resource_path + "bsSeeker.Mouse.GRCm38.fasta_bowtie2.tar.gz"
    bam_file = resource_path + "bsSeeker.Mouse.GRCm38.bam"
    home = os.path.expanduser('~')

    bsmc = bs_seeker_methylation_caller.bssMethylationCallerTool()
    bs_files, bs_meta = bsmc.run(
        [
            bam_file,
            genome_fa_file
        ],
        [],
        {
            'bss_path': home + "/lib/BSseeker2"
        }
    )

    assert os.path.isfile(bs_files[0]) is True
    assert os.path.getsize(bs_files[0]) > 0
    assert os.path.isfile(bs_files[1]) is True
    # Blank file for this small dataset
    #assert os.path.getsize(bs_files[1]) > 0
    assert os.path.isfile(bs_files[2]) is True
    assert os.path.getsize(bs_files[2]) > 0
