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

from tool import bs_seeker_aligner

@pytest.mark.wgbs
def test_bs_seeker_aligner():
    """
    Test to ensure bs-Seeker aligner works
    """

    home = os.path.expanduser('~')
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genomefa_file = resource_path + "bsSeeker.Mouse.GRCm38.fasta"
    genomeidx_file = resource_path + "bsSeeker.Mouse.GRCm38.fasta_bowtie2.tar.gz"
    fastq1_file = resource_path + "bsSeeker.Mouse.GRCm38_1.fastq"
    fastq2_file = resource_path + "bsSeeker.Mouse.GRCm38_2.fastq"

    bsa = bs_seeker_aligner.bssAlignerTool()
    bs_files, bs_meta = bsa.run(
        [
            genomefa_file,
            genomeidx_file,
            fastq1_file,
            fastq2_file,
        ],
        [],
        {
            "aligner" : "bowtie2",
            "aligner_path" : home + "/lib/bowtie2-2.3.2",
            "bss_path" : home + "/lib/BSseeker2"
        }
    )

    assert os.path.isfile(bs_files[0]) is True
    assert os.path.getsize(bs_files[0]) > 0
