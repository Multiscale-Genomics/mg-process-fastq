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
import pysam
import pytest # pylint: disable=unused-import

from tool import bs_seeker_methylation_caller

@pytest.mark.wgbs
def test_bs_seeker_methylation_caller():
    """
    Test that it is possible to call the methylation called by BS seeker
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genome_fa_file = resource_path + "bsSeeker.Mouse.GRCm38.fasta_bowtie2"
    bam_file = resource_path + "bsSeeker.Mouse.GRCm38.bam"
    home = os.path.expanduser('~')

    pysam.sort("-o", str(bam_file+".sorted"), str(bam_file))

    bsmc = bs_seeker_methylation_caller.bssMethylationCallerTool()
    bsmc.run(
        [bam_file+".sorted"],
        [],
        {
            'bss_path': home + "/lib/BSseeker2",
            'index_path' : genome_fa_file
        }
    )
