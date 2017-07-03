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

import random
import os.path
import pytest

from tool import macs2


def test_macs2():
    """
    Function to test MACS2
    """
    m = macs2.macs2()
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    bam_file = resource_path + "biobambam.Human.DRR000150.22_output.bam"
    if not os.path.isfile(bam_file):
        bam_file = resource_path + "macs2.Human.DRR000150.22.bam"

    summits_bed = resource_path + "_summits.bed"
    narrow_peak = resource_path + "_narrowPeak"
    broad_peak = resource_path + "_broadPeak"
    gapped_peak = resource_path + "_gappedPeak"

    m.run([bam_file], {}, [summits_bed, narrow_peak, broad_peak, gapped_peak])
