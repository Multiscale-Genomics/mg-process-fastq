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
import pytest

from tool import inps

@pytest.mark.py3
@pytest.mark.mnaseseq
def test_inps():
    """
    Function to test INPS works.
    """

    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    bam_file = resource_path + "inps.Mouse.DRR000386.bam"
    peak_bed = bam_file.replace('.bam', '.bed')

    inps_obj = inps.inps()
    inps_obj.run([bam_file, peak_bed], {})

    assert os.path.isfile(resource_path + "inps.Mouse.DRR000386.bam.bed") is True
    assert os.path.getsize(resource_path + "inps.Mouse.DRR000386.bam.bed") > 0
    assert os.path.isfile(resource_path + "inps.Mouse.DRR000386.bed_19.like_bed") is True
    assert os.path.getsize(resource_path + "inps.Mouse.DRR000386.bed_19.like_bed") > 0
    assert os.path.isfile(resource_path + "inps.Mouse.DRR000386.bed_19.like_wig") is True
    assert os.path.getsize(resource_path + "inps.Mouse.DRR000386.bed_19.like_wig") > 0

    # All results final file
    assert os.path.isfile(resource_path + "inps.Mouse.DRR000386.bed_Gathering.like_bed") is True
    assert os.path.getsize(resource_path + "inps.Mouse.DRR000386.bed_Gathering.like_bed") > 0
