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

import process_wgbs  # from mg-process-fastq

@pytest.mark.wgbs
def test_single_splitter():
    """
    Function to test single splitter
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    fastq_1file = resource_path + "bsSeeker.Mouse.GRCm38_1.fastq"

    ssp = process_wgbs.process_wgbs()
    ssp.single_splitter(fastq_1file)
