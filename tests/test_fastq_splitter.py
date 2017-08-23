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

from tool.fastq_splitter import fastq_splitter

@pytest.mark.wgbs
def test_paired_splitter():
    """
    Function to test paired splitter
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    fastq_1file = resource_path + "bsSeeker.Mouse.GRCm38_1.fastq"
    fastq_2file = resource_path + "bsSeeker.Mouse.GRCm38_2.fastq"

    fqs_handle = fastq_splitter()
    results = fqs_handle.run([fastq_1file, fastq_2file], [], {})

    print("WGBS - PAIRED RESULTS:", results)

    assert os.path.isfile(results[0]) is True
    assert os.path.getsize(results[0]) > 0


@pytest.mark.wgbs
def test_single_splitter():
    """
    Function to test single splitter
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    fastq_2file = resource_path + "bsSeeker.Mouse.GRCm38_2.fastq"

    fqs_handle = fastq_splitter()
    results = fqs_handle.run([fastq_2file], [], {})

    print("WGBS - SINGLE RESULTS:", results)

    assert os.path.isfile(results[0]) is True
    assert os.path.getsize(results[0]) > 0
