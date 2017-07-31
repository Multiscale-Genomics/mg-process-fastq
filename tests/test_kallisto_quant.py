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
from __future__ import print_function

import os.path
import pytest # pylint: disable=unused-import

from tool.kallisto_quant import kallistoQuantificationTool

@pytest.mark.ranseq
def test_kallisto_quant():
    """
    Function to test Kallisto quantifier
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    kqft = kallistoQuantificationTool()
    fastq1 = resource_path + "kallisto.Human.ERR030872_1.fastq"
    fastq2 = resource_path + "kallisto.Human.ERR030872_2.fastq"
    kqft.run(
        [resource_path + "kallisto.Human.GRCh38.idx", fastq1, fastq2], {}, )

    print(__file__)

    assert os.path.isfile(resource_path + "abundance.h5") is True
    assert os.path.getsize(resource_path + "abundance.h5") > 0
    assert os.path.isfile(resource_path + "abundance.tsv") is True
    assert os.path.getsize(resource_path + "abundance.tsv") > 0
    assert os.path.isfile(resource_path + "run_info.json") is True
    assert os.path.getsize(resource_path + "run_info.json") > 0
