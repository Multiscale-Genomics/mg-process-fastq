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

from tool.tb_full_mapping import tbFullMappingTool

def test_tb_full_mapping():
    """
    Test case to ensure that the BWA indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    gem_file = resource_path + "tb.Human.GCA_000001405.22.fasta.gem"

    fastq_file_1 = resource_path + "tb.Human.SRR1658573_1.fastq"

    files = [
        gem_file,
        fastq_file_1
    ]

    metadata = {
        'assembly' : 'test',
        'windows' : ((1, 25), (1, 50), (1, 75), (1, 100))
    }

    gem_file = files[1]

    print(gem_file)

    tfm1 = tbFullMappingTool()
    tfm1_files, tfm1_meta = tfm1.run(files, [], metadata)


    #assert os.path.isfile(gem_file) is True
    #assert os.path.getsize(gem_file) > 0
