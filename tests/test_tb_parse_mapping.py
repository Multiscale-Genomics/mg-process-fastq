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

from tool.tb_parse_mapping import tbParseMappingTool

@pytest.mark.hic
def test_tb_parse_mapping_frag():
    """
    Test case to ensure that the BWA indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genome_fa = resource_path + "tb.Human.GCA_000001405.22_gem.fasta"

    map_frag_1 = resource_path + "tb.Human.SRR1658573_1_frag.map"
    map_full_1 = resource_path + "tb.Human.SRR1658573_1_full.map"

    map_frag_2 = resource_path + "tb.Human.SRR1658573_2_frag.map"
    map_full_2 = resource_path + "tb.Human.SRR1658573_2_full.map"

    files = [
        genome_fa,
        map_full_1, map_frag_1,
        map_full_2, map_frag_2
    ]

    metadata = {
        'expt_name' : 'tb.Human.SRR1658573',
        'enzyme_name' : 'MboI',
        'mapping' : ['frag', 'frag']
    }

    tpm = tbParseMappingTool()
    tpm_files, tpm_meta = tpm.run(files, [], metadata)

    output_frag = resource_path + "tb.Human.SRR1658573_frag.tsv"

    assert os.path.isfile(output_frag) is True
    assert os.path.getsize(output_frag) > 0

@pytest.mark.hic
def test_tb_parse_mapping_iter():
    """
    Test case to ensure that the BWA indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genome_fa = resource_path + "tb.Human.GCA_000001405.22_gem.fasta"

    map25_1 = resource_path + "tb.Human.SRR1658573_1_full_1-25.map"
    map50_1 = resource_path + "tb.Human.SRR1658573_1_full_1-50.map"
    map75_1 = resource_path + "tb.Human.SRR1658573_1_full_1-75.map"
    map100_1 = resource_path + "tb.Human.SRR1658573_1_full_1-100.map"

    map25_2 = resource_path + "tb.Human.SRR1658573_2_full_1-25.map"
    map50_2 = resource_path + "tb.Human.SRR1658573_2_full_1-50.map"
    map75_2 = resource_path + "tb.Human.SRR1658573_2_full_1-75.map"
    map100_2 = resource_path + "tb.Human.SRR1658573_2_full_1-100.map"

    files = [
        genome_fa,
        map25_1, map50_1, map75_1, map100_1,
        map25_2, map50_2, map75_2, map100_2,
    ]

    metadata = {
        'assembly' : 'test',
        'expt_name' : 'tb.Human.SRR1658573',
        'enzyme_name' : 'MboI',
        'windows' : ((1, 25), (1, 50), (1, 75), (1, 100)),
        'mapping' : ['iter', 'iter']
    }

    tpm = tbParseMappingTool()
    tpm_files, tpm_meta = tpm.run(files, [], metadata)

    output_iter = resource_path + "tb.Human.SRR1658573_iter.tsv"

    assert os.path.isfile(output_iter) is True
    assert os.path.getsize(output_iter) > 0
