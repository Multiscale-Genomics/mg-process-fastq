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

from tool.tb_filter import tbFilterTool

@pytest.mark.hic
def test_tb_filter_frag_01():
    """
    Test case to ensure that the BWA indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    reads_tsv = resource_path + "tb.Human.SRR1658573_frag.tsv"

    metadata = {
        'assembly' : 'test',
        'expt_name' : 'tb.Human.SRR1658573_frag_01',
        'enzyme_name' : 'MboI',
        'windows' : ((1, 'end')),
        'mapping' : ['frag', 'frag']
    }

    tpm = tbFilterTool()
    tpm_files, tpm_meta = tpm.run([reads_tsv], [], metadata)

    reads_tsv = resource_path + metadata['expt_name'] + "_filtered_map.tsv"

    assert os.path.isfile(resource_path + "tb.Human.SRR1658573_frag_01_filtered_map.tsv") is True
    assert os.path.getsize(resource_path + "tb.Human.SRR1658573_frag_01_filtered_map.tsv") > 0
    assert os.path.isfile(reads_tsv + '_dangling-end.tsv') is True
    assert os.path.getsize(reads_tsv + '_dangling-end.tsv') > 0
    assert os.path.isfile(reads_tsv + '_duplicated.tsv') is True
    assert os.path.getsize(reads_tsv + '_duplicated.tsv') > 0
    assert os.path.isfile(reads_tsv + '_error.tsv') is True
    assert os.path.getsize(reads_tsv + '_error.tsv') > 0
    assert os.path.isfile(reads_tsv + '_extra_dangling-end.tsv') is True
    assert os.path.getsize(reads_tsv + '_extra_dangling-end.tsv') > 0
    assert os.path.isfile(reads_tsv + '_over-represented.tsv') is True
    assert os.path.getsize(reads_tsv + '_over-represented.tsv') > 0
    assert os.path.isfile(reads_tsv + '_random_breaks.tsv') is True
    assert os.path.getsize(reads_tsv + '_random_breaks.tsv') > 0
    assert os.path.isfile(reads_tsv + '_self-circle.tsv') is True
    assert os.path.getsize(reads_tsv + '_self-circle.tsv') > 0
    assert os.path.isfile(reads_tsv + '_too_close_from_RES.tsv') is True
    assert os.path.getsize(reads_tsv + '_too_close_from_RES.tsv') > 0
    assert os.path.isfile(reads_tsv + '_too_large.tsv') is True
    #assert os.path.getsize(reads_tsv + '_too_large.tsv') > 0
    assert os.path.isfile(reads_tsv + '_too_short.tsv') is True
    assert os.path.getsize(reads_tsv + '_too_short.tsv') > 0

@pytest.mark.hic
def test_tb_filter_frag_02():
    """
    Test case to ensure that the BWA indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    reads_tsv = resource_path + "tb.Human.SRR1658573_frag.tsv"

    metadata = {
        'assembly' : 'test',
        'expt_name' : 'tb.Human.SRR1658573_frag_02',
        'enzyme_name' : 'MboI',
        'windows' : ((1, 'end')),
        'mapping' : ['frag', 'frag'],
        'conservative_filtering' : True
    }

    tpm = tbFilterTool()
    tpm_files, tpm_meta = tpm.run([reads_tsv], [], metadata)

    reads_tsv = resource_path + metadata['expt_name'] + "_filtered_map.tsv"

    assert os.path.isfile(resource_path + "tb.Human.SRR1658573_frag_02_filtered_map.tsv") is True
    assert os.path.getsize(resource_path + "tb.Human.SRR1658573_frag_02_filtered_map.tsv") > 0
    assert os.path.isfile(reads_tsv + '_dangling-end.tsv') is True
    assert os.path.getsize(reads_tsv + '_dangling-end.tsv') > 0
    assert os.path.isfile(reads_tsv + '_duplicated.tsv') is True
    assert os.path.getsize(reads_tsv + '_duplicated.tsv') > 0
    assert os.path.isfile(reads_tsv + '_error.tsv') is True
    assert os.path.getsize(reads_tsv + '_error.tsv') > 0
    assert os.path.isfile(reads_tsv + '_extra_dangling-end.tsv') is True
    assert os.path.getsize(reads_tsv + '_extra_dangling-end.tsv') > 0
    assert os.path.isfile(reads_tsv + '_over-represented.tsv') is True
    assert os.path.getsize(reads_tsv + '_over-represented.tsv') > 0
    assert os.path.isfile(reads_tsv + '_random_breaks.tsv') is True
    assert os.path.getsize(reads_tsv + '_random_breaks.tsv') > 0
    assert os.path.isfile(reads_tsv + '_self-circle.tsv') is True
    assert os.path.getsize(reads_tsv + '_self-circle.tsv') > 0
    assert os.path.isfile(reads_tsv + '_too_close_from_RES.tsv') is True
    assert os.path.getsize(reads_tsv + '_too_close_from_RES.tsv') > 0
    assert os.path.isfile(reads_tsv + '_too_large.tsv') is True
    #assert os.path.getsize(reads_tsv + '_too_large.tsv') > 0
    assert os.path.isfile(reads_tsv + '_too_short.tsv') is True
    assert os.path.getsize(reads_tsv + '_too_short.tsv') > 0

@pytest.mark.hic
def test_tb_filter_iter_01():
    """
    Test case to ensure that the BWA indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    reads_tsv = resource_path + "tb.Human.SRR1658573_iter.tsv"

    metadata = {
        'assembly' : 'test',
        'expt_name' : 'tb.Human.SRR1658573_iter_01',
        'enzyme_name' : 'MboI',
        'windows' : ((1, 'end')),
        'mapping' : ['iter', 'iter']
    }

    tpm = tbFilterTool()
    tpm_files, tpm_meta = tpm.run([reads_tsv], [], metadata)

    reads_tsv = resource_path + metadata['expt_name'] + "_filtered_map.tsv"

    assert os.path.isfile(resource_path + "tb.Human.SRR1658573_iter_01_filtered_map.tsv") is True
    assert os.path.getsize(resource_path + "tb.Human.SRR1658573_iter_01_filtered_map.tsv") > 0
    assert os.path.isfile(reads_tsv + '_dangling-end.tsv') is True
    assert os.path.getsize(reads_tsv + '_dangling-end.tsv') > 0
    assert os.path.isfile(reads_tsv + '_duplicated.tsv') is True
    assert os.path.getsize(reads_tsv + '_duplicated.tsv') > 0
    assert os.path.isfile(reads_tsv + '_error.tsv') is True
    assert os.path.getsize(reads_tsv + '_error.tsv') > 0
    assert os.path.isfile(reads_tsv + '_extra_dangling-end.tsv') is True
    assert os.path.getsize(reads_tsv + '_extra_dangling-end.tsv') > 0
    assert os.path.isfile(reads_tsv + '_over-represented.tsv') is True
    assert os.path.getsize(reads_tsv + '_over-represented.tsv') > 0
    assert os.path.isfile(reads_tsv + '_random_breaks.tsv') is True
    assert os.path.getsize(reads_tsv + '_random_breaks.tsv') > 0
    assert os.path.isfile(reads_tsv + '_self-circle.tsv') is True
    assert os.path.getsize(reads_tsv + '_self-circle.tsv') > 0
    assert os.path.isfile(reads_tsv + '_too_close_from_RES.tsv') is True
    assert os.path.getsize(reads_tsv + '_too_close_from_RES.tsv') > 0
    assert os.path.isfile(reads_tsv + '_too_large.tsv') is True
    #assert os.path.getsize(reads_tsv + '_too_large.tsv') > 0
    assert os.path.isfile(reads_tsv + '_too_short.tsv') is True
    assert os.path.getsize(reads_tsv + '_too_short.tsv') > 0

@pytest.mark.hic
def test_tb_filter_iter_02():
    """
    Test case to ensure that the BWA indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    reads_tsv = resource_path + "tb.Human.SRR1658573_iter.tsv"

    metadata = {
        'assembly' : 'test',
        'expt_name' : 'tb.Human.SRR1658573_iter_02',
        'enzyme_name' : 'MboI',
        'windows' : ((1, 'end')),
        'mapping' : ['iter', 'iter'],
        'conservative_filtering' : True
    }

    tpm = tbFilterTool()
    tpm_files, tpm_meta = tpm.run([reads_tsv], [], metadata)

    reads_tsv = resource_path + metadata['expt_name'] + "_filtered_map.tsv"

    assert os.path.isfile(resource_path + "tb.Human.SRR1658573_iter_02_filtered_map.tsv") is True
    assert os.path.getsize(resource_path + "tb.Human.SRR1658573_iter_02_filtered_map.tsv") > 0
    assert os.path.isfile(reads_tsv + '_dangling-end.tsv') is True
    assert os.path.getsize(reads_tsv + '_dangling-end.tsv') > 0
    assert os.path.isfile(reads_tsv + '_duplicated.tsv') is True
    assert os.path.getsize(reads_tsv + '_duplicated.tsv') > 0
    assert os.path.isfile(reads_tsv + '_error.tsv') is True
    assert os.path.getsize(reads_tsv + '_error.tsv') > 0
    assert os.path.isfile(reads_tsv + '_extra_dangling-end.tsv') is True
    assert os.path.getsize(reads_tsv + '_extra_dangling-end.tsv') > 0
    assert os.path.isfile(reads_tsv + '_over-represented.tsv') is True
    assert os.path.getsize(reads_tsv + '_over-represented.tsv') > 0
    assert os.path.isfile(reads_tsv + '_random_breaks.tsv') is True
    assert os.path.getsize(reads_tsv + '_random_breaks.tsv') > 0
    assert os.path.isfile(reads_tsv + '_self-circle.tsv') is True
    assert os.path.getsize(reads_tsv + '_self-circle.tsv') > 0
    assert os.path.isfile(reads_tsv + '_too_close_from_RES.tsv') is True
    assert os.path.getsize(reads_tsv + '_too_close_from_RES.tsv') > 0
    assert os.path.isfile(reads_tsv + '_too_large.tsv') is True
    #assert os.path.getsize(reads_tsv + '_too_large.tsv') > 0
    assert os.path.isfile(reads_tsv + '_too_short.tsv') is True
    assert os.path.getsize(reads_tsv + '_too_short.tsv') > 0
