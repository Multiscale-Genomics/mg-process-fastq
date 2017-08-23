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

from __future__ import print_function

import os.path
import pytest # pylint: disable=unused-import

from tool.tb_save_hdf5_matrix import tbSaveAdjacencyHDF5Tool

@pytest.mark.hic
def test_tb_save_hdf5_matrix_frag_01():
    """
    Test case to ensure that the BWA indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    reads_tsv = resource_path + "tb.Human.SRR1658573_frag_01_filtered_map.tsv"
    genome_gem_fa = resource_path + "tb.Human.GCA_000001405.22_gem.fasta"

    metadata = {
        'assembly' : 'test',
        'expt_name' : 'tb.Human.SRR1658573_frag_01',
        'enzyme_name' : 'MboI',
        'windows' : ((1, 'end')),
        'mapping' : ['frag', 'frag'],
        'resolutions' : [10000, 100000],
        'normalized' : False
    }

    tgt_handle = tbSaveAdjacencyHDF5Tool()
    tgt_files, tgt_meta = tgt_handle.run([reads_tsv, genome_gem_fa], [], metadata)

    print(tgt_files)

    assert os.path.isfile(tgt_files[0]) is True
    assert os.path.getsize(tgt_files[0]) > 0

@pytest.mark.hic
def test_tb_save_hdf5_matrix_frag_02():
    """
    Test case to ensure that the BWA indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    reads_tsv = resource_path + "tb.Human.SRR1658573_frag_02_filtered_map.tsv"
    genome_gem_fa = resource_path + "tb.Human.GCA_000001405.22_gem.fasta"

    metadata = {
        'assembly' : 'test',
        'expt_name' : 'tb.Human.SRR1658573_frag_02',
        'enzyme_name' : 'MboI',
        'windows' : ((1, 'end')),
        'mapping' : ['frag', 'frag'],
        'resolutions' : [10000, 100000],
        'normalized' : False
    }

    tgt_handle = tbSaveAdjacencyHDF5Tool()
    tgt_files, tgt_meta = tgt_handle.run([reads_tsv, genome_gem_fa], [], metadata)

    print(tgt_files)

    assert os.path.isfile(tgt_files[0]) is True
    assert os.path.getsize(tgt_files[0]) > 0

@pytest.mark.hic
def test_tb_save_hdf5_matrix_iter_01():
    """
    Test case to ensure that the BWA indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    reads_tsv = resource_path + "tb.Human.SRR1658573_iter_01_filtered_map.tsv"
    genome_gem_fa = resource_path + "tb.Human.GCA_000001405.22_gem.fasta"

    metadata = {
        'assembly' : 'test',
        'expt_name' : 'tb.Human.SRR1658573_frag_01',
        'enzyme_name' : 'MboI',
        'windows' : ((1, 'end')),
        'mapping' : ['iter', 'iter'],
        'resolutions' : [10000, 100000],
        'normalized' : False
    }

    tgt_handle = tbSaveAdjacencyHDF5Tool()
    tgt_files, tgt_meta = tgt_handle.run([reads_tsv, genome_gem_fa], [], metadata)

    print(tgt_files)

    assert os.path.isfile(tgt_files[0]) is True
    assert os.path.getsize(tgt_files[0]) > 0

@pytest.mark.hic
def test_tb_save_hdf5_matrix_iter_02():
    """
    Test case to ensure that the BWA indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    reads_tsv = resource_path + "tb.Human.SRR1658573_iter_02_filtered_map.tsv"
    genome_gem_fa = resource_path + "tb.Human.GCA_000001405.22_gem.fasta"

    metadata = {
        'assembly' : 'test',
        'expt_name' : 'tb.Human.SRR1658573_iter_02',
        'enzyme_name' : 'MboI',
        'windows' : ((1, 'end')),
        'mapping' : ['iter', 'iter'],
        'resolutions' : [10000, 100000],
        'normalized' : False
    }

    tgt_handle = tbSaveAdjacencyHDF5Tool()
    tgt_files, tgt_meta = tgt_handle.run([reads_tsv, genome_gem_fa], [], metadata)

    print(tgt_files)

    assert os.path.isfile(tgt_files[0]) is True
    assert os.path.getsize(tgt_files[0]) > 0
