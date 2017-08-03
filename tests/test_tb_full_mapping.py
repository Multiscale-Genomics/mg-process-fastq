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
import gzip
import pytest # pylint: disable=unused-import

from tool.tb_full_mapping import tbFullMappingTool

#@pytest.mark.hic
def test_tb_extract_fastq():
    """
    Extract the compressed FASTQ files
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    fastq_file_1 = resource_path + "tb.Human.SRR1658573_1.fastq"
    fastq_file_2 = resource_path + "tb.Human.SRR1658573_2.fastq"

    with gzip.open(fastq_file_1 + '.gz', 'rb') as fgz_in:
        with open(fastq_file_1, 'w') as f_out:
            f_out.write(fgz_in.read())

    with gzip.open(fastq_file_2 + '.gz', 'rb') as fgz_in:
        with open(fastq_file_2, 'w') as f_out:
            f_out.write(fgz_in.read())

    assert os.path.isfile(fastq_file_1) is True
    assert os.path.getsize(fastq_file_1) > 0
    assert os.path.isfile(fastq_file_2) is True
    assert os.path.getsize(fastq_file_2) > 0

@pytest.mark.hic
def test_tb_full_mapping_frag_01():
    """
    Test case to ensure that the fragment based full mapping works as expected
    for the first paired end
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    gem_file = resource_path + "tb.Human.GCA_000001405.22_gem.fasta.gem"

    fastq_file_1 = resource_path + "tb.Human.SRR1658573_1.fastq"

    files = [
        gem_file,
        fastq_file_1
    ]

    metadata = {
        'assembly' : 'test',
        'enzyme_name' : 'MboI',
        'windows' : None
    }

    gem_file = files[1]

    print(gem_file)

    tfm1 = tbFullMappingTool()
    tfm1_files, tfm1_meta = tfm1.run(files, [], metadata)

    map_frag = resource_path + "tb.Human.SRR1658573_1_frag.map"
    map_full = resource_path + "tb.Human.SRR1658573_1_full.map"

    assert os.path.isfile(map_frag) is True
    assert os.path.getsize(map_frag) > 0
    assert os.path.isfile(map_full) is True
    assert os.path.getsize(map_full) > 0

@pytest.mark.hic
def test_tb_full_mapping_frag_02():
    """
    Test case to ensure that the fragment based full mapping works as expected
    for the second paired end
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    gem_file = resource_path + "tb.Human.GCA_000001405.22_gem.fasta.gem"

    fastq_file_2 = resource_path + "tb.Human.SRR1658573_2.fastq"

    files = [
        gem_file,
        fastq_file_2
    ]

    metadata = {
        'assembly' : 'test',
        'enzyme_name' : 'MboI',
        'windows' : None
    }

    gem_file = files[1]

    print(gem_file)

    tfm2 = tbFullMappingTool()
    tfm2_files, tfm2_meta = tfm2.run(files, [], metadata)

    map_frag = resource_path + "tb.Human.SRR1658573_2_frag.map"
    map_full = resource_path + "tb.Human.SRR1658573_2_full.map"

    assert os.path.isfile(map_frag) is True
    assert os.path.getsize(map_frag) > 0
    assert os.path.isfile(map_full) is True
    assert os.path.getsize(map_full) > 0

@pytest.mark.hic
def test_tb_full_mapping_iter_01():
    """
    Test case to ensure that the iterative based full mapping works as expected
    for the first paired end
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    gem_file = resource_path + "tb.Human.GCA_000001405.22_gem.fasta.gem"

    fastq_file_1 = resource_path + "tb.Human.SRR1658573_1.fastq"

    files = [
        gem_file,
        fastq_file_1
    ]

    metadata = {
        'assembly' : 'test',
        # 'enzyme_name' : 'MboI',
        'windows' : ((1, 25), (1, 50), (1, 75), (1, 100))
    }

    gem_file = files[1]

    print(gem_file)

    tfm1 = tbFullMappingTool()
    tfm1_files, tfm1_meta = tfm1.run(files, [], metadata)

    map25 = resource_path + "tb.Human.SRR1658573_1_full_1-25.map"
    map50 = resource_path + "tb.Human.SRR1658573_1_full_1-50.map"
    map75 = resource_path + "tb.Human.SRR1658573_1_full_1-75.map"
    map100 = resource_path + "tb.Human.SRR1658573_1_full_1-100.map"

    assert os.path.isfile(map25) is True
    assert os.path.getsize(map25) > 0
    assert os.path.isfile(map50) is True
    assert os.path.getsize(map50) > 0
    assert os.path.isfile(map75) is True
    assert os.path.getsize(map75) > 0
    assert os.path.isfile(map100) is True
    assert os.path.getsize(map100) > 0

@pytest.mark.hic
def test_tb_full_mapping_iter_02():
    """
    Test case to ensure that the iterative based full mapping works as expected
    for the second paired end
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    gem_file = resource_path + "tb.Human.GCA_000001405.22_gem.fasta.gem"

    fastq_file_2 = resource_path + "tb.Human.SRR1658573_2.fastq"

    files = [
        gem_file,
        fastq_file_2
    ]

    metadata = {
        'assembly' : 'test',
        # 'enzyme_name' : 'MboI',
        'windows' : ((1, 25), (1, 50), (1, 75), (1, 100))
    }

    gem_file = files[1]

    print(gem_file)

    tfm2 = tbFullMappingTool()
    tfm2_files, tfm2_meta = tfm2.run(files, [], metadata)

    map25 = resource_path + "tb.Human.SRR1658573_2_full_1-25.map"
    map50 = resource_path + "tb.Human.SRR1658573_2_full_1-50.map"
    map75 = resource_path + "tb.Human.SRR1658573_2_full_1-75.map"
    map100 = resource_path + "tb.Human.SRR1658573_2_full_1-100.map"

    assert os.path.isfile(map25) is True
    assert os.path.getsize(map25) > 0
    assert os.path.isfile(map50) is True
    assert os.path.getsize(map50) > 0
    assert os.path.isfile(map75) is True
    assert os.path.getsize(map75) > 0
    assert os.path.isfile(map100) is True
    assert os.path.getsize(map100) > 0
