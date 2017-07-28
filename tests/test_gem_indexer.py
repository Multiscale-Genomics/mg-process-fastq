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

from tool.gem_indexer import gemIndexerTool

def test_gem_indexer():
    """
    Test case to ensure that the GEM indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genome_fa = resource_path + "tb.Human.GCA_000001405.22.fasta"
    genome_gem_fa = resource_path + "tb.Human.GCA_000001405.22_gem.fasta"

    with gzip.open(genome_fa + '.gz', 'rb') as fgz_in:
        with open(genome_fa, 'w') as f_out:
            f_out.write(fgz_in.read())


    genome_gem = resource_path + "tb.Human.GCA_000001405.22_gem.fasta.gem"

    gem_it = gemIndexerTool()
    gem_it.run([genome_fa], {'assembly' : 'test'})

    assert os.path.isfile(genome_gem_fa) is True
    assert os.path.getsize(genome_gem_fa) > 0
    assert os.path.isfile(genome_gem) is True
    assert os.path.getsize(genome_gem) > 0
