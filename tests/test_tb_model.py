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
import pytest  # pylint: disable=unused-import

from tool.tb_model import tbModelTool


@pytest.mark.hic
def test_tb_model():
    """
    Test case to build a model from a matrix file
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    hic_contacts_matrix_norm = resource_path + "tb_nrm_chr21_100kb.abc"

    files = [hic_contacts_matrix_norm]

    metadata = {
        'resolution': "100000",
        'optimize_only': False,
        'gen_pos_chrom_name': "chr21",
        'gen_pos_begin': '15400000',
        'gen_pos_end': '20000000',
        'max_dist': '10000',
        'upper_bound': '0',
        'lower_bound': '-0.6',
        'cutoff': '2',
        'num_mod_comp': '4',
        'num_mod_keep': '4',
        "species": "Homo sapiens",
        "assembly": "Unknown",
        'ncpus': 4
    }

    tm = tbModelTool()
    tm_files, tm_meta = tm.run(files, [], metadata)

    assert len(tm_files) == 2
