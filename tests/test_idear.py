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

from basic_modules.metadata import Metadata

from tool.idear import idearTool

@pytest.mark.idamidseq
def test_idear():
    """
    Function to test forging BSgenomes
    """

    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "bsgenome": resource_path + "idear.Human.GCA_000001405.22.22.bsgenome.tar.gz",
        "bam_1": resource_path + "idear.Human.SRR3714775.bam",
        "bam_2": resource_path + "idear.Human.SRR3714776.bam",
        "bg_bam_1": resource_path + "idear.Human.SRR3714777.bam",
        "bg_bam_2": resource_path + "idear.Human.SRR3714778.bam",
    }

    output_files = {
        "bigwig": resource_path + "idear.Human.Nup98-GFP.bw"
    }

    metadata = {
        "bsgenome": Metadata(
            "data_damid_seq", "bsgenome", [], None,
            {'assembly' : 'test'}, 9606),
        "bam_1": Metadata(
            "data_damid_seq", "bam", [], None,
            {'assembly' : 'test'}, 9606),
        "bam_2": Metadata(
            "data_damid_seq", "bam", [], None,
            {'assembly' : 'test'}, 9606),
        "bg_bam_1": Metadata(
            "data_damid_seq", "bam", [], None,
            {'assembly' : 'test'}, 9606),
        "bg_bam_2": Metadata(
            "data_damid_seq", "bam", [], None,
            {'assembly' : 'test'}, 9606),
    }

    config = {
        "idear_common_name": "Human",
        "idear_sample_param": "Nup98",
        "idear_background_param": "GFP"
    }

    idear_handle = idearTool(config)
    idear_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(resource_path + "idear.Human.Nup98-GFP.bw") is True
    assert os.path.getsize(resource_path + "idear.Human.Nup98-GFP.bw") > 0
