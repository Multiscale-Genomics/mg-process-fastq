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

from tool.forge_bsgenome import bsgenomeTool

@pytest.mark.idamidseq
def test_bsgenome():
    """
    Function to test forging BSgenomes
    """

    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "genome": resource_path + "idear.Human.GCA_000001405.22.fasta"
    }

    output_files = {
        "bsgenome": resource_path + "idear.Human.GCA_000001405.22.22.bsgenome.tar.gz",
        "chrom_size": resource_path + "chrom.size",
        "genome_2bit": resource_path + "idear.Human.GCA_000001405.22.2bit",
        "seed_file": resource_path + "idear.Human.GCA_000001405.22.seed"
    }

    metadata = {
        "genome": Metadata(
            "data_damid_seq", "fastq", [], None,
            {'assembly' : 'test'}),
    }

    config = {
        "idear_title": "Full genome sequences for Homo sapiens (GRCh38)",
        "idear_description": "Full genome sequences for Homo sapiens (GRCh38)",
        "idear_common_name": "Human",
        "idear_organism": "Homo sapiens",
        "idear_provider": "ENA",
        "idear_release_date": "2013",
        #"idear_circ_chrom": "",
    }

    bsg_handle = bsgenomeTool(config)
    bsg_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(resource_path + "idear.Human.GCA_000001405.22.22.bsgenome.tar.gz") is True
    assert os.path.getsize(resource_path + "idear.Human.GCA_000001405.22.22.bsgenome.tar.gz") > 0
    assert os.path.isfile(resource_path + "chrom.size") is True
    assert os.path.getsize(resource_path + "chrom.size") > 0
    assert os.path.isfile(resource_path + "idear.Human.GCA_000001405.22.2bit") is True
    assert os.path.getsize(resource_path + "idear.Human.GCA_000001405.22.2bit") > 0
    assert os.path.isfile(resource_path + "idear.Human.GCA_000001405.22.seed") is True
    assert os.path.getsize(resource_path + "idear.Human.GCA_000001405.22.seed") > 0
