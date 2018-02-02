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

from process_idear import process_idear
from basic_modules.metadata import Metadata

@pytest.mark.idamidseq
@pytest.mark.pipeline
def test_idear_pipeline():
    """
    Test case to ensure that the iDEAR pipeline code works.

    Running the pipeline with the test data from the command line:
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    files = {
        'bsgenome': resource_path + "idear.Human.GCA_000001405.22.22.bsgenome.tar.gz",
        'bam_1': resource_path + 'idear.Human.SRR3714775.bam',
        'bam_2': resource_path + 'idear.Human.SRR3714776.bam',
        'bg_bam_1': resource_path + 'idear.Human.SRR3714777.bam',
        'bg_bam_2': resource_path + 'idear.Human.SRR3714778.bam',
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

    config_param = {
        "idear_title": "Full genome sequences for Homo sapiens (GRCh38)",
        "idear_description": "Full genome sequences for Homo sapiens (GRCh38)",
        "idear_common_name": "Human",
        "idear_organism": "Homo sapiens",
        "idear_provider": "ENA",
        "idear_release_date": "2013",
        "idear_sample_param": "Nup98",
        "idear_background_param": "GFP",
    }

    damidseq_handle = process_idear(config_param)
    damidseq_files, damidseq_meta = damidseq_handle.run(files, metadata, output_files)

    print(damidseq_files)

    # Add tests for all files created
    for f_out in damidseq_files:
        assert os.path.isfile(damidseq_files[f_out]) is True
        assert os.path.getsize(damidseq_files[f_out]) > 0
