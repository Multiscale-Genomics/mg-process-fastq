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
import pytest

from basic_modules.metadata import Metadata
from process_sleuth import process_sleuth


@pytest.mark.sleuth
@pytest.mark.pipeline
def test_sleuth_pipeline():
    """
    Test case to ensure that the Sleuth pipeline code works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data", "sleuth")

    # The files need moving to a file structure like;
    #     "results/ERR030903/sleuth.Human.ERR030903.tar.gz"
    # Then that needs to be tar'ed up and passed to a @task for decompression
    # and analysis

    input_files = {
        "kallisto": [
            os.path.join(resource_path, "sleuth.Human.ERR030856.tar.gz"),
            os.path.join(resource_path, "sleuth.Human.ERR030857.tar.gz"),
            os.path.join(resource_path, "sleuth.Human.ERR030858.tar.gz"),
            os.path.join(resource_path, "sleuth.Human.ERR030872.tar.gz"),
            os.path.join(resource_path, "sleuth.Human.ERR030903.tar.gz")
        ],
    }

    output_files = {
        "sleuth_object": os.path.join(resource_path, "sleuth.Rbin"),
        "sleuth_sig_genes_table": os.path.join(resource_path, "sleuth_sig_genes.tsv"),
        "sleuth_image_tar": os.path.join(resource_path, "sleuth_images.tar.gz")
    }

    metadata = {
        "kallisto": [
            Metadata(
                "data_rna_seq", "TAR",
                os.path.join(resource_path, "sleuth.Human.ERR030856.tar.gz"), None,
                {'assembly': 'test', "dataset": "ERR030856", "condition:tissue": "mixture"}, 9606),
            Metadata(
                "data_rna_seq", "TAR",
                os.path.join(resource_path, "sleuth.Human.ERR030857.tar.gz"), None,
                {'assembly': 'test', "dataset": "ERR030857", "condition:tissue": "mixture"}, 9606),
            Metadata(
                "data_rna_seq", "TAR",
                os.path.join(resource_path, "sleuth.Human.ERR030858.tar.gz"), None,
                {'assembly': 'test', "dataset": "ERR030858", "condition:tissue": "mixture"}, 9606),
            Metadata(
                "data_rna_seq", "TAR",
                os.path.join(resource_path, "sleuth.Human.ERR030872.tar.gz"), None,
                {'assembly': 'test', "dataset": "ERR030872", "condition:tissue": "thyroid"}, 9606),
            Metadata(
                "data_rna_seq", "TAR",
                os.path.join(resource_path, "sleuth.Human.ERR030903.tar.gz"), None,
                {'assembly': 'test', "dataset": "ERR030903", "condition:tissue": "thyroid"}, 9606),
        ]
    }

    sleuth_config = {
        "sleuth_sig_level": 1.0,
        "sleuth_tag": "test"
    }

    sleuth_handle = process_sleuth(sleuth_config)
    sleuth_files, sleuth_meta = sleuth_handle.run(  # pylint: disable=unused-variable
        input_files, metadata, output_files)

    # Checks that the returned files matches the expected set of results
    assert len(sleuth_files) == 3

    # Add tests for all files created
    for f_out in sleuth_files:
        print("Sleuth RESULTS FILE:", f_out)
        assert sleuth_files[f_out] == output_files[f_out]
        assert os.path.isfile(sleuth_files[f_out]) is True
        assert os.path.getsize(sleuth_files[f_out]) > 0
