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
import tarfile
import pytest

from basic_modules.metadata import Metadata
from process_sleuth import process_sleuth


@pytest.mark.sleuth
@pytest.mark.pipeline
def test_sleuth_pipeline():
    """
    Test case to ensure that the Sleuth pipeline code works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/sleuth/")

    tar = tarfile.open(resource_path + "results.tar.gz", "w")
    tar.add(
        resource_path + "sleuth.Human.ERR030856.tar.gz",
        arcname="results/ERR030856/sleuth.Human.ERR030856.tar.gz")
    tar.add(
        resource_path + "sleuth.Human.ERR030857.tar.gz",
        arcname="results/ERR030857/sleuth.Human.ERR030857.tar.gz")
    tar.add(
        resource_path + "sleuth.Human.ERR030858.tar.gz",
        arcname="results/ERR030858/sleuth.Human.ERR030858.tar.gz")
    tar.add(
        resource_path + "sleuth.Human.ERR030872.tar.gz",
        arcname="results/ERR030872/sleuth.Human.ERR030872.tar.gz")
    tar.add(
        resource_path + "sleuth.Human.ERR030903.tar.gz",
        arcname="results/ERR030903/sleuth.Human.ERR030903.tar.gz")
    tar.close()

    files = {
        'kallisto_tar': resource_path + 'results.tar.gz'
    }

    metadata = {
        "kallisto_tar": Metadata(
            "data_rna_seq", "TAR", [
                "results/ERR030856/sleuth.Human.ERR030856.tar.gz",
                "results/ERR030857/sleuth.Human.ERR030857.tar.gz",
                "results/ERR030858/sleuth.Human.ERR030858.tar.gz",
                "results/ERR030872/sleuth.Human.ERR030872.tar.gz",
                "results/ERR030903/sleuth.Human.ERR030903.tar.gz"
            ], None,
            {'assembly': 'GRCh38'}),
    }

    files_out = {
        "sleuth_object": resource_path + "sleuth.Rbin",
        "sleuth_sig_genes_table": resource_path + "sleuth_sig_genes.tsv",
        "sleuth_image_tar": resource_path + "sleuth_images.tar.gz"
    }

    sleuth_config = {
        "kallisto_tar_config": {
            "ERR030856": {"tissue": "mixture"},
            "ERR030857": {"tissue": "mixture"},
            "ERR030858": {"tissue": "mixture"},
            "ERR030872": {"tissue": "thyroid"},
            "ERR030903": {"tissue": "thyroid"}
        },
        "sleuth_sig_level": 1.0,
        "sleuth_tag": "test"
    }

    sleuth_handle = process_sleuth(sleuth_config)
    sleuth_files, sleuth_meta = sleuth_handle.run(  # pylint: disable=unused-variable
        files, metadata, files_out)

    # Checks that the returned files matches the expected set of results
    assert len(sleuth_files) == 3

    # Add tests for all files created
    for f_out in sleuth_files:
        print("Sleuth RESULTS FILE:", f_out)
        assert sleuth_files[f_out] == files_out[f_out]
        assert os.path.isfile(sleuth_files[f_out]) is True
        assert os.path.getsize(sleuth_files[f_out]) > 0
