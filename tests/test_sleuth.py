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
import tarfile
import pytest  # pylint: disable=unused-import

from basic_modules.metadata import Metadata

from tool.sleuth import sleuthTool


@pytest.mark.sleuth
def test_sleuth():
    """
    Function to test the sleuth analysis tool
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

    input_files = {
        "kallisto_tar": resource_path + "results.tar.gz",
    }

    output_files = {
        "sleuth_object": resource_path + "sleuth.Rbin"
    }

    metadata = {
        "kallisto_tar": Metadata(
            "data_rna_seq", "TAR", [
                "results/ERR030856/sleuth.Human.ERR030856.tar.gz",
                "results/ERR030856/sleuth.Human.ERR030857.tar.gz",
                "results/ERR030856/sleuth.Human.ERR030858.tar.gz",
                "results/ERR030856/sleuth.Human.ERR030872.tar.gz",
                "results/ERR030856/sleuth.Human.ERR030903.tar.gz",
            ], None,
            {'assembly': 'test'}, 9606),
    }

    config = {
        "kallisto_tar_config": {
            "ERR030856": "mixture",
            "ERR030857": "mixture",
            "ERR030858": "mixture",
            "ERR030872": "thyroid",
            "ERR030903": "thyroid"
        }
    }

    sleuth_handle = sleuthTool(config)
    sleuth_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(resource_path + "sleuth.Rbin") is True
    assert os.path.getsize(resource_path + "sleuth.Rbin") > 0
