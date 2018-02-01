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

from process_bs_seeker_index import process_bs_seeker_index
from basic_modules.metadata import Metadata

@pytest.mark.wgbs
@pytest.mark.pipeline
def test_wgbs_pipeline_index():
    """
    Test case to ensure that the RNA-seq pipeline code works.

    Running the pipeline with the test data from the command line:
    """
    home = os.path.expanduser('~')
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    genomefa_file = resource_path + "bsSeeker.Mouse.GRCm38.fasta"

    files = {
        "genome" : genomefa_file
    }

    metadata = {
        "genome": Metadata(
            "Assembly", "fasta", files['genome'], None,
            {'assembly' : 'GRCm38'}
        )
    }

    files_out = {
        "index" : resource_path + "wgbs.Mouse.GRCm38.fasta.bt2.tar.gz"
    }

    print("WGBS TEST FILES:", files)
    rs_handle = process_bs_seeker_index(
        configuration={
            "bss_path" : home + "/lib/BSseeker2",
            "aligner" : "bowtie2",
            "aligner_path" : home + "/lib/bowtie2-2.3.2"
        }
    )
    rs_files, rs_meta = rs_handle.run(files, metadata, files_out)

    assert len(rs_files) == 1

    # Add tests for all files created
    for f_out in rs_files:
        print("WGBS RESULTS FILE:", f_out)
        assert rs_files[f_out] == files_out[f_out]
        assert os.path.isfile(rs_files[f_out]) is True
        assert os.path.getsize(rs_files[f_out]) > 0
