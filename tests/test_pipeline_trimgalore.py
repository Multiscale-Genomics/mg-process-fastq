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

from process_trim_galore import process_trim_galore
from basic_modules.metadata import Metadata

@pytest.mark.rnaseq
@pytest.mark.pipeline
def test_trim_galore_pipeline():
    """
    Test case to ensure that the trimgalore pipeline code works.

    Running the pipeline with the test data from the command line:

    .. code-block:: none

       runcompss                                                         \\
          --lang=python                                                  \\
          --library_path=${HOME}/bin                                     \\
          --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \\
          --log_level=debug                                              \\
          process_trim_galore.py                                         \\
             --taxon_id 9606                                             \\
             --file /<dataset_dir>/bsSeeker.Mouse.SRR892982_1.fastq.gz

    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    files = {
        'fastq': resource_path + 'bsSeeker.Mouse.SRR892982_1.fastq.gz'
    }

    metadata = {
        "fastq1": Metadata(
            "data_wgbs", "fastq", files['fastq1'], None,
            {'assembly' : 'SRR892982_1'}
        )
    }

    files_out = {
        "fastq.trimmed": 'tests/data/bsSeeker.Mouse.SRR892982_1.trimmed.fq'
    }

    rs_handle = process_rnaseq()
    rs_files, rs_meta = rs_handle.run(files, metadata, files_out)

    # Checks that the returned files matches the expected set of results
    assert len(rs_files) == 4

    # Add tests for all files created
    for f_out in rs_files:
        print("TRIM GALORE RESULTS FILE:", f_out)
        assert rs_files[f_out] == files_out[f_out]
        assert os.path.isfile(rs_files[f_out]) is True
        assert os.path.getsize(rs_files[f_out]) > 0