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
import gzip
import pytest # pylint: disable=unused-import

from process_trim_galore import process_trim_galore
from basic_modules.metadata import Metadata

@pytest.mark.trimgalore
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
             --fastq /<dataset_dir>/bsSeeker.Mouse.SRR892982_1.fastq.gz

    """

    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    fastq_file_1 = resource_path + "bsSeeker.Mouse.SRR892982_1.fastq"
    fastq_file_2 = resource_path + "bsSeeker.Mouse.SRR892982_2.fastq"

    with gzip.open(fastq_file_1 + '.gz', 'rb') as fgz_in:
        with open(fastq_file_1, 'w') as f_out:
            f_out.write(fgz_in.read())

    with gzip.open(fastq_file_2 + '.gz', 'rb') as fgz_in:
        with open(fastq_file_2, 'w') as f_out:
            f_out.write(fgz_in.read())

    files = {
        'fastq1': resource_path + 'bsSeeker.Mouse.SRR892982_1.fastq'
    }

    metadata = {
        "fastq1": Metadata(
            "data_wgbs", "fastq", files['fastq1'], None,
        )
    }

    files_out = {
        "fastq1_trimmed": 'tests/data/bsSeeker.Mouse.SRR892982_1_trimmed.fq'
    }

    tg_handle = process_trim_galore()
    #print ("****files_out : *****", files_out )
    tg_files, tg_meta = tg_handle.run(files, metadata, files_out)
    print ("****files : *****", files )
    print ("****metadata : *****", metadata )
    print ("****files_out : *****", files_out )
    print ("****tg_files : *****", tg_files )

    # Checks that the returned files matches the expected set of results
    assert len(tg_files) == 1

    # Add tests for all files created
    for f_out in tg_files:
        print("TRIM GALORE RESULTS FILE:", f_out)
        assert tg_files[f_out] == files_out[f_out]
        assert os.path.isfile(tg_files[f_out]) is True
        assert os.path.getsize(tg_files[f_out]) > 0

def test_trim_galore_pipeline_02():
    """
    Test case to ensure that the trimgalore pipeline code works 
    for paired end data.

    Running the pipeline with the test data from the command line:

    .. code-block:: none

       runcompss                                                         \\
          --lang=python                                                  \\
          --library_path=${HOME}/bin                                     \\
          --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \\
          --log_level=debug                                              \\
          process_trim_galore.py                                         \\
             --taxon_id 9606                                             \\
             --fastq1 /<dataset_dir>/bsSeeker.Mouse.SRR892982_1.fastq.gz \\
             --fastq2 /<dataset_dir>/bsSeeker.Mouse.SRR892982_2.fastq.gz

    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    files = {
        'fastq1': resource_path + 'bsSeeker.Mouse.SRR892982_1.fastq',
        'fastq2': resource_path + 'bsSeeker.Mouse.SRR892982_2.fastq'
    }

    metadata = {
        "fastq1": Metadata(
            "data_wgbs", "fastq", files['fastq1'], None,
        ),

        "fastq2": Metadata(
            "data_wgbs", "fastq", files['fastq2'], None,
        )
    }

    files_out = {
        "fastq1_trimmed": 'tests/data/bsSeeker.Mouse.SRR892982_1.trimmed.fq',
        "fastq2_trimmed": 'tests/data/bsSeeker.Mouse.SRR892982_2.trimmed.fq'
    }

    tg_handle = process_trim_galore()
    tg_files, tg_meta = tg_handle.run(files, metadata, files_out)

    # Checks that the returned files matches the expected set of results
    assert len(tg_files) == 2

    # Add tests for all files created
    for f_out in tg_files:
        print("TRIM GALORE RESULTS FILE:", f_out)
        assert tg_files[f_out] == files_out[f_out]
        assert os.path.isfile(tg_files[f_out]) is True
        assert os.path.getsize(tg_files[f_out]) > 0
