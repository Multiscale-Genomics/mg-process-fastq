"""
.. Copyright 2017 EMBL-European Bioinformatics Institute

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

from process_rnaseq import process_rnaseq

@pytest.mark.rnaseq
@pytest.mark.pipeline
def test_rnaseq_pipeline():
    """
    Test case to ensure that the RNA-seq pipeline code works.

    Running the pipeline with the test data from the command line:

    .. code-block:: none

       runcompss                                                         \\
          --lang=python                                                  \\
          --library_path=${HOME}/bin                                     \\
          --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \\
          --log_level=debug                                              \\
          process_rnaseq.py                                              \\
             --taxon_id 9606                                             \\
             --genome /<dataset_dir>/Human.GRCh38.fasta                  \\
             --assembly GRCh38                                           \\
             --file /<dataset_dir>/ERR030872_1.fastq                     \\
             --file2 /<dataset_dir>/ERR030872_2.fastq

    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    rs_handle = process_rnaseq()
    rs_files, rs_meta = rs_handle.run(
        [
            resource_path + 'kallisto.Human.GRCh38.fasta',
            resource_path + 'kallisto.Human.ERR030872_1.fastq',
            resource_path + 'kallisto.Human.ERR030872_2.fastq'
        ],
        {'assembly' : 'GRCh38'},
        []
    )

    print(rs_files)

    # Add tests for all files created
    for f_out in rs_files:
        print("RNA SEQ RESULTS FILE:", f_out)
        assert os.path.isfile(f_out) is True
        assert os.path.getsize(f_out) > 0
