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

from process_hic import process_hic

@pytest.mark.hic
@pytest.mark.pipeline
def test_tb_pipeline():
    """
    Test case to ensure that the Hi-C pipeline code works.

    Running the pipeline with the test data from the command line:

    .. code-block:: none

       runcompss                                                                 \\
          --lang=python                                                          \\
          --library_path=/home/compss/bin                                        \\
          --pythonpath=/<pyenv_virtenv_dir>//lib/python2.7/site-packages/        \\
          --log_level=debug                                                      \\
          process_hic.py                                                         \\
             --taxon_id 9606                                                     \\
             --genome /<dataset_dir>/tb.Human.GCA_000001405.22_gem.fasta         \\
             --assembly GRCh38                                                   \\
             --file1 /<dataset_dir>/tb.Human.SRR1658573_1.fastq                  \\
             --file2 /<dataset_dir>/tb.Human.SRR1658573_2.fastq                  \\
             --genome_gem /<dataset_dir>/tb.Human.GCA_000001405.22_gem.fasta.gem \\
             --taxon_id 9606                                                     \\
             --enzyme_name MboI                                                  \\
             --resolutions 10000,100000                                          \\
             --windows1 1,100                                                    \\
             --windows2 1,100                                                    \\
             --normalized 1                                                      \\
             --tag tb.Human.SRR1658573                                           \\
             --window_type frag                                                  \\

    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    files = [
        resource_path + 'tb.Human.GCA_000001405.22_gem.fasta',
        resource_path + 'tb.Human.GCA_000001405.22_gem.fasta.gem',
        resource_path + 'tb.Human.SRR1658573_1.fastq',
        resource_path + 'tb.Human.SRR1658573_2.fastq'
    ]

    metadata = {
        'assembly' : 'GRCh38',
        'expt_name' : 'tb.Human.SRR1658573',
        'enzyme_name' : 'MboI',
        'windows1' : ((1, '100')),
        'windows2' : ((1, '100')),
        'window_type' :  'frag',
        'resolutions' : [10000, 100000],
        'normalized' : False,
        'hdf5' : True,
    }

    hic_handle = process_hic()
    hic_files, hic_meta = hic_handle.run(files, metadata, [])

    print(hic_files)

    # Add tests for all files created
    assert os.path.isfile(hic_files[0]) is True
    assert os.path.getsize(hic_files[0]) > 0
