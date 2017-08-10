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

@pytest.mark.wgbs
@pytest.mark.pipeline
def test_wgbs_pipeline():
    """
    Test case to ensure that the RNA-seq pipeline code works.

    Running the pipeline with the test data from the command line:

    .. code-block:: none
       runcompss                                  \
          --lang=python                           \
          --library_path=${HOME}/bin              \
          --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \
          --log_level=debug                       \
          process_wgbs.py                         \
             --taxon_id 10090                     \
             --genome /<dataset_dir>/Mouse.GRCm38.fasta \
             --assembly GRCm38                    \
             --fastq1 /<dataset_dir>/expt_1.fastq \
             --fastq2 /<dataset_dir>/expt_2.fastq \
             --aligner bowtie2                    \
             --aligner_path ${HOME}/lib/bowtie2-2.3.2 \
             --bss_path ${HOME}/lib/BSseeker2
    """
    home = os.path.expanduser('~')
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    genomefa_file = resource_path + "bsSeeker.Mouse.GRCm38.fasta"
    fastq1_file = resource_path + "bsSeeker.Mouse.GRCm38_1.fastq"
    fastq2_file = resource_path + "bsSeeker.Mouse.GRCm38_2.fastq"

    rs_handle = process_rnaseq()
    rs_files, rs_meta = rs_handle.run(
        [
            genomefa_file,
            fastq1_file,
            fastq2_file
        ],
        {
            'assembly' : 'GRCh38',
            'aligner' : 'bowtie2',
            'aligner_path' : home + '/lib/bowtie2-2.3.2',
            'bss_path' : home + '/lib/BSseeker2'
        },
        []
    )

    print(rs_files)

    # Add tests for all files created
    for f_out in rs_files:
        print("WGBS RESULTS FILE:", f_out)
        assert os.path.isfile(f_out) is True
        assert os.path.getsize(f_out) > 0
