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

from process_genome import process_genome
from process_chipseq import process_chipseq

@pytest.mark.chipseq
@pytest.mark.pipeline
def test_chipseq_pipeline():
    """
    Test case to ensure that the ChIP-seq pipeline code works.

    Running the pipeline with the test data from the command line:

    .. code-block:: none
       runcompss                     \
          --lang=python              \
          --library_path=${HOME}/bin \
          --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \
          --log_level=debug          \
          process_chipseq.py         \
             --taxon_id 9606         \
             --genome /<dataset_dir>/Human.GCA_000001405.22.fasta \
             --assembly GRCh38       \
             --file /<dataset_dir>/DRR000150.22.fastq
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    genome_handle = process_genome()
    genome_files, genome_meta = genome_handle.run(
        [resource_path + 'macs2.Human.GCA_000001405.22.fasta'],
        {'assembly' : 'GRCh38'},
        []
    )

    files = [resource_path + 'macs2.Human.GCA_000001405.22.fasta']
    files += genome_files[6:11]
    files += [
        resource_path + 'macs2.Human.DRR000150.22.fastq',
        None
    ]

    metadata = {
        'assembly' : 'GRCh38',
        'expt_name' : 'macs.Human.SRR1658573'
    }

    chipseq_handle = process_chipseq()
    chipseq_files, chipseq_meta = chipseq_handle.run(files, metadata, [])

    print(chipseq_files)

    # Add tests for all files created
    for f_out in chipseq_files:
        print("CHIP-SEQ RESULTS FILE:", f_out)
        assert os.path.isfile(f_out) is True
        assert os.path.getsize(f_out) > 0
