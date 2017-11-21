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

from basic_modules.metadata import Metadata

from process_genome import process_genome

@pytest.mark.genome
@pytest.mark.chipseq
@pytest.mark.pipeline
def test_genome_pipeline():
    """
    Test case to ensure that the Genome indexing pipeline code works.

    Running the pipeline with the test data from the command line:

    .. code-block:: none

       runcompss                                                         \\
          --lang=python                                                  \\
          --library_path=${HOME}/bin                                     \\
          --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \\
          --log_level=debug                                              \\
          process_genome.py                                              \\
             --taxon_id 9606                                             \\
             --genome /<dataset_dir>/Human.GCA_000001405.22.fasta        \\
             --assembly GRCh38                                           \\
             --file /<dataset_dir>/DRR000150.22.fastq
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    files = {
        'genome': resource_path + 'macs2.Human.GCA_000001405.22.fasta',
    }

    metadata = {
        "genome": Metadata(
            "Assembly", "fasta", files['genome'], None,
            {'assembly' : 'GCA_000001405.22'}),
    }

    files_out = {
        "bwa_index": resource_path + 'macs2.Human.GCA_000001405.22.fasta.bwa.tar.gz',
        "bwt_index": resource_path + 'macs2.Human.GCA_000001405.22.fasta.bt2.tar.gz',
        "gem_index": resource_path + 'macs2.Human.GCA_000001405.22.gem.fasta.gem.gz',
        "genome_gem": resource_path + 'macs2.Human.GCA_000001405.22.gem.fasta'
    }

    genome_handle = process_genome()
    genome_files, genome_meta = genome_handle.run(files, metadata, files_out)

    print(genome_files)

    # Add tests for all files created
    for f_out in genome_files:
        print("GENOME RESULTS FILE:", f_out)
        assert genome_files[f_out] == files_out[f_out]
        assert os.path.isfile(genome_files[f_out]) is True
        assert os.path.getsize(genome_files[f_out]) > 0
