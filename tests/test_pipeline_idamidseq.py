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

from process_damidseq import process_damidseq
from basic_modules.metadata import Metadata

@pytest.mark.idamidseq
@pytest.mark.pipeline
def test_idamidseq_pipeline():
    """
    Test case to ensure that the ChIP-seq pipeline code works.

    Running the pipeline with the test data from the command line:

    .. code-block:: none

       runcompss                                                         \\
          --lang=python                                                  \\
          --library_path=${HOME}/bin                                     \\
          --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \\
          --log_level=debug                                              \\
          process_chipseq.py                                             \\
             --taxon_id 9606                                             \\
             --genome /<dataset_dir>/Human.GCA_000001405.22.fasta        \\
             --assembly GRCh38                                           \\
             --file /<dataset_dir>/DRR000150.22.fastq

    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    files = {
        'genome': resource_path + 'idear.Human.GCA_000001405.22.fasta',
        'index': resource_path + 'idear.Human.GCA_000001405.22.fasta.bwa.tar.gz',
        'fastq_1': resource_path + 'idear.Human.SRR3714775.fastq',
        'fastq_2': resource_path + 'idear.Human.SRR3714776.fastq',
        'bg_fastq_1': resource_path + 'idear.Human.SRR3714777.fastq',
        'bg_fastq_2': resource_path + 'idear.Human.SRR3714778.fastq',
    }

    metadata = {
        "genome": Metadata(
            "Assembly", "fasta", files['genome'], None,
            {'assembly' : 'GCA_000001405.22'}),
        "index": Metadata(
            "Index", "bwa_index", files['index'], files['genome'],
            {'assembly': 'GCA_000001405.22', "tool": "bwa_indexer"}),
        "fastq_1": Metadata(
            "data_idamid_seq", "fastq", files['fastq_1'], None,
            {'assembly' : 'GCA_000001405.22'}
        ),
        "fastq_2": Metadata(
            "data_idamid_seq", "fastq", files['fastq_2'], None,
            {'assembly' : 'GCA_000001405.22'}
        ),
        "bg_fastq_1": Metadata(
            "data_idamid_seq", "fastq", files['bg_fastq_1'], None,
            {'assembly' : 'GCA_000001405.22'}
        ),
        "bg_fastq_2": Metadata(
            "data_idamid_seq", "fastq", files['bg_fastq_2'], None,
            {'assembly' : 'GCA_000001405.22'}
        ),
    }

    config_param = {
        "idear_title": "Full genome sequences for Homo sapiens (GRCh38)",
        "idear_description": "Full genome sequences for Homo sapiens (GRCh38)",
        "idear_common_name": "Human",
        "idear_organism": "Homo sapiens",
        "idear_provider": "ENA",
        "idear_release_date": "2013",
        "idear_sample_param": "Nup98",
        "idear_background_param": "GFP",
    }

    files_out = {
        "bam_1": files['fastq_1'].replace(".fastq", ".bam"),
        "bam_2": files['fastq_2'].replace(".fastq", ".bam"),
        "bg_bam_1": files['bg_fastq_1'].replace(".fastq", ".bam"),
        "bg_bam_2": files['bg_fastq_2'].replace(".fastq", ".bam"),
        "bam_1_filtered": files['fastq_1'].replace(".fastq", ".filtered.bam"),
        "bam_2_filtered": files['fastq_2'].replace(".fastq", ".filtered.bam"),
        "bg_bam_1_filtered": files['bg_fastq_1'].replace(".fastq", ".filtered.bam"),
        "bg_bam_2_filtered": files['bg_fastq_2'].replace(".fastq", ".filtered.bam"),
        "bsgenome": resource_path + "idear.Human.GCA_000001405.22.22.bsgenome.tar.gz",
        "chrom_size": resource_path + "chrom.size",
        "genome_2bit": resource_path + "idear.Human.GCA_000001405.22.2bit",
        "seed_file": resource_path + "idear.Human.GCA_000001405.22.seed",
        "bigwig": resource_path + "idear.Human.Nup98-GFP.bw"
    }

    damidseq_handle = process_damidseq(config_param)
    damidseq_files, damidseq_meta = damidseq_handle.run(files, metadata, files_out)

    print(damidseq_files)

    # Add tests for all files created
    for f_out in damidseq_files:
        print("iDamID-SEQ RESULTS FILE:", f_out)
        #assert damidseq_files[f_out] == files_out[f_out]
        assert os.path.isfile(damidseq_files[f_out]) is True
        assert os.path.getsize(damidseq_files[f_out]) > 0
