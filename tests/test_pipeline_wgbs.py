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

from process_wgbs import process_wgbs
from basic_modules.metadata import Metadata

@pytest.mark.wgbs
@pytest.mark.pipeline
def test_wgbs_pipeline():
    """
    Test case to ensure that the RNA-seq pipeline code works.

    Running the pipeline with the test data from the command line:

    .. code-block:: none

       runcompss                                                         \\
          --lang=python                                                  \\
          --library_path=${HOME}/bin                                     \\
          --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \\
          --log_level=debug                                              \\
          process_wgbs.py                                                \\
             --taxon_id 10090                                            \\
             --genome /<dataset_dir>/Mouse.GRCm38.fasta                  \\
             --assembly GRCm38                                           \\
             --fastq1 /<dataset_dir>/expt_1.fastq                        \\
             --fastq2 /<dataset_dir>/expt_2.fastq                        \\
             --aligner bowtie2                                           \\
             --aligner_path ${HOME}/lib/bowtie2-2.3.2                    \\
             --bss_path ${HOME}/lib/BSseeker2

    """
    home = os.path.expanduser('~')
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    genomefa_file = resource_path + "bsSeeker.Mouse.GRCm38.fasta"
    fastq1_file = resource_path + "bsSeeker.Mouse.GRCm38_1.fastq"
    fastq2_file = resource_path + "bsSeeker.Mouse.GRCm38_2.fastq"

    files = {
        "genome" : genomefa_file,
        "fastq1" : fastq1_file,
        "fastq2" : fastq2_file
    }

    metadata = {
        "genome": Metadata(
            "Assembly", "fasta", files['genome'], None,
            {'assembly' : 'GRCm38'}
        ),
        "fastq1": Metadata(
            "data_wgbs_seq", "fastq", files['fastq1'], None,
            {"assembly" : "GRCm38"}
        ),
        "fastq2": Metadata(
            "data_wgbs_seq", "fastq", files['fastq2'], None,
            {"assembly" : "GRCm38"}
        ),
    }

    files_out = {
        "index" : resource_path + "wgbs.Mouse.GRCm38.fasta.bt2.tar.gz",
        "fastq1_filtered" : resource_path + 'bsSeeker.Mouse.GRCm38_1_filtered.fastq',
        "fastq2_filtered" : resource_path + 'bsSeeker.Mouse.GRCm38_2_filtered.fastq',
        "bam" : resource_path + "bsSeeker.Mouse.GRCm38_1_filtered.bam",
        "bai" : resource_path + "bsSeeker.Mouse.GRCm38_1_filtered.bai",
        "wig_file" : resource_path + "bsSeeker.Mouse.GRCm38_1.wig",
        "cgmap_file" : resource_path + "bsSeeker.Mouse.GRCm38_1.cgmap",
        "atcgmap_file" : resource_path + "bsSeeker.Mouse.GRCm38_1.atcgmap"
    }

    print("WGBS TEST FILES:", files)
    rs_handle = process_wgbs(
        configuration={
            "bss_path" : home + "/lib/BSseeker2",
            "aligner" : "bowtie2",
            "aligner_path" : home + "/lib/bowtie2-2.3.2"
        }
    )
    rs_files, rs_meta = rs_handle.run(files, metadata, files_out)

    print("WGBS RESULTS FILES:", len(rs_files), rs_files)
    # Checks that the returned files matches the expected set of results
    assert len(rs_files) == 8

    # Add tests for all files created
    for f_out in rs_files:
        print("WGBS RESULTS FILE:", f_out)
        assert rs_files[f_out] == files_out[f_out]
        assert os.path.isfile(rs_files[f_out]) is True
        assert os.path.getsize(rs_files[f_out]) > 0
