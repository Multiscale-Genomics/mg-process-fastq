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

import os
import tarfile
import pytest # pylint: disable=unused-import

from tool import bs_seeker_aligner

@pytest.mark.wgbs
def test_bs_seeker_aligner():
    """
    Test to ensure bs-Seeker aligner works
    """

    home = os.path.expanduser('~')
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genomefa_file = resource_path + "bsSeeker.Mouse.GRCm38.fasta"
    genomeidx_file = resource_path + "bsSeeker.Mouse.GRCm38.fasta_bowtie2.tar.gz"
    fastq_gz = resource_path + "bsSeeker.Mouse.GRCm38_1.fastq.filtered.fastq.tar.gz"
    fastq1_file = "bsSeeker.Mouse.GRCm38_1.fastq.filtered.fastq"
    fastq2_file = "bsSeeker.Mouse.GRCm38_2.fastq.filtered.fastq"

    tar = tarfile.open(fastq_gz, "w:gz")
    tar.add(resource_path + fastq1_file, arcname='tmp/' + fastq1_file)
    tar.add(resource_path + fastq2_file, arcname='tmp/' + fastq2_file)
    tar.close()

    bsa = bs_seeker_aligner.bssAlignerTool()
    bs_files, bs_meta = bsa.run(
        [
            genomefa_file,
            genomeidx_file,
            fastq_gz
        ],
        [],
        {
            "aligner" : "bowtie2",
            "aligner_path" : home + "/lib/bowtie2-2.3.2",
            "bss_path" : home + "/lib/BSseeker2",
            "expt_name" : "bsSeeker.Mouse.GRCm38",
            "fastq_list" : [[fastq1_file, fastq2_file]]
        }
    )

    assert os.path.isfile(fastq_gz) is True
    assert os.path.getsize(fastq_gz) > 0
    assert os.path.isfile(bs_files[0]) is True
    assert os.path.getsize(bs_files[0]) > 0
    assert os.path.isfile(bs_files[1]) is True
    assert os.path.getsize(bs_files[1]) > 0
