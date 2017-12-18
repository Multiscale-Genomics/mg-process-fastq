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
import pytest

from basic_modules.metadata import Metadata

from tool import inps

@pytest.mark.py3
@pytest.mark.mnaseseq
def test_inps():
    """
    Function to test INPS works.
    """

    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genome_fa = resource_path + "inps.Mouse.GRCm38.fasta"
    fastq_file = resource_path + "inps.Mouse.DRR000386.fastq"
    bam_file = resource_path + "inps.Mouse.DRR000386.bam"
    peak_bed = bam_file.replace('.bam', '.bed')

    input_files = {
        "bam": bam_file,
    }

    output_files = {
        "bed": peak_bed
    }

    metadata = {
        "bam": Metadata(
            data_type='data_mnaseseq',
            file_type="BAM",
            file_path=bam_file,
            sources=[genome_fa, fastq_file],
            taxon_id=10090,
            meta_data={
                "assembly": "GRCm38",
                "tool": "bwa_aligner"
            }
        )
    }

    inps_obj = inps.inps()
    inps_files, inps_meta = inps_obj.run(
        input_files,
        metadata,
        output_files
    )

    # Add tests for all files created
    for f_out in inps_files:
        print("iNPS RESULTS FILE:", f_out)
        assert os.path.isfile(inps_files[f_out]) is True
        assert os.path.getsize(inps_files[f_out]) > 0
