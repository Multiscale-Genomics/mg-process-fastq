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

from tool.bwa_aligner import bwaAlignerTool

@pytest.mark.chipseq
def test_bwa_aligner_chipseq():
    """
    Function to test BWA Aligner
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genome_fa = resource_path + "macs2.Human.GCA_000001405.22.fasta"
    fastq_file = resource_path + "macs2.Human.DRR000150.22.fastq"

    input_files = {
        "genome": genome_fa,
        "index": genome_fa + ".bwa.tar.gz",
        "loc": fastq_file
    }

    output_files = {
        "output": fastq_file.replace(".fastq", ".bam")
    }

    metadata = {
        "genome": Metadata(
            "Assembly", "fasta", genome_fa, None,
            {"assembly": "test"}),
        "index": Metadata(
            "index_bwa", "", [genome_fa],
            {
                "assembly": "test",
                "tool": "bwa_indexer"
            }
        ),
        "loc": Metadata(
            "data_chip_seq", "fastq", fastq_file, None,
            {"assembly": "test"}
        )
    }

    bwa_t = bwaAlignerTool()
    bwa_t.run(input_files, metadata, output_files)

    print(__file__)

    assert os.path.isfile(resource_path + "macs2.Human.DRR000150.22.bam") is True
    assert os.path.getsize(resource_path + "macs2.Human.DRR000150.22.bam") > 0

@pytest.mark.idamidseq
def test_bwa_aligner_idamidseq():
    """
    Function to test BWA Aligner
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genome_fa = resource_path + "idear.Human.GCA_000001405.22.fasta"
    fastq_files = [
        resource_path + "idear.Human.SRR3714775.fastq",
        resource_path + "idear.Human.SRR3714776.fastq",
        resource_path + "idear.Human.SRR3714777.fastq",
        resource_path + "idear.Human.SRR3714778.fastq"
    ]

    for fastq_file in fastq_files:
        print(fastq_file.replace(".fastq", ".bam"))

        input_files = {
            "genome": genome_fa,
            "index": genome_fa + ".bwa.tar.gz",
            "loc": fastq_file
        }

        output_files = {
            "output": fastq_file.replace(".fastq", ".bam")
        }

        metadata = {
            "genome": Metadata(
                "Assembly", "fasta", genome_fa, None,
                {"assembly": "test"}),
            "index": Metadata(
                "index_bwa", "", [genome_fa],
                {
                    "assembly": "test",
                    "tool": "bwa_indexer"
                }
            ),
            "loc": Metadata(
                "data_damid_seq", "fastq", fastq_file, None,
                {"assembly": "test"}
            )
        }

        bwa_t = bwaAlignerTool()
        bwa_t.run(input_files, metadata, output_files)

        assert os.path.isfile(fastq_file.replace(".fastq", ".bam")) is True
        assert os.path.getsize(fastq_file.replace(".fastq", ".bam")) > 0

@pytest.mark.mnaseseq
def test_bwa_aligner_mnaseseq():
    """
    Function to test BWA Aligner for MNase seq data
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genome_fa = resource_path + "inps.Mouse.GRCm38.fasta"
    fastq_file = resource_path + "inps.Mouse.DRR000386.fastq"

    input_files = {
        "genome": genome_fa,
        "index": genome_fa + ".bwa.tar.gz",
        "loc": fastq_file
    }

    output_files = {
        "output": fastq_file.replace(".fastq", ".bam")
    }

    metadata = {
        "genome": Metadata(
            "Assembly", "fasta", genome_fa, None,
            {"assembly": "test"}),
        "index": Metadata(
            "index_bwa", "", [genome_fa],
            {
                "assembly": "test",
                "tool": "bwa_indexer"
            }
        ),
        "loc": Metadata(
            "data_chip_seq", "fastq", fastq_file, None,
            {"assembly": "test"}
        )
    }

    bwa_t = bwaAlignerTool()
    bwa_t.run(input_files, metadata, output_files)

    print(__file__)

    assert os.path.isfile(resource_path + "inps.Mouse.DRR000386.bam") is True
    assert os.path.getsize(resource_path + "inps.Mouse.DRR000386.bam") > 0
