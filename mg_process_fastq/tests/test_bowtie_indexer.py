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

import os.path
import pytest

from basic_modules.metadata import Metadata

from mg_process_fastq.tool import bowtie_indexer


@pytest.mark.chipseq
@pytest.mark.genome
@pytest.mark.bowtie2
def test_bowtie_indexer_chipseq():
    """
    Test to ensure Bowtie indexer is working for macs data set
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    genome_fa = resource_path + "macs2.Human.GCA_000001405.22.fasta"

    input_files = {
        "genome": genome_fa
    }

    output_files = {
        "index": genome_fa + ".bt2.tar.gz"
    }

    metadata = {
        "genome": Metadata(
            "Assembly", "fasta", genome_fa, None,
            {'assembly': 'test'}),
    }

    bti = bowtie_indexer.bowtieIndexerTool({"execution": resource_path})
    bti.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["index"]) is True
    assert os.path.getsize(output_files["index"]) > 0


@pytest.mark.mnaseseq
@pytest.mark.genome
def test_bowtie_indexer_mnaseseq():
    """
    Test to ensure Bowtie indexer is working for macs data set
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    genome_fa = resource_path + "inps.Mouse.GRCm38.fasta"

    input_files = {
        "genome": genome_fa
    }

    output_files = {
        "index": genome_fa + ".bt2.tar.gz"
    }

    metadata = {
        "genome": Metadata(
            "Assembly", "fasta", genome_fa, None,
            {'assembly': 'test'}),
    }

    bti = bowtie_indexer.bowtieIndexerTool({"execution": resource_path})
    bti.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["index"]) is True
    assert os.path.getsize(output_files["index"]) > 0


@pytest.mark.genome
@pytest.mark.bowtie2
def test_bowtie_indexer_wgbs():
    """
    Test to ensure Bowtie indexer is working for macs data set
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    genome_fa = resource_path + "bsSeeker.Mouse.GRCm38.fasta"

    input_files = {
        "genome": genome_fa
    }

    output_files = {
        "index": genome_fa + ".bt2.tar.gz"
    }

    metadata = {
        "genome": Metadata(
            "Assembly", "fasta", genome_fa, None,
            {'assembly': 'test'}),
    }

    bti = bowtie_indexer.bowtieIndexerTool({"execution": resource_path})
    bti.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["index"]) is True
    assert os.path.getsize(output_files["index"]) > 0
