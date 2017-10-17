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

from tool.kallisto_quant import kallistoQuantificationTool

@pytest.mark.rnaseq
def test_kallisto_quant():
    """
    Function to test Kallisto quantifier
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "cdna": resource_path + "kallisto.Human.GRCh38.fasta",
        "index": resource_path + "kallisto.Human.GRCh38.idx",
        "fastq1": resource_path + "kallisto.Human.ERR030872_1.fastq",
        "fastq2": resource_path + "kallisto.Human.ERR030872_2.fastq"
    }

    output_files = {
        "abundance_h5_file": resource_path + "kallisto.Human.ERR030872.abundance.h5",
        "abundance_tsv_file": resource_path + "kallisto.Human.ERR030872.abundance.tsv",
        "run_info_file": resource_path + "kallisto.Human.ERR030872.run_info.json"
    }

    metadata = {
        "cdna": Metadata(
            "data_cdna", "fasta", [], None,
            {'assembly' : 'test'}),
        "index": Metadata(
            "data_cdna", "fasta", [], None,
            {'assembly' : 'test'}),
        "fastq1": Metadata(
            "data_rnaseq", "fastq", [], None,
            {'assembly' : 'test'}),
        "fastq2": Metadata(
            "data_rnaseq", "fastq", [], None,
            {'assembly' : 'test'}),
    }

    kqft = kallistoQuantificationTool()
    kqft.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["abundance_h5_file"]) is True
    assert os.path.getsize(output_files["abundance_h5_file"]) > 0
    assert os.path.isfile(output_files["abundance_tsv_file"]) is True
    assert os.path.getsize(output_files["abundance_tsv_file"]) > 0
    assert os.path.isfile(output_files["run_info_file"]) is True
    assert os.path.getsize(output_files["run_info_file"]) > 0
