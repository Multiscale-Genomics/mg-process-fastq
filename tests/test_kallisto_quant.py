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
import pytest

from basic_modules.metadata import Metadata

from tool.kallisto_quant import kallistoQuantificationTool


@pytest.mark.rnaseq
def test_kallisto_quant_paired():
    """
    Function to test Kallisto quantifier
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "cdna": resource_path + "kallisto.Human.GRCh38.fasta",
        "index": resource_path + "kallisto.Human.GRCh38.idx",
        "fastq1": resource_path + "kallisto.Human.ERR030872_1.fastq",
        "fastq2": resource_path + "kallisto.Human.ERR030872_2.fastq",
        "gff": resource_path + "kallisto.Human.GRCh38.gff3"
    }

    output_files = {
        "abundance_h5_file": resource_path + "kallisto.Human.ERR030872.paired.abundance.h5",
        "abundance_tsv_file": resource_path + "kallisto.Human.ERR030872.paired.abundance.tsv",
        "abundance_gff_file": resource_path + "kallisto.Human.ERR030872.paired.abundance.gff",
        "run_info_file": resource_path + "kallisto.Human.ERR030872.paired.run_info.json"
    }

    metadata = {
        "cdna": Metadata(
            "data_cdna", "fasta", [], None,
            {"assembly": "GCA_000001405.22", "ensembl": True}),
        "index": Metadata(
            "data_cdna", "fasta", [], None,
            {"assembly": "GCA_000001405.22", "ensembl": True}),
        "fastq1": Metadata(
            "data_rnaseq", "fastq", [], None,
            {"assembly": "GCA_000001405.22", "ensembl": True}),
        "fastq2": Metadata(
            "data_rnaseq", "fastq", [], None,
            {"assembly": "GCA_000001405.22", "ensembl": True}),
        "gff": Metadata(
            "data_seq", "gff", [], None,
            {"assembly": "GCA_000001405.22", "ensembl": True}),
    }

    kqft = kallistoQuantificationTool({"execution": resource_path})
    rs_files, rs_meta = kqft.run(input_files, metadata, output_files)

    # Checks that the returned files matches the expected set of results
    assert len(rs_meta) == 4

    # Add tests for all files created
    for f_out in rs_files:
        print("RNA SEQ RESULTS FILE:", f_out)
        assert rs_files[f_out] == output_files[f_out]
        assert os.path.isfile(rs_files[f_out]) is True
        assert os.path.getsize(rs_files[f_out]) > 0
        os.remove(rs_files[f_out])


@pytest.mark.rnaseq
def test_kallisto_quant_single():
    """
    Function to test Kallisto quantifier
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "cdna": resource_path + "kallisto.Human.GRCh38.fasta",
        "index": resource_path + "kallisto.Human.GRCh38.idx",
        "fastq1": resource_path + "kallisto.Human.ERR030872_1.fastq",
        "gff": resource_path + "kallisto.Human.GRCh38.gff3"
    }

    output_files = {
        "abundance_h5_file": resource_path + "kallisto.Human.ERR030872.single.abundance.h5",
        "abundance_tsv_file": resource_path + "kallisto.Human.ERR030872.single.abundance.tsv",
        "abundance_gff_file": resource_path + "kallisto.Human.ERR030872.single.abundance.gff",
        "run_info_file": resource_path + "kallisto.Human.ERR030872.single.run_info.json"
    }

    metadata = {
        "cdna": Metadata(
            "data_cdna", "fasta", [], None,
            {"assembly": "GCA_000001405.22", "ensembl": True}),
        "index": Metadata(
            "data_cdna", "fasta", [], None,
            {"assembly": "GCA_000001405.22", "ensembl": True}),
        "fastq1": Metadata(
            "data_rnaseq", "fastq", [], None,
            {"assembly": "GCA_000001405.22", "ensembl": True}),
        "gff": Metadata(
            "data_seq", "gff", [], None,
            {"assembly": "GCA_000001405.22", "ensembl": True}),
    }

    kqft = kallistoQuantificationTool({"execution": resource_path})
    rs_files, rs_meta = kqft.run(input_files, metadata, output_files)

    # Checks that the returned files matches the expected set of results
    assert len(rs_meta) == 4

    # Add tests for all files created
    for f_out in rs_files:
        print("RNA SEQ RESULTS FILE:", f_out)
        assert rs_files[f_out] == output_files[f_out]
        assert os.path.isfile(rs_files[f_out]) is True
        assert os.path.getsize(rs_files[f_out]) > 0
        os.remove(rs_files[f_out])
