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
from tool.kallisto_indexer import kallistoIndexerTool


@pytest.mark.rnaseq
def test_kallisto_indexer_00():
    """
    Function to test Kallisto indexer
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "cdna": resource_path + "kallisto.Human.GRCh38.fasta"
    }

    output_files = {
        "index": resource_path + "kallisto.Human.GRCh38.idx"
    }

    metadata = {
        "cdna": Metadata(
            "data_cdna", "fasta", [], None,
            {'assembly': 'test'}),
    }

    ki_handle = kallistoIndexerTool({"execution": resource_path})
    ki_handle.run(input_files, metadata, output_files)

    print(__file__)

    assert os.path.isfile(resource_path + "kallisto.Human.GRCh38.idx") is True
    assert os.path.getsize(resource_path + "kallisto.Human.GRCh38.idx") > 0


@pytest.mark.sleuth
def test_kallisto_indexer_01():
    """
    Function to test Kallisto indexer
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/sleuth/")

    input_files = {
        "cdna": resource_path + "sleuth.Human.GRCh38.cdna.fasta"
    }

    output_files = {
        "index": resource_path + "sleuth.Human.GRCh38.cdna.idx"
    }

    metadata = {
        "cdna": Metadata(
            "data_cdna", "fasta", [], None,
            {'assembly': 'test'}),
    }

    ki_handle = kallistoIndexerTool()
    ki_handle.run(input_files, metadata, output_files)

    print(__file__)

    assert os.path.isfile(resource_path + "sleuth.Human.GRCh38.cdna.idx") is True
    assert os.path.getsize(resource_path + "sleuth.Human.GRCh38.cdna.idx") > 0
