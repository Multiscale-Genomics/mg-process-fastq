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

from tool.validate_fastqc import fastqcTool

@pytest.mark.chipseq
def test_fastqc_chipseq_0():
    """
    Test case to ensure that the GEM indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "fastq": resource_path + "kallisto.Human.ERR030872_1.fastq"
    }

    output_files = {
        "report": resource_path + "kallisto.Human.ERR030872_1.report",
    }

    metadata = {
        "fastq": Metadata(
            "data_rnaseq", "fastq", [], None,
            {'assembly' : 'test'})
    }

    fastqc_handle = fastqcTool()
    fastqc_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["report"]) is True
    assert os.path.getsize(output_files["report"]) > 0

@pytest.mark.hic
def test_fastqc_hic_0():
    """
    Test case to ensure that the GEM indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "fastq": resource_path + "tb.Human.SRR1658573_1.fastq.gz"
    }

    output_files = {
        "report": resource_path + "tb.Human.SRR1658573_1.report",
    }

    metadata = {
        "fastq": Metadata(
            "data_rnaseq", "fastq", [], None,
            {'assembly' : 'test'})
    }

    fastqc_handle = fastqcTool()
    fastqc_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["report"]) is True
    assert os.path.getsize(output_files["report"]) > 0

@pytest.mark.hic
def test_fastqc_hic_1():
    """
    Test case to ensure that the GEM indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "fastq": resource_path + "tb.Human.SRR1658573_2.fastq.gz"
    }

    output_files = {
        "report": resource_path + "tb.Human.SRR1658573_2.report",
    }

    metadata = {
        "fastq": Metadata(
            "data_rnaseq", "fastq", [], None,
            {'assembly' : 'test'})
    }

    fastqc_handle = fastqcTool()
    fastqc_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["report"]) is True
    assert os.path.getsize(output_files["report"]) > 0

@pytest.mark.idamidseq
def test_fastqc_idamseq_0():
    """
    Test case to ensure that the GEM indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "fastq": resource_path + "idear.Human.SRR3714775.fastq.gz"
    }

    output_files = {
        "report": resource_path + "idear.Human.SRR3714775.report",
    }

    metadata = {
        "fastq": Metadata(
            "data_idamseq", "fastq", [], None,
            {'assembly' : 'test'})
    }

    fastqc_handle = fastqcTool()
    fastqc_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["report"]) is True
    assert os.path.getsize(output_files["report"]) > 0

@pytest.mark.idamidseq
def test_fastqc_idamseq_1():
    """
    Test case to ensure that the GEM indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "fastq": resource_path + "idear.Human.SRR3714776.fastq.gz"
    }

    output_files = {
        "report": resource_path + "idear.Human.SRR3714776.report",
    }

    metadata = {
        "fastq": Metadata(
            "data_idamseq", "fastq", [], None,
            {'assembly' : 'test'})
    }

    fastqc_handle = fastqcTool()
    fastqc_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["report"]) is True
    assert os.path.getsize(output_files["report"]) > 0

@pytest.mark.idamidseq
def test_fastqc_idamseq_2():
    """
    Test case to ensure that the GEM indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "fastq": resource_path + "idear.Human.SRR3714777.fastq.gz"
    }

    output_files = {
        "report": resource_path + "idear.Human.SRR3714777.report",
    }

    metadata = {
        "fastq": Metadata(
            "data_idamseq", "fastq", [], None,
            {'assembly' : 'test'})
    }

    fastqc_handle = fastqcTool()
    fastqc_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["report"]) is True
    assert os.path.getsize(output_files["report"]) > 0

@pytest.mark.idamidseq
def test_fastqc_idamseq_3():
    """
    Test case to ensure that the GEM indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "fastq": resource_path + "idear.Human.SRR3714778.fastq.gz"
    }

    output_files = {
        "report": resource_path + "idear.Human.SRR3714778.report",
    }

    metadata = {
        "fastq": Metadata(
            "data_idamseq", "fastq", [], None,
            {'assembly' : 'test'})
    }

    fastqc_handle = fastqcTool()
    fastqc_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["report"]) is True
    assert os.path.getsize(output_files["report"]) > 0

@pytest.mark.inps
def test_fastqc_inps_0():
    """
    Test case to ensure that the GEM indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "fastq": resource_path + "inps.Mouse.DRR000386.fastq"
    }

    output_files = {
        "report": resource_path + "inps.Mouse.DRR000386.report",
    }

    metadata = {
        "fastq": Metadata(
            "data_rnaseq", "fastq", [], None,
            {'assembly' : 'test'})
    }

    fastqc_handle = fastqcTool()
    fastqc_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["report"]) is True
    assert os.path.getsize(output_files["report"]) > 0

@pytest.mark.rnaseq
def test_fastqc_rnaseq_0():
    """
    Test case to ensure that the GEM indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "fastq": resource_path + "kallisto.Human.ERR030872_1.fastq"
    }

    output_files = {
        "report": resource_path + "kallisto.Human.ERR030872_1.report",
    }

    metadata = {
        "fastq": Metadata(
            "data_rnaseq", "fastq", [], None,
            {'assembly' : 'test'})
    }

    fastqc_handle = fastqcTool()
    fastqc_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["report"]) is True
    assert os.path.getsize(output_files["report"]) > 0

@pytest.mark.rnaseq
def test_fastqc_rnaseq_1():
    """
    Test case to ensure that the GEM indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "fastq": resource_path + "kallisto.Human.ERR030872_2.fastq"
    }

    output_files = {
        "report": resource_path + "kallisto.Human.ERR030872_2.report",
    }

    metadata = {
        "fastq": Metadata(
            "data_rnaseq", "fastq", [], None,
            {'assembly' : 'test'})
    }

    fastqc_handle = fastqcTool()
    fastqc_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["report"]) is True
    assert os.path.getsize(output_files["report"]) > 0

@pytest.mark.wgbs
def test_fastqc_wgbs_0():
    """
    Test case to ensure that the GEM indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "fastq": resource_path + "bsSeeker.Mouse.GRCm38_1.fastq"
    }

    output_files = {
        "report": resource_path + "bsSeeker.Mouse.GRCm38_1.report",
    }

    metadata = {
        "fastq": Metadata(
            "data_rnaseq", "fastq", [], None,
            {'assembly' : 'test'})
    }

    fastqc_handle = fastqcTool()
    fastqc_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["report"]) is True
    assert os.path.getsize(output_files["report"]) > 0

@pytest.mark.wgbs
def test_fastqc_wgbs_1():
    """
    Test case to ensure that the GEM indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "fastq": resource_path + "bsSeeker.Mouse.GRCm38_2.fastq"
    }

    output_files = {
        "report": resource_path + "bsSeeker.Mouse.GRCm38_2.report",
    }

    metadata = {
        "fastq": Metadata(
            "data_rnaseq", "fastq", [], None,
            {'assembly' : 'test'})
    }

    fastqc_handle = fastqcTool()
    fastqc_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["report"]) is True
    assert os.path.getsize(output_files["report"]) > 0
