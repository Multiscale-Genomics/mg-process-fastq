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

import os
import shutil
import gzip
import subprocess
import pytest

from basic_modules.metadata import Metadata

from mg_process_fastq.tool.bowtie_aligner import bowtie2AlignerTool


@pytest.mark.bowtie2
def test_bowtie2_aligner_00():
    """
    Extract the compressed FASTQ files
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    fastq_file_1 = resource_path + "bsSeeker.Mouse.SRR892982_1.fastq"
    fastq_file_2 = resource_path + "bsSeeker.Mouse.SRR892982_2.fastq"

    files = [fastq_file_1, fastq_file_2]

    for fastq_gz in files:
        with gzip.open(fastq_gz + '.gz', 'rb') as fgz_in:
            with open(fastq_gz, 'w') as f_out:
                f_out.write(fgz_in.read())

    assert os.path.isfile(fastq_file_1) is True
    assert os.path.getsize(fastq_file_1) > 0
    assert os.path.isfile(fastq_file_2) is True
    assert os.path.getsize(fastq_file_2) > 0


@pytest.mark.bowtie2
def test_bowtie2_aligner_single():
    """
    Function to test Bowtie Aligner
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genome_fa = resource_path + "macs2.Human.GCA_000001405.22.fasta"
    fastq_file = resource_path + "macs2.Human.DRR000150.22.fastq"

    input_files = {
        "genome": genome_fa,
        "index": genome_fa + ".bt2.tar.gz",
        "loc": fastq_file
    }

    output_files = {
        "output": fastq_file.replace(".fastq", "_bt2.bam")
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

    bowtie2_handle = bowtie2AlignerTool()
    bowtie2_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(resource_path + "macs2.Human.DRR000150.22_bt2.bam") is True
    assert os.path.getsize(resource_path + "macs2.Human.DRR000150.22_bt2.bam") > 0

    cmdl = "samtools view -c -f 0 {}macs2.Human.DRR000150.22_bt2.bam".format(resource_path)
    process = subprocess.Popen(cmdl, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.wait()
    proc_out, proc_err = process.communicate()  # pylint: disable=unused-variable
    assert int(proc_out.strip()) > 0

    try:
        os.remove(resource_path + "macs2.Human.DRR000150.22_bt2.bam")
    except OSError, ose:
        print("Error: %s - %s." % (ose.filename, ose.strerror))


@pytest.mark.bowtie2
def test_bowtie2_aligner_paired():
    """
    Function to test Bowtie Aligner
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genome_fa = resource_path + "bsSeeker.Mouse.GRCm38.fasta"
    fastq_file_1 = resource_path + "bsSeeker.Mouse.SRR892982_1.fastq"
    fastq_file_2 = resource_path + "bsSeeker.Mouse.SRR892982_2.fastq"

    input_files = {
        "genome": genome_fa,
        "index": genome_fa + ".bt2.tar.gz",
        "loc": fastq_file_1,
        "fastq2": fastq_file_2
    }

    output_files = {
        "output": fastq_file_1.replace(".fastq", "_bt2.bam")
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
            "data_wgbs", "fastq", fastq_file_1, None,
            {"assembly": "test"}
        ),
        "fastq2": Metadata(
            "data_wgbs", "fastq", fastq_file_2, None,
            {"assembly": "test"}
        )
    }

    bowtie2_handle = bowtie2AlignerTool()
    bowtie2_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(resource_path + "bsSeeker.Mouse.SRR892982_1_bt2.bam") is True
    assert os.path.getsize(resource_path + "bsSeeker.Mouse.SRR892982_1_bt2.bam") > 0

    cmdl = "samtools view -c -f 0 {}bsSeeker.Mouse.SRR892982_1_bt2.bam".format(resource_path)
    process = subprocess.Popen(cmdl, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.wait()
    proc_out, proc_err = process.communicate()  # pylint: disable=unused-variable
    assert int(proc_out.strip()) > 0

    try:
        os.remove(resource_path + "bsSeeker.Mouse.SRR892982_1_bt2.bam")
    except OSError, ose:
        print("Error: %s - %s." % (ose.filename, ose.strerror))

    try:
        shutil.rmtree(resource_path + "tmp")
    except OSError, ose:
        print("Error: %s - %s." % (ose.filename, ose.strerror))
