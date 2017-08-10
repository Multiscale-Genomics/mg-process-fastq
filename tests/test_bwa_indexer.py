"""
.. Copyright 2017 EMBL-European Bioinformatics Institute

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

from tool.bwa_indexer import bwaIndexerTool

@pytest.mark.chipseq
@pytest.mark.genome
def test_bwa_indexer():
    """
    Test case to ensure that the BWA indexer works.
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genome_fa = resource_path + "macs2.Human.GCA_000001405.22.fasta"

    files = [
        genome_fa,
        genome_fa + ".amb",
        genome_fa + ".ann",
        genome_fa + ".bwt",
        genome_fa + ".pac",
        genome_fa + ".sa"
    ]

    bwa_amb = files[1]
    bwa_ann = files[2]
    bwa_bwt = files[3]
    bwa_pac = files[4]
    bwa_sa = files[5]

    print(bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_sa)

    bwa_it = bwaIndexerTool()
    bwa_it.run([genome_fa], {'assembly' : 'test'})

    assert os.path.isfile(resource_path + "macs2.Human.GCA_000001405.22.fasta.amb") is True
    assert os.path.getsize(resource_path + "macs2.Human.GCA_000001405.22.fasta.amb") > 0
    assert os.path.isfile(resource_path + "macs2.Human.GCA_000001405.22.fasta.ann") is True
    assert os.path.getsize(resource_path + "macs2.Human.GCA_000001405.22.fasta.ann") > 0
    assert os.path.isfile(resource_path + "macs2.Human.GCA_000001405.22.fasta.bwt") is True
    assert os.path.getsize(resource_path + "macs2.Human.GCA_000001405.22.fasta.bwt") > 0
    assert os.path.isfile(resource_path + "macs2.Human.GCA_000001405.22.fasta.pac") is True
    assert os.path.getsize(resource_path + "macs2.Human.GCA_000001405.22.fasta.pac") > 0
    assert os.path.isfile(resource_path + "macs2.Human.GCA_000001405.22.fasta.sa") is True
    assert os.path.getsize(resource_path + "macs2.Human.GCA_000001405.22.fasta.sa") > 0


@pytest.mark.mnaseseq
@pytest.mark.genome
def test_bwa_indexer_02():
    """
    Test case to ensure that the BWA indexer works
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    genome_fa = resource_path + "inps.Mouse.GRCm38.fasta"

    files = [
        genome_fa,
        genome_fa + ".amb",
        genome_fa + ".ann",
        genome_fa + ".bwt",
        genome_fa + ".pac",
        genome_fa + ".sa"
    ]

    bwa_amb = files[1]
    bwa_ann = files[2]
    bwa_bwt = files[3]
    bwa_pac = files[4]
    bwa_sa = files[5]

    print(bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_sa)

    bwa_it = bwaIndexerTool()
    bwa_it.run([genome_fa], {'assembly' : 'test'})

    assert os.path.isfile(resource_path + "inps.Mouse.GRCm38.fasta.amb") is True
    assert os.path.getsize(resource_path + "inps.Mouse.GRCm38.fasta.amb") > 0
    assert os.path.isfile(resource_path + "inps.Mouse.GRCm38.fasta.ann") is True
    assert os.path.getsize(resource_path + "inps.Mouse.GRCm38.fasta.ann") > 0
    assert os.path.isfile(resource_path + "inps.Mouse.GRCm38.fasta.bwt") is True
    assert os.path.getsize(resource_path + "inps.Mouse.GRCm38.fasta.bwt") > 0
    assert os.path.isfile(resource_path + "inps.Mouse.GRCm38.fasta.pac") is True
    assert os.path.getsize(resource_path + "inps.Mouse.GRCm38.fasta.pac") > 0
    assert os.path.isfile(resource_path + "inps.Mouse.GRCm38.fasta.sa") is True
    assert os.path.getsize(resource_path + "inps.Mouse.GRCm38.fasta.sa") > 0
