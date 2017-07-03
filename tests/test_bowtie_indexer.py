#!usr/bin/python

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

import random
import os.path
import pytest

from tool import bowtie_indexer


def test_bowtie_indexer():
    """
    Test to ensure Bowtie indexer is working for macs data set
    """
    bti = bowtie_indexer.bowtieIndexerTool()
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    bti.run([resource_path + "macs2.Human.GCA_000001405.22.fasta"], {})


def test_bowtie_indexer_02():
    """
    Test to ensure Bowtie indexer is working for macs data set
    """
    bti = bowtie_indexer.bowtieIndexerTool()
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    bti.run([resource_path + "inps.Mouse.GRCm38.fasta"], {})
