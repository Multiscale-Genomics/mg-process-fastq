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
import pytest  # pylint: disable=unused-import
import mock
import subprocess
import pysam

from tool.bam_utils import bamUtils


def touch(path):
    with open(path, 'a'):
        os.utime(path, None)


@pytest.mark.code
def test_bam_index(mocker):
    """
    Test the bam_index function code
    """
    cmd_view = ' '.join([
        'samtools index',
        '-b',
        'example.bam',
        'example.bam_tmp.bai'
    ])
    mocker.patch('subprocess.Popen')

    touch('example.bam')
    touch('example.bam_tmp.bai')
    result = bamUtils.bam_index('example.bam', 'example.bam.bai')
    subprocess.Popen.assert_called_once_with(cmd_view, shell=True)
    assert result is True


@pytest.mark.code
def test_bam_sort(mocker):
    """
    Test the bam_sort function code
    """
    mocker.patch('pysam.sort')
    result = bamUtils.bam_sort('example.bam')
    pysam.sort.assert_called_once_with('-o', 'example.bam', '-T', 'example.bam' + '_sort', 'example.bam')
    assert result is True


@pytest.mark.code
def test_bam_merge_list(mocker):
    """
    """
    mocker.patch('pysam.merge')
    touch('example_1.bam')
    touch('example_2.bam')
    result = bamUtils.bam_merge(['example_1.bam', 'example_2.bam'])
    pysam.merge.assert_called_once_with('-f', 'example_1.bam_merge.bam', 'example_1.bam', 'example_2.bam')
    assert result is False

@pytest.mark.code
def test_bam_merge(mocker):
    """
    """
    mocker.patch('pysam.merge')
    touch('example_1.bam')
    touch('example_2.bam')
    result = bamUtils.bam_merge('example_1.bam', 'example_2.bam')
    pysam.merge.assert_called_once_with('-f', 'example_1.bam_merge.bam', 'example_1.bam', 'example_2.bam')
    assert result is False
