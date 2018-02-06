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

import os
import pytest # pylint: disable=unused-import

from basic_modules.metadata import Metadata

from tool import bs_seeker_methylation_caller

@pytest.mark.wgbs
def test_bs_seeker_methylation_caller():
    """
    Test that it is possible to call the methylation called by BS seeker
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")
    home = os.path.expanduser('~')

    input_files = {
        "genome": resource_path + "bsSeeker.Mouse.GRCm38.fasta",
        "index": resource_path + "bsSeeker.Mouse.GRCm38.fasta.bt2.tar.gz",
        "bam": resource_path + "bsSeeker.Mouse.SRR892982_1.filtered.bam",
        "bai": resource_path + "bsSeeker.Mouse.SRR892982_1.filtered.bai",
    }

    output_files = {
        "wig_file": resource_path + "bsSeeker.Mouse.SRR892982_1.wig",
        "cgmap_file": resource_path + "bsSeeker.Mouse.SRR892982_1.cgmap",
        "atcgmap_file": resource_path + "bsSeeker.Mouse.SRR892982_1.atcgmap"
    }

    metadata = {
        "genome": Metadata(
            "Assembly", "fasta", input_files["genome"], None,
            {'assembly' : 'test'}),
        "index": Metadata(
            "index_bowtie", "index", input_files["genome"], None,
            {'assembly' : 'test'}),
        "bam": Metadata(
            "data_wgbs", "bam", input_files["bam"], None,
            {'assembly' : 'test'}),
        "bai": Metadata(
            "data_wgbs", "bai", input_files["bai"], None,
            {'assembly' : 'test'}),
    }

    config_param = {
        "aligner" : "bowtie2",
        "aligner_path" : home + "/lib/bowtie2-2.3.4-linux-x86_64",
        "bss_path" : home + "/lib/BSseeker2"
    }

    bsmc = bs_seeker_methylation_caller.bssMethylationCallerTool(config_param)
    bsmc.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["wig_file"]) is True
    assert os.path.getsize(output_files["wig_file"]) > 0
    assert os.path.isfile(output_files["cgmap_file"]) is True
    assert os.path.getsize(output_files["cgmap_file"]) > 0
    assert os.path.isfile(output_files["atcgmap_file"]) is True
    assert os.path.getsize(output_files["atcgmap_file"]) > 0
