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

from os import path, remove
import pytest # pylint: disable=unused-import

from basic_modules.metadata import Metadata
from tadbit_map_parse_filter_wrapper import tadbit_map_parse_filter
from tadbit_normalize_wrapper import tadbit_normalize
from tadbit_segment_wrapper import tadbit_segment
from tadbit_bin_wrapper import tadbit_bin

@pytest.mark.hic
@pytest.mark.pipeline
def test_tadbit_tools():
    """
    Test case to ensure that the tadbit tools code works.
    """
    
    resource_path = path.join(path.dirname(__file__), "data/")
    
    files = {
        "parsing:refGenome": resource_path + 'tb.Human.GCA_000001405.22_gem.fasta',
        "mapping:refGenome": resource_path + 'tb.Human.GCA_000001405.22_gem.fasta.gem',
        "read1": resource_path + 'tb.Human.SRR1658573_1.fastq',
        "read2": resource_path + 'tb.Human.SRR1658573_2.fastq'
    }
    metadata_1 = {}
    metadata_1["read1"] = Metadata(
                    data_type = "hic_reads",
                    file_type = "fastq",
                    file_path=files["read1"],
                    sources=[""],
                    meta_data={
                        "visible": True,
                        "assembly": "hg38"
                    },
                    taxon_id=9606)

    metadata_1["read2"] = Metadata(
                    data_type = "hic_reads",
                    file_type = "fastq",
                    file_path=files["read1"],
                    sources=[""],
                    meta_data={
                        "visible": True,
                        "assembly": "hg38"
                    },
                    taxon_id=9606)

    metadata_1["parsing:refGenome"] = Metadata(
                    data_type = "sequence_genomic",
                    file_type = "fasta",
                    file_path=files["parsing:refGenome"],
                    sources=[""],
                    meta_data={
                        "visible": True,
                        "assembly": "hg38"
                    },
                    taxon_id=9606)

    metadata_1["mapping:refGenome"] = Metadata(
                    data_type = "sequence_mapping_index",
                    file_type = "gem",
                    file_path=files["mapping:refGenome"],
                    sources=[""],
                    meta_data={
                        "visible": True,
                        "assembly": "hg38"
                    },
                    taxon_id=9606)

    config = {
        "project": resource_path,
        "execution": resource_path,
        "description": "TADbit tools test",
        "mapping:rest_enzyme": "MboI",
        "mapping:iterative_mapping": False,
        "filtering:filters": [
                "1",
                "2",
                "3",
                "4",
                "9",
                "10"
             ],
        "chromosomes": "",
        "mapping:windows": "1:20 1:40",
        "filtering:min_dist_RE": "500",
        "filtering:min_fragment_size": "50",
        "filtering:max_fragment_size": "100000"
    }
    if path.isfile(resource_path + "/data/paired_reads.bam"):
        remove(resource_path + "/data/paired_reads.bam")
    if path.isfile(resource_path + "/data/paired_reads.bam.bai"):
        remove(resource_path + "/data/paired_reads.bam.bai")
    if path.isfile(resource_path + "/data/map_parse_filter_stats.tar.gz"):
        remove(resource_path + "/data/map_parse_filter_stats.tar.gz")
    
    tb_handle = tadbit_map_parse_filter(configuration=config)
    tb_1_files, tb_1_meta = tb_handle.run(files, metadata_1, [])

    print(tb_1_files)

    for tb_file in tb_1_files:
        assert path.isfile(tb_1_files[tb_file]) is True
        assert path.getsize(tb_1_files[tb_file]) > 0

    files = {
        "bamin": tb_1_files["paired_reads"]
    }
    
    metadata_2 = {}
    metadata_2["bamin"] = tb_1_meta["paired_reads"]
    config = {
        "project": resource_path,
        "execution": resource_path,
        "description": "TADbit tools test",
        "resolution": "100000",
        "segmentation:chromosome_names": "",
        "segmentation:callers": [
                "1",
                "2"
            ]
    }
    tb_handle = tadbit_segment(configuration=config)
    tb_2_files, tb_2_meta = tb_handle.run(files, metadata_2, [])

    print(tb_2_files)

    for tb_file in tb_2_files:
        assert path.isfile(tb_2_files[tb_file]) is True
        assert path.getsize(tb_2_files[tb_file]) > 0
    
    files = {
        "bamin": tb_1_files["paired_reads"],
        "refGenomes_folder": resource_path + 'tb.Human.GCA_000001405.22_gem.fasta'
    }
    metadata_3 = {}
    metadata_3["bamin"] = tb_1_meta["paired_reads"]
    config = {
        "project": resource_path,
        "execution": resource_path,
        "description": "TADbit tools test",
        "resolution": "100000",
        "min_perc": "2",
        "max_perc": "99.8",
        "normalization": "Vanilla"
    }
    tb_handle = tadbit_normalize(configuration=config)
    tb_3_files, tb_3_meta = tb_handle.run(files, metadata_3, [])
 
    print(tb_3_files)
 
    for tb_file in tb_3_files:
        assert path.isfile(tb_3_files[tb_file]) is True
        assert path.getsize(tb_3_files[tb_file]) > 0
        
    files = {
        "bamin": tb_1_files["paired_reads"],
        "hic_biases": tb_3_files["hic_biases"]
    }
    
    metadata_4 = {}
    metadata_4["bamin"] = tb_1_meta["paired_reads"]
    metadata_4["hic_biases"] = tb_3_meta["hic_biases"]
    config = {
        "project": resource_path,
        "execution": resource_path,
        "description": "TADbit tools test",
        "resolution": "100000",
        "coord1": ""
    }
    tb_handle = tadbit_bin(configuration=config)
    tb_4_files, tb_4_meta = tb_handle.run(files, metadata_4, [])

    print(tb_4_files)

    for tb_file in tb_4_files:
        assert path.isfile(tb_4_files[tb_file]) is True
        assert path.getsize(tb_4_files[tb_file]) > 0