#!/bin/bash

# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

rc=0

python tests/test_toolchains.py --pipeline genome --verbose 1
tc=$?
rc=$(($rc + $tc))
# ls tests/data/
rm -rf tests/data/inps.Mouse.GRCm38.fasta.bt2
rm -rf tests/data/inps.Mouse.GRCm38.fasta.bwa
rm tests/data/inps.Mouse.GRCm38.fasta.bt2.tar.gz
rm tests/data/inps.Mouse.GRCm38.fasta.bwa.tar.gz
rm -rf tests/data/macs2.Human.GCA_000001405.22.fasta.bt2
rm -rf tests/data/macs2.Human.GCA_000001405.22.fasta.bwa
rm tests/data/macs2.Human.GCA_000001405.22.fasta.bt2.tar.gz
rm tests/data/macs2.Human.GCA_000001405.22.fasta.bwa.tar.gz
rm tests/data/tb.Human.GCA_000001405.22.fasta
rm tests/data/tb.Human.GCA_000001405.22_gem*

python tests/test_toolchains.py --pipeline chipseq
tc=$?
rc=$(($rc + $tc))
# ls tests/data/
rm tests/data/macs2.Human.DRR000150.22.bam
rm tests/data/macs2.Human.DRR000150.22.bam.filtered.tmp.bam
rm tests/data/macs2.Human.DRR000150.22.fastq.out.bam
rm tests/data/macs2.Human.DRR000150.22.fastq.sai
rm tests/data/macs2.Human.DRR000150.22.fastq.sam
rm tests/data/macs2.Human.DRR000150.22_filtered.bam
rm tests/data/macs2.Human.DRR000150.22_filtered_out_peaks.narrowPeak
rm tests/data/macs2.Human.DRR000150.22_filtered_out_peaks.xls
rm tests/data/macs2.Human.DRR000150.22_filtered_out_summits.bed
rm tests/data/macs2.Human.DRR000150.22_peaks.broadPeak
rm tests/data/macs2.Human.DRR000150.22_peaks.gappedPeak
rm tests/data/macs2.Human.DRR000150.22_peaks.narrowPeak
rm tests/data/macs2.Human.DRR000150.22_peaks.summits.bed
rm tests/data/macs2.Human.GCA_000001405.22.fasta.bwa.tar.gz
rm -r tests/data/macs2.Human.GCA_000001405.22.fasta.bwa

python tests/test_toolchains.py --pipeline rnaseq
tc=$?
rc=$(($rc + $tc))
# ls tests/data/
rm tests/data/abundance.h5
rm tests/data/abundance.tsv
rm tests/data/run_info.json
rm tests/data/kallisto.Human.ERR030872.abundance.h5
rm tests/data/kallisto.Human.ERR030872.abundance.tsv
rm tests/data/kallisto.Human.ERR030872.run_info.json
rm tests/data/kallisto.Human.GRCh38.idx

if [[ $python_version != *"3."* ]]; then
    python tests/test_toolchains.py --pipeline wgbs
    tc=$?
    rc=$(($rc + $tc))
    # ls tests/data/
    rm tests/data/bsSeeker.Mouse.GRCm38.fasta.bt2.tar.gz
    rm -r tests/data/bsSeeker.Mouse.GRCm38.fasta_bowtie2/
    rm tests/data/bsSeeker.Mouse.GRCm38_1.atcgmap
    rm tests/data/bsSeeker.Mouse.GRCm38_1.cgmap
    rm tests/data/bsSeeker.Mouse.GRCm38_1.wig
    rm tests/data/bsSeeker.Mouse.GRCm38_1_filtered.bai
    rm tests/data/bsSeeker.Mouse.GRCm38_1_filtered.bam
    rm tests/data/bsSeeker.Mouse.GRCm38_1_filtered.bam.bai
    rm tests/data/bsSeeker.Mouse.GRCm38_1_filtered.bam.call_methylation_log
    rm tests/data/bsSeeker.Mouse.GRCm38_1_filtered.bam_tmp.bai
    rm tests/data/bsSeeker.Mouse.GRCm38_1_filtered.fastq
    rm tests/data/bsSeeker.Mouse.GRCm38_1_filtered.fastq.tar.gz
    rm tests/data/bsSeeker.Mouse.GRCm38_1_filtered.fastq.tmp
    rm tests/data/bsSeeker.Mouse.GRCm38_2_filtered.fastq
    rm tests/data/bsSeeker.Mouse.GRCm38_2_filtered.fastq.tmp
fi

# if [[ $python_version == *"3."* ]]; then
#     python tests/test_toolchains.py --pipeline mnaseseq
#     tc=$?
#     rc=$(($rc + $tc))
#     ls tests/data/
# fi

if [[ $rc != 0 ]]; then exit $rc; fi
