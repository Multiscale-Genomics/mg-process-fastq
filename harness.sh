#!/bin/bash

python tests/test_toolchains.py --pipeline genome --verbose 1
tc1=$?
ls tests/data/
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
tc2=$?
ls tests/data/
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
tc3=$?
ls tests/data/
rm tests/data/abundance.h5
rm tests/data/abundance.tsv
rm tests/data/run_info.json
rm tests/data/kallisto.Human.ERR030872.abundance.h5
rm tests/data/kallisto.Human.ERR030872.abundance.tsv
rm tests/data/kallisto.Human.ERR030872.run_info.json
rm tests/data/kallisto.Human.GRCh38.idx

python tests/test_toolchains.py --pipeline wgbs
tc4=$?
ls tests/data/
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

rc=$(($tc1+$tc2+$tc3+$tc4))

if [[ $rc != 0 ]]; then exit $rc; fi

# pytest  tests/test_bwa_indexer.py
# pytest  tests/test_bwa_aligner.py
# pytest  tests/test_biobambam.py
# pytest  tests/test_bs_seeker_indexer.py
# pytest  tests/test_bs_seeker_aligner.py
# pytest  tests/test_bs_seeker_filter.py
# samtools
# ls /home/travis/build/Multiscale-Genomics/mg-process-fastq/tests/data
# pytest  tests/test_bs_seeker_methylation_caller.py
# pytest -s tests/test_bowtie_indexer.py
# pytest -s tests/test_kallisto_indexer.py
# pytest -s tests/test_kallisto_quant.py
# pytest -s tests/test_macs2.py
# pytest -s tests/test_paired_splitter.py
# pytest -s tests/test_single_splitter.py

# python_version=$(python --version 2>&1)

# echo $python_version
# if [[ $python_version == *"3."* ]]; then
#     pytest -s tests/test_inps.py
# fi
