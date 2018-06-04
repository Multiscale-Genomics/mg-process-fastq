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

# libtbb-dev did not seem to be installing correctly without using sudo
sudo apt-get install libtbb-dev

# FastQC
cd ${HOME}/lib
if [ ! -d "FastQC" ]; then
    wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
    unzip fastqc_v0.11.5.zip
    cd FastQC/
    chmod 755 fastqc
fi

# htslib
cd ${HOME}/lib
if [ ! -d "htslib" ]; then
    git clone https://github.com/samtools/htslib.git
    cd htslib
    autoheader
    autoconf
    ./configure --prefix=${HOME}/lib/htslib
    make
    make install
fi

# SAMtools
cd ${HOME}/lib
if [ ! -d "samtools" ]; then
    git clone https://github.com/samtools/samtools.git
    cd samtools
    autoheader
    autoconf -Wno-syntax
    ./configure --prefix=${HOME}/lib/samtools
    make
    make install
fi

# UCSC Tools
cd ${HOME}/lib
if [ ! -f "bedToBigBed" ]; then
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
fi

# Bowtie2 Aligner
cd ${HOME}/lib
if [ ! -d "bowtie2-2.3.4-linux-x86_64" ]; then
    wget --max-redirect 1 https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4/bowtie2-2.3.4-linux-x86_64.zip
    unzip bowtie2-2.3.4-linux-x86_64.zip
fi

# BWA Sequence Aligner
cd ${HOME}/lib
if [ ! -d "bwa" ]; then
    git clone https://github.com/lh3/bwa.git
    cd bwa
    make
fi

# GEM Indexer
cd ${HOME}/lib
if [ ! -d "GEMTools-static-core2-1.7.1" ]; then
    wget http://barnaserver.com/gemtools/releases/GEMTools-static-core2-1.7.1.tar.gz
    tar -xzf GEMTools-static-core2-1.7.1.tar.gz
fi

# iNPS Peak Caller
cd ${HOME}/lib
if [ ! -d "iNPS" ]; then
    mkdir iNPS
    cd iNPS
    wget http://www.picb.ac.cn/hanlab/files/iNPS_V1.2.2.zip
    unzip iNPS_V1.2.2.zip
fi

# bedTools
cd ${HOME}/lib
if [ ! -d "bedtools2" ]; then
    wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz
    tar -zxvf bedtools-2.27.1.tar.gz
    cd bedtools2
    make
fi

cd ${HOME}/lib
if [ ! -d "TrimGalore-0.4.3" ]; then
    wget -O trim_galore.tar.gz https://github.com/FelixKrueger/TrimGalore/archive/0.4.3.tar.gz
    tar -xzf trim_galore.tar.gz
fi

# Install R packages required by iDEAR
# cd ${HOME}/build/Multiscale-Genomics/mg-process-fastq
# chmod +x scripts/travis/install_packages.R
# sudo Rscript scripts/install_packages.R

# Post Installation Tidyup
cd ${HOME}/lib
rm -f *.zip *.tar.gz

# cd ${HOME}/build/Multiscale-Genomics/mg-process-fastq
