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
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip fastqc_v0.11.5.zip
cd FastQC/
chmod 755 fastqc

# htslib
cd ${HOME}/lib
git clone https://github.com/samtools/htslib.git
cd htslib
autoheader
autoconf
./configure --prefix=${HOME}/lib/htslib
make
make install

# SAMtools
cd ${HOME}/lib
git clone https://github.com/samtools/samtools.git
cd samtools
autoheader
autoconf -Wno-syntax
./configure --prefix=${HOME}/lib/samtools
make
make install

# UCSC Tools
cd ${HOME}/lib
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig

# Bowtie2 Aligner
cd ${HOME}/lib
wget --max-redirect 1 https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4/bowtie2-2.3.4-linux-x86_64.zip
unzip bowtie2-2.3.4-linux-x86_64.zip

# BWA Sequence Aligner
cd ${HOME}/lib
git clone https://github.com/lh3/bwa.git
cd bwa
make

# GEM Indexer
cd ${HOME}/lib
wget http://barnaserver.com/gemtools/releases/GEMTools-static-core2-1.7.1.tar.gz
tar -xzf GEMTools-static-core2-1.7.1.tar.gz

# iNPS Peak Caller
cd ${HOME}/lib
mkdir iNPS
cd iNPS
wget http://www.picb.ac.cn/hanlab/files/iNPS_V1.2.2.zip
unzip iNPS_V1.2.2.zip

# Kallisto
cd ${HOME}/lib
wget https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz
tar -xzf kallisto_linux-v0.43.1.tar.gz

# bedTools
wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz
tar -zxvf bedtools-2.26.0.tar.gz
cd bedtools2
make
# cd ${HOME}/build/Multiscale-Genomics/mg-process-fastq

# Install MACS2
cd ${HOME}/build/Multiscale-Genomics/mg-process-fastq
chmod +x scripts/travis/includeMAC2.sh
./scripts/travis/includeMAC2.sh

# Install R packages required by iDEAR
# cd ${HOME}/build/Multiscale-Genomics/mg-process-fastq
# chmod +x scripts/travis/includeMAC2.sh
# sudo Rscript scripts/install_packages.R

# Install TADbit (only on Python 2.7)
cd ${HOME}/build/Multiscale-Genomics/mg-process-fastq
chmod +x scripts/travis/includeTADbit.sh
./scripts/travis/includeTADbit.sh

# Post Installation Tidyup
cd ${HOME}/lib
rm *.zip *.tar.gz