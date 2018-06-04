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

# FastQC
cd ${HOME}/lib
if [ ! -d "FastQC" ]; then
    wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
    unzip fastqc_v0.11.5.zip
    cd FastQC/
    chmod 755 fastqc
fi

# Kallisto
cd ${HOME}/lib
if [ ! -d "kallisto_linux-v0.43.1" ]; then
    wget https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz
    tar -xzf kallisto_linux-v0.43.1.tar.gz
fi
