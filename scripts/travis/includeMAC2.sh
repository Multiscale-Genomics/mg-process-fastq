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

python_version=$(python --version 2>&1)
echo $python_version
if [[ $python_version != *"3."* ]]; then
    cd ${HOME}/code
    pip install MACS2
    ln -s ${HOME}/.pyenv/versions/mg-process-fastq/bin/macs2 ${HOME}/bin/macs2
else
    cd ${HOME}/code
    pip install MACS2 == MACS2p3
    ln -s ${HOME}/.pyenv/versions/mg-process-fastq/bin/macs2p3 ${HOME}/bin/macs2

fi
