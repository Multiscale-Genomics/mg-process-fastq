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

pv=$(python -c 'import platform; v=platform.python_version_tuple(); print("{}.{}".format(v[0], v[1]))')
echo $pv
cd ${HOME}/lib
if [[ $pv == "2.7" ]]; then
    pip install MACS2
else
    if [ ! -d "MACS" ]; then
        git clone https://github.com/taoliu/MACS.git
        cd MACS
        git checkout MACS2p3

        cython MACS2/*.pyx
        cython MACS2/IO/*.pyx
        python setup_w_cython.py install

        pip install .
        alias macs2="macs2p3"
    fi
fi
