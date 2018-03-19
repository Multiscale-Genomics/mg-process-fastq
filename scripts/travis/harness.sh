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
pv=$(python -c 'import platform; print(platform.python_version())')

python tests/test_toolchains.py --pipeline genome
tc=$?
rc=$(($rc + $tc))
bash tidy_data.sh

python tests/test_pipelines.py --pipeline genome
tc=$?
rc=$(($rc + $tc))
bash tidy_data.sh

# # Not included for the moment due to installation issues with versions or R and python
# python tests/test_toolchains.py --pipeline idamidseq
# tc=$?
# rc=$(($rc + $tc))
# sh tidy_data.sh

python tests/test_toolchains.py --pipeline bowtie2
tc=$?
rc=$(($rc + $tc))
bash tidy_data.sh

python tests/test_toolchains.py --pipeline bwa
tc=$?
rc=$(($rc + $tc))
bash tidy_data.sh

if [[ $pv == "2.7.12" ]]; then
    python tests/test_toolchains.py --pipeline chipseq
    tc=$?
    rc=$(($rc + $tc))
    bash tidy_data.sh

    python tests/test_pipelines.py --pipeline chipseq
    tc=$?
    rc=$(($rc + $tc))
    bash tidy_data.sh
fi

python tests/test_toolchains.py --pipeline rnaseq
tc=$?
rc=$(($rc + $tc))
bash tidy_data.sh

python tests/test_pipelines.py --pipeline rnaseq
tc=$?
rc=$(($rc + $tc))
bash tidy_data.sh

if [[ $pv == "2.7.12" ]]; then
    python tests/test_toolchains.py --pipeline wgbs
    tc=$?
    rc=$(($rc + $tc))
    bash tidy_data.sh

    python tests/test_pipelines.py --pipeline wgbs
    tc=$?
    rc=$(($rc + $tc))
    bash tidy_data.sh
fi

# if [[ $python_version == *"3."* ]]; then
#     python tests/test_toolchains.py --pipeline mnaseseq
#     tc=$?
#     rc=$(($rc + $tc))
#     bash tidy_data.sh
#
#     python tests/test_pipelines.py --pipeline mnaseseq
#     tc=$?
#     rc=$(($rc + $tc))
#     bash tidy_data.sh
# fi

if [[ $rc != 0 ]]; then exit $rc; fi
