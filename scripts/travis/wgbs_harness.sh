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
pv=$(python -c 'import platform; v=platform.python_version_tuple(); print("{}.{}".format(v[0], v[1]))')

if [[ $pv == "2.7" ]]; then
    if [[ $TESTENV == "wgbs_code_1" ]]; then
        python mg_process_fastq/tests/test_toolchains.py --pipeline wgbs
        tc=$?
        rc=$(($rc + $tc))
        bash tidy_data.sh
    else
        python mg_process_fastq/tests/test_pipelines.py --pipeline wgbs
        tc=$?
        rc=$(($rc + $tc))
        bash tidy_data.sh
    fi
fi

if [[ $rc != 0 ]]; then exit $rc; fi
