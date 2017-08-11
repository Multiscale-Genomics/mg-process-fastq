#!/bin/bash

python_version=$(python --version 2>&1)
echo $python_version
if [[ $python_version != *"3."* ]]; then
     cd ${HOME}/code
     pip install MACS2 
     ln -s ${HOME}/.pyenv/versions/mg-process-fastq/bin/macs2 ${HOME}/bin/macs2
fi
