#!/bin/bash

python_version=$(python --version 2>&1)
echo $python_version
if [[ $python_version != *"3."* ]]; then
     cd ${HOME}/code
     pip install MACS2    
fi
