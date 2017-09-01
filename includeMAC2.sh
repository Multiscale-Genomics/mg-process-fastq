#!/bin/bash

python_version=$(python --version 2>&1)
echo $python_version
if [[ $python_version != *"3."* ]]; then
     cd ${HOME}/code
     pip install MACS2 
     ln -s ${HOME}/.pyenv/versions/mg-process-fastq/bin/macs2 ${HOME}/bin/macs2
     
     cd ${HOME}/lib
     wget https://github.com/3DGenomes/tadbit/archive/master.zip -O tadbit.zip
     unzip tadbit.zip
     cd TADbit-master
     yes | python setup.py install --install-lib=${HOME}/.pyenv/versions/mg-process-fastq/lib/python2.7/site-packages/ --install-scripts=${HOME}/bin 
     
else 
     cd ${HOME}/code
     pip install MACS2 == MACS2p3
     ln -s ${HOME}/.pyenv/versions/mg-process-fastq/bin/macs2p3 ${HOME}/bin/macs2
    
fi
