.. Copyright 2017 EMBL-European Bioinformatics Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Requirements and Installation
=============================

`Source (github) <https://github.com/Multiscale-Genomics/mg-process-fastq>`_

Requirements
------------

Software
^^^^^^^^
- Python 2.7.10+
- bedtools (specifically bamtobed)
- bedToBigBed - http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
- wigToBigWig - http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
- BioBamBam2
- Bowtie2
- BWA
- GEMtools
- HDF5
- iNPS
- Kallisto
- libmaus2
- pyenv
- R 2.9.1+
- SAMtools
- MCL
- pigz

Python Modules
^^^^^^^^^^^^^^

- numpy
- h5py
- pysam
- pyBigWig
- scipy
- matplotlib
- rpy2

Instructions for the follownig modules are listed in the installation section.
All other python modules should be installed with pip prior to the following
libraries.

- BS-Seeker2
- TADbit

Installation
------------

For a guide to the full installation procedure the :doc:`full_installation`.

Directly from GitHub:

.. code-block:: none
   :linenos:

   cd ${HOME}/code

   git clone https://github.com/Multiscale-Genomics/mg-process-fastq.git

   cd mg-process-fastq

Create the Python environment

.. code-block:: none
   :linenos:

   pyenv-virtualenv 2.7.10 mg-process-fastq
   pip install --editable .

Install the pyTADbit modules

.. code-block:: none
   :linenos:

   cd ${HOME}/lib
   wget https://github.com/3DGenomes/tadbit/archive/master.zip -O tadbit.zip
   unzip tadbit.zip
   cd tadbit-master

   pyenv activate mg-process-fastq
   pip install .

Check out the following software for use by the process_wgbs.py pipeline:

.. code-block:: none
   :linenos:

   cd cd ${HOME}/lib
   gti clone https://github.com/BSSeeker/BSseeker2.git

   cd ${HOME}/code
   cd mg-process-fastq
   ln -s $code_root/bs_align bs_align
   ln -s $code_root/bs_index bs_index
   ln -s $code_root/bs_utils bs_utils

   cd cd ${HOME}/code/mg-process-fastq/tool
   ln -s $code_root/FilterReads.py FilterReads.py


Documentation
-------------
To build the documentation:

.. code-block:: none
   :linenos:

   pip install Sphinx
   pip install sphinx-autobuild
   cd docs
   make html
