# mg-process-fastq

[![Documentation Status](https://readthedocs.org/projects/mg-process-fastq/badge/?version=latest)](http://mg-process-fastq.readthedocs.org/en/latest/) [![Build Status](https://travis-ci.org/Multiscale-Genomics/mg-process-fastq.svg?branch=master)](https://travis-ci.org/Multiscale-Genomics/mg-process-fastq) [![Code Health](https://landscape.io/github/Multiscale-Genomics/mg-process-fastq/master/landscape.svg?style=flat)](https://landscape.io/github/Multiscale-Genomics/mg-process-fastq/master)

Scripts required for the processing of FASTQ files (eg generating adjacency lists for Hi-C data)

# Requirements
- Python 2.7.12 (required for MACS2 in ChIP-Seq pipeline)
- Python 3.5.2 (required for iNPS in MNase-Seq pipeline)
- Python Modules:
  - numpy
  - h5py
  - scipy
  - matplotlib
  - TADbit
  - pysam
  - MACS2 - can be installed with pip, but runs on command line
  - rpy2
- GEMtools
- HDF5
- Bowtie2
- BWA
- SAMtools
- BS-Seeker2
- libmaus2
- BioBamBam2
- imp (for 3D modelling with TADbit)
- mcl
- R (2.9.1)
- iNPS


Installation
------------

For a guide to the full installation procedure the see [ReadTheDocs](http://mg-process-fastq.readthedocs.io).

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



