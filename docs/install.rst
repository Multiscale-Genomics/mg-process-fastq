Requirements and Installation
=============================

Requirements
============

Software
--------
- bedtools
- bedToBigBed - http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
- wigToBigWig - http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
- biobambam2
- Bowtie2
- BWA
- GEM
- HDF5
- iNPS
- Kallisto
- Python 2.7.10+
- pyenv

Pytho Modules
-------------

- numpy
- h5py
- pysam

Installation
============

Directly from GitHub:

.. code-block:: none
   :linenos:
   
   git clone https://github.com/Multiscale-Genomics/mg-process-fastq.git
   
   cd mg-process-fastq

Create the Python environment

.. code-block:: none
   pyenv-virtualenv 2.7.10 mg-process-fastq
   pip install --editable .
   
   
