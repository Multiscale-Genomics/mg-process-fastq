Requirements and Installation
=============================

`Source (github) <https://github.com/Multiscale-Genomics/mg-process-fastq>`_

Requirements
------------

Software
^^^^^^^^

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
- Python 2.7.10+
- pyenv
- R 2.9.1+
- SAMtools
- MCL

Python Modules
^^^^^^^^^^^^^^

- numpy
- h5py
- pysam
- scipy
- matplotlib
- rpy2
- BS-Seeker2
- TADbit

Installation
------------

Directly from GitHub:

.. code-block:: none
   :linenos:
   
   code_root = pwd

   git clone https://github.com/Multiscale-Genomics/mg-process-fastq.git
   
   cd mg-process-fastq

Create the Python environment

.. code-block:: none
   :linenos:
   
   pyenv-virtualenv 2.7.10 mg-process-fastq
   pip install --editable .

Check out the following software for use by the process_wgbs.py pipeline:

.. code-block:: none
	:linenos:

	cd $code_root
	gti clone https://github.com/BSSeeker/BSseeker2.git

	cd mg-process-fastq
	ln -s $code_root/bs_align bs_align
	ln -s $code_root/bs_index bs_index
	ln -s $code_root/bs_utils bs_utils
	
	cd tool
	ln -s $code_root/FilterReads.py FilterReads.py
   

