Full Installation
=================

The following document is for the full installation of all software required by
the mg-process-fastq module and all programmes that it uses. The document has
been written with Ubuntu Linux in mind, although many of the commands are cross
platform (\*nix) complient.

If you already have certain packages installed feel free to skip over certain
steps. Likewise the bin, lib and code directories are relative to the home dir;
if this is not the case for your system then make the required changes when
running these commands.

Setup the System Environment
----------------------------

.. code-block:: none
   :linenos:

   sudo apt-get install -y make build-essential libssl-dev zlib1g-dev       //
   libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev //
   libncursesw5-dev xz-utils tk-dev unzip mcl libgtk2.0-dev r-base-core     //
   libcurl4-gnutls-dev python-rpy2 git

   cd ${HOME}
   mkdir bin lib code
   echo 'export PATH="${HOME}/bin:$PATH"' >> ~/.bash_profile

Setup pyenv and pyenv-virtualenv
--------------------------------

This is required for managing the version of Python and the installation
environment for the Python modules so that they can be installed in the user
space.

.. code-block:: none
   :linenos:

   git clone https://github.com/pyenv/pyenv.git ~/.pyenv
   echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bash_profile
   echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bash_profile
   echo 'eval "$(pyenv init -)"' >> ~/.bash_profile

   git clone https://github.com/pyenv/pyenv-virtualenv.git ${PYENV_ROOT}/plugins/pyenv-virtualenv

   pyenv install 2.7.12
   pyenv virtualenv 2.7.12 mg-process-fastq

Installation Process
--------------------

UCSC Tools
^^^^^^^^^^

.. code-block:: none
   :linenos:

   cd ${HOME}/lib
   wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
   wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig

BioBamBam2
^^^^^^^^^^

BioBamBam is used for the filtering of aligned reads as part of the ChIP-seq
pipeline. It also requires the libmaus2 package to be installed.

.. code-block:: none
   :linenos:

   cd ${HOME}/lib
   git clone https://github.com/gt1/libmaus2.git
   cd libmaus2
   libtoolize
   aclocal
   autoheader
   automake --force-missing --add-missing
   autoconf
   ./configure --prefix=${HOME}/lib/libmaus2
   make
   make install

   cd ${HOME}/lib
   git clone https://github.com/gt1/biobambam2.git
   cd biobambam2
   autoreconf -i -f
   ./configure --with-libmaus2=${HOME}/lib/libmaus2 --prefix=${HOME}/lib/biobambam2
   make install

Bowtie2 Aligner
^^^^^^^^^^^^^^^

.. code-block:: none
   :linenos:

   cd ${HOME}/lib
   wget --max-redirect 1 https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.2/bowtie2-2.3.2-linux-x86_64.zip
   unzip bowtie2-2.3.2-linux-x86_64.zip

BWA Sequence Aligner
^^^^^^^^^^^^^^^^^^^^

.. code-block:: none
   :linenos:

   cd ${HOME}/lib
   git clone https://github.com/lh3/bwa.git
   cd bwa
   make

GEM Sequence Aligner
^^^^^^^^^^^^^^^^^^^^

.. code-block:: none
   :linenos:

   cd ${HOME}/lib
   wget http://barnaserver.com/gemtools/releases/GEMTools-static-core2-1.7.1.tar.gz
   tar -xzf GEMTools-static-core2-1.7.1.tar.gz

iNPS Peak Caller
^^^^^^^^^^^^^^^^

.. code-block:: none
   :linenos:

   cd ${HOME}/lib
   mkdir iNPS
   cd iNPS
   wget http://www.picb.ac.cn/hanlab/files/iNPS_V1.2.2.zip
   unzip iNPS_V1.2.2.zip

Kallisto
^^^^^^^^

.. code-block:: none
   :linenos:

   cd ${HOME}/lib
   wget https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz
   tar -xzf kallisto_linux-v0.43.1.tar.gz

SAMtools
^^^^^^^^

.. code-block:: none
   :linenos:

   cd ${HOME}/lib
   git clone https://github.com/samtools/htslib.git
   cd htslib
   autoheader
   autoconf
   ./configure --prefix=${HOME}/lib/htslib
   make
   make install

   cd ${HOME}/lib
   git clone https://github.com/samtools/samtools.git
   cd samtools
   autoheader
   autoconf -Wno-syntax
   ./configure --prefix=${HOME}/lib/samtools
   make
   make install

bedTools
^^^^^^^^

.. code-block:: none
   :linenos:

   wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz
   tar -zxvf bedtools-2.26.0.tar.gz
   cd bedtools2
   make


Setup the symlinks
------------------

.. code-block:: none
   :linenos:

   cd ${HOME}/bin

   ln -s ${HOME}/lib/bedtools2/bin/bedtools bedtools

   ln -s ${HOME}/lib/bedToBigBed bedToBigBed
   ln -s ${HOME}/lib/wigToBigWig wigToBigWig

   ln -s ${HOME}/lib/bwa/bwa bwa

   ln -s ${HOME}/lib/bowtie2-2.3.2/bowtie2 bowtie2
   ln -s ${HOME}/lib/bowtie2-2.3.2/bowtie2-align-s bowtie2-align-s
   ln -s ${HOME}/lib/bowtie2-2.3.2/bowtie2-align-l bowtie2-align-l
   ln -s ${HOME}/lib/bowtie2-2.3.2/bowtie2-build bowtie2-build
   ln -s ${HOME}/lib/bowtie2-2.3.2/bowtie2-build-s bowtie2-build-s
   ln -s ${HOME}/lib/bowtie2-2.3.2/bowtie2-build-l bowtie2-build-l
   ln -s ${HOME}/lib/bowtie2-2.3.2/bowtie2-inspect bowtie2-inspect
   ln -s ${HOME}/lib/bowtie2-2.3.2/bowtie2-inspect-s bowtie2-inspect-s
   ln -s ${HOME}/lib/bowtie2-2.3.2/bowtie2-inspect-l bowtie2-inspect-l

   ln -s ${HOME}/lib/gemtools-1.7.1-core2/bin/gem-2-bed gem-2-bed
   ln -s ${HOME}/lib/gemtools-1.7.1-core2/bin/gem-2-gem gem-2-gem
   ln -s ${HOME}/lib/gemtools-1.7.1-core2/bin/gem-2-sam gem-2-sam
   ln -s ${HOME}/lib/gemtools-1.7.1-core2/bin/gem-2-wig gem-2-wig
   ln -s ${HOME}/lib/gemtools-1.7.1-core2/bin/gem-indexer gem-indexer
   ln -s ${HOME}/lib/gemtools-1.7.1-core2/bin/gem-indexer_bwt-dna gem-indexer_bwt-dna
   ln -s ${HOME}/lib/gemtools-1.7.1-core2/bin/gem-indexer_fasta2meta+cont gem-indexer_fasta2meta+cont
   ln -s ${HOME}/lib/gemtools-1.7.1-core2/bin/gem-indexer_generate gem-indexer_generate
   ln -s ${HOME}/lib/gemtools-1.7.1-core2/bin/gem-info gem-info
   ln -s ${HOME}/lib/gemtools-1.7.1-core2/bin/gem-mapper gem-mapper
   ln -s ${HOME}/lib/gemtools-1.7.1-core2/bin/gemtools gemtools

   ln -s ${HOME}/lib/iNPS/iNPS_V1.2.2.py iNPS_V1.2.2.py
   ln -s ${HOME}/lib/kallisto_linux-v0.43.1/kallisto kallisto

   ln -s ${HOME}/lib/htslib/bin/bgzip bgzip
   ln -s ${HOME}/lib/htslib/bin/htsfile htsfile
   ln -s ${HOME}/lib/htslib/bin/tabix tabix

   ln -s ${HOME}/lib/samtools/bin/ace2sam ace2sam
   ln -s ${HOME}/lib/samtools/bin/blast2sam.pl blast2sam.pl
   ln -s ${HOME}/lib/samtools/bin/bowtie2sam.pl bowtie2sam.pl
   ln -s ${HOME}/lib/samtools/bin/export2sam.pl export2sam.pl
   ln -s ${HOME}/lib/samtools/bin/interpolate_sam.pl interpolate_sam.pl
   ln -s ${HOME}/lib/samtools/bin/maq2sam-long maq2sam-long
   ln -s ${HOME}/lib/samtools/bin/maq2sam-short maq2sam-short
   ln -s ${HOME}/lib/samtools/bin/md5fa md5fa
   ln -s ${HOME}/lib/samtools/bin/md5sum-lite md5sum-lite
   ln -s ${HOME}/lib/samtools/bin/novo2sam.pl novo2sam.pl
   ln -s ${HOME}/lib/samtools/bin/plot-bamstats plot-bamstats
   ln -s ${HOME}/lib/samtools/bin/psl2sam.pl psl2sam.pl
   ln -s ${HOME}/lib/samtools/bin/sam2vcf.pl sam2vcf.pl
   ln -s ${HOME}/lib/samtools/bin/samtools samtools
   ln -s ${HOME}/lib/samtools/bin/samtools.pl samtools.pl
   ln -s ${HOME}/lib/samtools/bin/seq_cache_populate.pl seq_cache_populate.pl
   ln -s ${HOME}/lib/samtools/bin/soap2sam.pl soap2sam.pl
   ln -s ${HOME}/lib/samtools/bin/varfilter.py varfilter.py
   ln -s ${HOME}/lib/samtools/bin/wgsim wgsim
   ln -s ${HOME}/lib/samtools/bin/wgsim_eval.pl wgsim_eval.pl
   ln -s ${HOME}/lib/samtools/bin/zoom2sam.pl zoom2sam.pl

   ln -s ${HOME}/lib/biobambam2/bin/bam12auxmerge bam12auxmerge
   ln -s ${HOME}/lib/biobambam2/bin/bam12split bam12split
   ln -s ${HOME}/lib/biobambam2/bin/bam12strip bam12strip
   ln -s ${HOME}/lib/biobambam2/bin/bamadapterclip bamadapterclip
   ln -s ${HOME}/lib/biobambam2/bin/bamadapterfind bamadapterfind
   ln -s ${HOME}/lib/biobambam2/bin/bamalignfrac bamalignfrac
   ln -s ${HOME}/lib/biobambam2/bin/bamauxmerge bamauxmerge
   ln -s ${HOME}/lib/biobambam2/bin/bamauxsort bamauxsort
   ln -s ${HOME}/lib/biobambam2/bin/bamcat bamcat
   ln -s ${HOME}/lib/biobambam2/bin/bamchecksort bamchecksort
   ln -s ${HOME}/lib/biobambam2/bin/bamclipreinsert bamclipreinsert
   ln -s ${HOME}/lib/biobambam2/bin/bamcollate bamcollate
   ln -s ${HOME}/lib/biobambam2/bin/bamcollate2 bamcollate2
   ln -s ${HOME}/lib/biobambam2/bin/bamdownsamplerandom bamdownsamplerandom
   ln -s ${HOME}/lib/biobambam2/bin/bamexplode bamexplode
   ln -s ${HOME}/lib/biobambam2/bin/bamfilteraux bamfilteraux
   ln -s ${HOME}/lib/biobambam2/bin/bamfilterflags bamfilterflags
   ln -s ${HOME}/lib/biobambam2/bin/bamfilterheader bamfilterheader
   ln -s ${HOME}/lib/biobambam2/bin/bamfilterheader2 bamfilterheader2
   ln -s ${HOME}/lib/biobambam2/bin/bamfilterlength bamfilterlength
   ln -s ${HOME}/lib/biobambam2/bin/bamfiltermc bamfiltermc
   ln -s ${HOME}/lib/biobambam2/bin/bamfilternames bamfilternames
   ln -s ${HOME}/lib/biobambam2/bin/bamfilterrg bamfilterrg
   ln -s ${HOME}/lib/biobambam2/bin/bamfixmateinformation bamfixmateinformation
   ln -s ${HOME}/lib/biobambam2/bin/bamflagsplit bamflagsplit
   ln -s ${HOME}/lib/biobambam2/bin/bamheap2 bamheap2
   ln -s ${HOME}/lib/biobambam2/bin/bamindex bamindex
   ln -s ${HOME}/lib/biobambam2/bin/bamintervalcomment bamintervalcomment
   ln -s ${HOME}/lib/biobambam2/bin/bamintervalcommenthist bamintervalcommenthist
   ln -s ${HOME}/lib/biobambam2/bin/bamlastfilter bamlastfilter
   ln -s ${HOME}/lib/biobambam2/bin/bammapdist bammapdist
   ln -s ${HOME}/lib/biobambam2/bin/bammarkduplicates bammarkduplicates
   ln -s ${HOME}/lib/biobambam2/bin/bammarkduplicates2 bammarkduplicates2
   ln -s ${HOME}/lib/biobambam2/bin/bammarkduplicatesopt bammarkduplicatesopt
   ln -s ${HOME}/lib/biobambam2/bin/bammaskflags bammaskflags
   ln -s ${HOME}/lib/biobambam2/bin/bammdnm bammdnm
   ln -s ${HOME}/lib/biobambam2/bin/bammerge bammerge
   ln -s ${HOME}/lib/biobambam2/bin/bamnumericalindex bamnumericalindex
   ln -s ${HOME}/lib/biobambam2/bin/bamrank bamrank
   ln -s ${HOME}/lib/biobambam2/bin/bamranksort bamranksort
   ln -s ${HOME}/lib/biobambam2/bin/bamrecalculatecigar bamrecalculatecigar
   ln -s ${HOME}/lib/biobambam2/bin/bamrecompress bamrecompress
   ln -s ${HOME}/lib/biobambam2/bin/bamreset bamreset
   ln -s ${HOME}/lib/biobambam2/bin/bamscrapcount bamscrapcount
   ln -s ${HOME}/lib/biobambam2/bin/bamseqchksum bamseqchksum
   ln -s ${HOME}/lib/biobambam2/bin/bamsormadup bamsormadup
   ln -s ${HOME}/lib/biobambam2/bin/bamsort bamsort
   ln -s ${HOME}/lib/biobambam2/bin/bamsplit bamsplit
   ln -s ${HOME}/lib/biobambam2/bin/bamsplitdiv bamsplitdiv
   ln -s ${HOME}/lib/biobambam2/bin/bamstreamingmarkduplicates bamstreamingmarkduplicates
   ln -s ${HOME}/lib/biobambam2/bin/bamtagconversion bamtagconversion
   ln -s ${HOME}/lib/biobambam2/bin/bamtofastq bamtofastq
   ln -s ${HOME}/lib/biobambam2/bin/bamvalidate bamvalidate
   ln -s ${HOME}/lib/biobambam2/bin/bamzztoname bamzztoname
   ln -s ${HOME}/lib/biobambam2/bin/fastaexplode fastaexplode
   ln -s ${HOME}/lib/biobambam2/bin/fastqtobam fastqtobam
   ln -s ${HOME}/lib/biobambam2/bin/fastqtobampar fastqtobampar
   ln -s ${HOME}/lib/biobambam2/bin/filtersam filtersam
   ln -s ${HOME}/lib/biobambam2/bin/kmerprob kmerprob
   ln -s ${HOME}/lib/biobambam2/bin/lasToBAM lasToBAM
   ln -s ${HOME}/lib/biobambam2/bin/normalisefasta normalisefasta


Prepare the Python Environment
------------------------------

Install APIs and Pipelines
^^^^^^^^^^^^^^^^^^^^^^^^^^

Checkout the code for the DM API and the mg-process-fastq pipelines:

.. code-block:: none
   :linenos:

   cd ${HOME}/code
   pyenv activate mg-process-fastq
   pip install git+https://github.com/Multiscale-Genomics/mg-dm-api.git
   pip install git+https://github.com/Multiscale-Genomics/mg-tool-api.git

   git clone https://github.com/Multiscale-Genomics/mg-process-fastq.git
   cd mg-process-fastq
   pip install --editable .


Install MACS2
^^^^^^^^^^^^^

This should get installed as part of the installation in the mg-process-fastq
package, if not then it will need to be installed separately:

.. code-block:: none
   :linenos:

   cd ${HOME}/code
   pyenv activate mg-process-fastq
   pip install MACS2


Install TADbit
^^^^^^^^^^^^^^

.. code-block:: none
   :linenos:

   cd ${HOME}/lib
   wget https://github.com/3DGenomes/tadbit/archive/master.zip -O tadbit.zip
   unzip tadbit.zip
   cd TADbit-master

   # If the pyenv env is not called mg-process-fastq then change this to match,
   # the sme is true for teh version of python
   python setup.py install --install-lib=${HOME}/.pyenv/versions/mg-process-fastq/lib/python2.7/site-packages/ --install-scripts=${HOME}/bin

Install BSseeker
^^^^^^^^^^^^^^^^

.. code-block:: none
   :linenos:

   cd ${HOME}/lib
   git clone https://github.com/BSSeeker/BSseeker2.git

   cd ${HOME}/code/mg-process-fastq
   ln -s ${HOME}/lib/BSseeker2/bs_align bs_align
   ln -s ${HOME}/lib/BSseeker2/bs_index bs_index
   ln -s ${HOME}/lib/BSseeker2/bs_utils bs_utils

   cd ${HOME}/code/mg-process-fastq/tool
   ln -s ${HOME}/lib/BSseeker2/FilterReads.py FilterReads.py

Post Installation Tidyup
------------------------

.. code-block:: none
   :linenos:

   cd ${HOME}/lib
   rm *.zip *.tar.gz