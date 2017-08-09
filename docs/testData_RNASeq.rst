Test Data for RNA-seq pipeline
===============================

The following document is for the preparation of data set required for testing the RNA-seq pipeline. The document has
been written with macOS Sierra in mind, although many of the commands are cross
platform (\*nix) complient.

You would need to have the tools listed in "Prerequisites" installed on your system.
For more details on installing the tools for this pipeline please refer to

http://multiscale-genomics.readthedocs.io/projects/mg-process-fastq/en/latest/full_installation.html

If you already have certain packages installed feel free to skip over certain
steps. Likewise the bin, lib and code directories are relative to the home dir;
if this is not the case for your system then make the required changes when
running these commands.

Prerequisites
-------------

.. code-block:: none
   :linenos:

   Kallisto
   Samtools


Data set for genome file
------------------------

Go to Ensemble website >> Human >> Example gene
 
.. code-block:: none

   http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000139618;r=13:32315474-32400266
   
Copy the chromosome number and coordinates given in the "location" field. Go to BioMart (top panel), and select Filters from the left panel. Expand Regions and enter the information retrieved above. 

Click on "Attributes" in the left panel. Select Gene stable ID, Transcript stable ID from Features. Select cDNA sequences from Sequences radio button.

Click on the Results button above the left panel. Export results to fasta file. 

Index this file using Kallisto indexer : 

.. code-block:: none

   kallisto index -i kallisto.Human.GRCh38.fasta.idx /path/to/file/exportSequences.fasta
   
Download the fastq files 

.. code-block:: none

   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR030/ERR030872/ERR030872_1.fastq.gz
   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR030/ERR030872/ERR030872_2.fastq.gz
   
Run the Kallisto quantifier using command : 

.. code-block:: none

   kallisto quant -i kallisto.Human.GRCh38.fasta.idx -o out --pseudobam /path/to/ERR030872_1.fastq.gz /path/to/ERR030872_2.fastq.gz  >kallisto.ERR030872.sam
   
Filter the aligned sequence entries from the above sam file : 

.. code-block:: none

   awk '$3 != "*"' kallisto.ERR030872.sam >kallisto.ERR030872.filtered.sam
   
Unzip the fastq files.

.. code-block:: none

   unzip ERR030872_1.fastq.gz
   unzip ERR030872_2.fastq.gz

   
Checkout https://github.com/Multiscale-Genomics/mg-misc-scripts/blob/master/RNASeq_Scripts/makeFastQFiles.py  and use the following command to generate the fastq files : 

.. code-block:: none

   python /path/to/makeFastQFiles.py --samfile kallisto.ERR030872.filtered.sam --fastQfile ERR030872_1.fastq --pathToOutput /path/to/make/fastqFile/ --fastqOut ERR030872_1.RNAseq.fastq
   python /path/to/makeFastQFiles.py --samfile kallisto.ERR030872.filtered.sam --fastQfile ERR030872_2.fastq --pathToOutput /path/to/make/fastqFile/ --fastqOut ERR030872_2.RNAseq.fastq
   
   
Shorten these files by running the script at https://github.com/Multiscale-Genomics/mg-misc-scripts/blob/master/RNASeq_Scripts/randomSeqSelector.py

using 

.. code-block:: none

   python PythonScripts/randomSeqSelector.py ERR030872_1.RNAseq.fastq kallisto.Human.ERR030872_1.fastq
   python PythonScripts/randomSeqSelector.py ERR030872_2.RNAseq.fastq kallisto.Human.ERR030872_2.fastq

Then zip them : 

.. code-block:: none

   gzip kallisto.Human.ERR030872_1.fastq
   gzip kallisto.Human.ERR030872_2.fastq
   
   
