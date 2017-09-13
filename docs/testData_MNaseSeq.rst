Test Data for MNase-seq pipeline
===============================

The following document is for the preparation of data set required for testing the MNase-seq pipeline. The document has
been written with macOS Sierra in mind, although many of the commands are cross
platform (\*nix) compliant.

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

   BWA   
   Samtools
   iNPS
   

.. note:: iNPS will only run within a Python 3 or above environment

Data set for genome file
------------------------

Filtering for required coverage
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Download the fasta file for Mouse chromosome 19 from

.. code-block:: none

   wget "http://www.ebi.ac.uk/ena/data/view/CM001012&display=fasta&download=fasta&filename=CM001012.fasta" -O Mouse.CM001012.2.fasta

Download the fastq file from

.. code-block:: none

   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR000/DRR000386/DRR000386.fastq.gz
   
Unzip this file.

.. code-block:: none

   gunzip DRR000386.fastq.gz


Index the fasta file

.. code-block:: none

   bwa index Mouse.CM001012.2.fasta

Align the fastq file

.. code-block:: none

   bwa aln Mouse.CM001012.2.fasta DRR000386.fastq >Mouse.CM001012.2.sai

And make the sam file

.. code-block:: none

   bwa samse Mouse.CM001012.2.fasta Mouse.CM001012.2.sai DRR000386.fastq >Mouse.CM001012.2.sam

Filter out the aligned reads from the above sam file.

.. code-block:: none

   awk '$3 != "*"' Mouse.CM001012.2.sam >Mouse.CM001012.2.filtered.sam

Sort the sam file

.. code-block:: none

   samtools sort Mouse.CM001012.2.filtered.sam >Mouse.CM001012.2.sorted.sam

Find the depths of coverage from the sorted file

.. code-block:: none

   samtools depth Mouse.CM001012.2.sorted.sam >Mouse.CM001012.2.dp


From the depth file, find regions with >= 70 depth, spanning over >=55 base pairs. You may get the script for this from https://github.com/Multiscale-Genomics/mg-misc-scripts/blob/master/MNaseSeq_Scripts/traverseForCoverageRegion_MNase.py. Run it using:

.. code-block:: none

   python traverseForCoverageRegion_MNase.py path/to/Mouse.CM001012.2.dp

Running this script would print the spanning regions. Running this script for this data set gives multiple regions. The output is in the format : start - end - depth.  The one at the end has a maximal coverage from this data set. Since it is a continuous region, you may take the first starting base pair and the last ending base pair, as inputs for the next step. (Take out 1000 and add in 1000 to these respectively to get upstream and downstream spanning bases)

Extract the corresponding fasta sequence from the chromosome file for the positions retrieved from the above step. Checkout file from https://github.com/Multiscale-Genomics/mg-misc-scripts/blob/master/MNaseSeq_Scripts/extractChromosomalRegion.py and run using command:

.. code-block:: none

   python extractChromosomalRegion.py path/to/original/fasta/file path/to/output/file/for/region/inps.Mouse.GRCm38.fasta starting_base_position ending_base_position

Making the Fastq file
^^^^^^^^^^^^^^^^^^^^^^

Index the fasta file for the selected region

.. code-block:: none

   bwa index inps.Mouse.GRCm38.fasta

Align the fastq file

.. code-block:: none

   bwa aln inps.Mouse.GRCm38.fasta DRR000386.fastq >inps.Mouse.GRCm38.sai

And make the sam file

.. code-block:: none

   bwa samse inps.Mouse.GRCm38.fasta inps.Mouse.GRCm38.sai DRR000386.fastq >inps.Mouse.GRCm38.sam

Filter this sam file for the reads which aligned with chromosome 19 using the following command:

.. code-block:: none

   awk '$3 != "*"' inps.Mouse.GRCm38.sam >inps.Mouse.GRCm38.sam.19.sam

From the filtered reads from the above output file, extract the corresponding entries in fastq file. You may do this using the file at :

.. code-block:: none

   https://github.com/Multiscale-Genomics/mg-misc-scripts/blob/master/MNaseSeq_Scripts/makeFastQFiles.py

and running it via command line:

.. code-block:: none

   python makeFastQFiles.py --samfile path/to/inps.Mouse.GRCm38.sam.19.sam --fastQfile /path/to/DRR000386.fastq --pathToOutput /path/to/save/output/fastq/file/to/ --fastqOut DRR000386.MNaseseq.fastq

Shorten this file by running the script at https://github.com/Multiscale-Genomics/mg-misc-scripts/blob/master/MNASeq_Scripts/randomSeqSelector.py

using 

.. code-block:: none

   python randomSeqSelector.py DRR000386.MNaseseq.fastq inps.Mouse.DRR000386.fastq

The fastq file in the above step and fasta file inps.Mouse.GRCm38.fasta together make up the data set for MNase-seq pipeline
