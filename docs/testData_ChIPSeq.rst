Test Data for ChIP-seq pipeline
===============================

The following document is for the preparation of data set required for testing the ChIP-seq pipeline. The document has
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

   BWA
   MACS 2
   Biobambam
   Samtools

Data set for genome file
------------------------

Filtering for required coverage
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Download the genome file from

.. code-block:: none

   wget http://www.ebi.ac.uk/ena/data/view/CM000663.2,CM000664.2,CM000665.2,CM000666.2,CM000667.2,CM000668.2,CM000669.2,CM000670.2,CM000671.2,CM000672.2,CM000673.2,CM000674.2,CM000675.2,CM000676.2,CM000677.2,CM000678.2,CM000679.2,CM000680.2,CM000681.2,CM000682.2,CM000683.2,CM000684.2,CM000685.2,CM000686.2,J01415.2&display=fasta&download=fasta&filename=entry.fasta -O GCA_000001405.22.fasta


Checkout https://github.com/Multiscale-Genomics/mg-misc-scripts/blob/master/ChIPSeq_Scripts/extract_chromosomeForChIP.py and extract chromosome 22 from the above file using the following command.

.. code-block:: none

   python extract_chromosomeForChIP.py path/to/your/input/file path/to/output/file

Download the fastq file from

.. code-block:: none

   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR000/DRR000150/DRR000150.fastq.gz
   
Unzip this file.

.. code-block:: none

   unzip DRR000150.fastq.gz


Index the fasta file

.. code-block:: none

   bwa index GCA_000001405.22.chr22.fa.fasta

Align the fastq file

.. code-block:: none

   bwa aln GCA_000001405.22.chr22.fa.fasta DRR000150.chr22.fastq >GCA_000001405.22.chr22.sai

And make the sam file

.. code-block:: none

   bwa samse GCA_000001405.22.chr22.fa.fasta GCA_000001405.22.chr22.sai DRR000150.chr22.fastq >GCA_000001405.22.chr22.sam

Sort the sam file

.. code-block:: none

   samtools sort GCA_000001405.22.chr22.sam >GCA_000001405.22.chr22.sorted.sam


Find the depths of coverage from the sorted file

.. code-block:: none

   samtools depth GCA_000001405.22.chr22.sorted.sam >GCA_000001405.22.chr22.dp


From the depth file, find regions with >= 70 depth, spanning over >=55 base pairs. You may get the script for this from https://github.com/Multiscale-Genomics/mg-misc-scripts/blob/master/ChIPSeq_Scripts/traverseForCoverageRegion_ChIP.py. Run it using:

.. code-block:: none

   python traverseForCoverageRegion_ChIP.py path/to/GCA_000001405.22.chr22.dp

Running this script would print the spanning regions. If it is a continuous region, you may only take the first starting base pair and the last ending base pair, as inputs for the next step. (Take out 1000 and add in 1000 to these respectively to get upstream and downstream spanning bases)

Extract the corresponding fasta sequence from the chromosome file for the positions retrieved from the above step. Checkout file from https://github.com/Multiscale-Genomics/mg-misc-scripts/blob/master/ChIPSeq_Scripts/extractChromosomalRegion.py and run using command:

.. code-block:: none

   python extractChromosomalRegion.py path/to/original/fasta/file path/to/output/file/for/region/macs2.Human.GCA_000001405.22.fasta starting_base_position ending_base_position

Index the fasta file for the selected region

.. code-block:: none

   bwa index macs2.Human.GCA_000001405.22.fasta

Align the fastq file

.. code-block:: none

   bwa aln macs2.Human.GCA_000001405.22.fasta DRR000150.chr22.fastq >macs2.Human.GCA_000001405.22.sai

And make the sam file

.. code-block:: none

   bwa samse macs2.Human.GCA_000001405.22.fasta macs2.Human.GCA_000001405.22.sai DRR000150.chr22.fastq >macs2.Human.GCA_000001405.22.sam

Filter this sam file for the reads which aligned with chromosome 22 using the following command:

.. code-block:: none

   awk '$3 == chr22' macs2.Human.GCA_000001405.22.sam >macs2.Human.GCA_000001405.22.22.sam

From the filtered reads from the above output file, extract the corresponding entries in fastq file. You may do this using the file at :

.. code-block:: none

   https://github.com/Multiscale-Genomics/mg-misc-scripts/blob/master/ChIPSeq_Scripts/makeFastQFiles.py

and running it via command line:

.. code-block:: none

   python makeFastQFiles.py --samfile macs2.Human.GCA_000001405.22.22.sam --fastQfile /path/to/DRR000150.chr22.fastq --pathToOutput /path/to/save/output/fastq/file/to/ --fastqOut macs2.Human.DRR000150.22.fastq

The fastq file in the above step and fasta file macs2.Human.GCA_000001405.22.fasta together make up the data set for ChIP-seq pipeline
