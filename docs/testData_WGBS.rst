Test Data for WGBS pipeline
===============================

The following document is for the preparation of data set required for testing the WGBS pipeline. The document has
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

   BS Seeker
   Samtools
   Bowtie indexer


Data set for genome file
------------------------

Filtering for required coverage
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Download the genome file from

.. code-block:: none

     wget "http://www.ebi.ac.uk/ena/data/view/CM001012&display=fasta&download=fasta&filename=CM001012.fasta" -O Mouse.CM001012.2.fasta

Download the fastq files from

.. code-block:: none

   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR892/SRR892982/SRR892982_1.fastq.gz
   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR892/SRR892982/SRR892982_2.fastq.gz
   
Unzip these files.

.. code-block:: none

   gunzip SRR892982_1.fastq.gz
   gunzip SRR892982_2.fastq.gz
   

Index the fasta file using the Bs seeker indexer. You can either do this by running a command for Bs Seeker, or simply change the path and file name in test_bs_seeker_indexer.py in the following lines 

.. code-block:: none

   resource_path = os.path.join(os.path.dirname(__file__), "data/") # change to your file path 
   genomefa_file = resource_path + "bsSeeker.Mouse.GRCm38.fasta"" # change to name of the fasta file 
   
   
Then run 

.. code-block:: none
   pytest test_bs_seeker_indexer.py
   

Filter the fastq files. You can again do this by changing the following lines within test_bs_seeker_filter.py for fastq 1 in test_bs_seeker_filter_01():

.. code-block:: none

   resource_path = os.path.join(os.path.dirname(__file__), "data/") # change to your file path 
   genomefa_file = resource_path + "bsSeeker.Mouse.GRCm38_1.fastq" # change to name of the fasta file 
   
   
You would need to repeat this for the 2nd fastq file within test_bs_seeker_filter_02():

Then run 

.. code-block:: none
   pytest test_bs_seeker_filter.py
   

Split the filtered fastq.s. Alter the following lines in test_fastq_splitter.py : 

.. code-block:: none

    resource_path = os.path.join(os.path.dirname(__file__), "data/") # change to your file path
    fastq_1file = resource_path + "SRR892982_1.filtered.fastq" # change to your file path for the filtered file
    fastq_2file = resource_path + "SRR892982_2.filtered.fastq" # change to your file path for the filtered file

Then run 

.. code-block:: none
   pytest test_fastq_splitter.py
       
    
Align the fastq files to the indexed fasta. Change the following lines in test_bs_seeker_aligner.py

.. code-block:: none

    resource_path = os.path.join(os.path.dirname(__file__), "data/") # change to your file path 
    genomefa_file = resource_path + "bsSeeker.Mouse.GRCm38.fasta" # change to name of the fasta file 
    genomeidx_file = resource_path + "bsSeeker.Mouse.GRCm38.fasta_bowtie2.tar.gz" # change to your file path for the indexed files
    fastq_gz = resource_path + "bsSeeker.Mouse.GRCm38_1.fastq.filtered.fastq.tar.gz" # change to your file path for the filtered/splitted files
    fastq1_file = "bsSeeker.Mouse.GRCm38_1.fastq.filtered.fastq" # change to your file path for the filtered file
    fastq2_file = "bsSeeker.Mouse.GRCm38_2.fastq.filtered.fastq" # change to your file path for the filtered file
    
        
Then run 

.. code-block:: none
   pytest test_bs_seeker_aligner.py
   

Aligning would give out the bam and bai files. 

Run the methylation caller. Alter the following lines within test_bs_seeker_methylation_caller.py

.. code-block:: none

    resource_path = "/nfs/nobackup/ensembl/reham_ens/BS_seeker_tests/" #os.path.join(os.path.dirname(__file__), "data/")
    genome_fa_file = resource_path + "Mouse.selected_region.fasta_bowtie2.tar.gz"#"bsSeeker.Mouse.GRCm38.fasta_bowtie2.tar.gz"
    bam_file = resource_path + "bsSeeker.Mouse.GRCm38.bam"
    

This would give the wig file. 

Traverse wig file for a suitable region 

Select chromosomal region corresponding to the above

Re run pipeline till aligner.     
    

And make the sam file

.. code-block:: none

   bwa samse GCA_000001405.22.chr22.fa.fasta GCA_000001405.22.chr22.sai DRR000150.chr22.fastq >GCA_000001405.22.chr22.sam


===========



.. code-block:: none

   python traverseForCoverageRegion_ChIP.py path/to/GCA_000001405.22.chr22.dp

Running this script would print the spanning regions. If it is a continuous region, you may only take the first starting base pair and the last ending base pair, as inputs for the next step. (Take out 1000 and add in 1000 to these respectively to get upstream and downstream spanning bases)

Extract the corresponding fasta sequence from the chromosome file for the positions retrieved from the above step. Checkout file from https://github.com/Multiscale-Genomics/mg-misc-scripts/blob/master/ChIPSeq_Scripts/extractChromosomalRegion.py and run using command:

.. code-block:: none

   python extractChromosomalRegion.py path/to/original/fasta/file path/to/output/file/for/region/macs2.Human.GCA_000001405.22.fasta starting_base_position (39112298) ending_base_position (39112402)


Making the Fastq file
^^^^^^^^^^^^^^^^^^^^^^

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
