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

    resource_path = os.path.join(os.path.dirname(__file__), "data/") # change to your file path 
    genome_fa_file = resource_path + "bsSeeker.Mouse.GRCm38.fasta_bowtie2.tar.gz" # change to your file path for the indexed files
    bam_file = resource_path + "bsSeeker.Mouse.GRCm38.bam" # change to your file for an intermediate bam file
    

This would give the wig file. Traverse this wig file for a suitable region. Download the script : 

.. code-block:: none

   https://github.com/Multiscale-Genomics/mg-misc-scripts/blob/master/WGBS_Scripts/regionsFromWig.py
   
Run it using : 

.. code-block:: none

   python regionsFromWig.py /path/to/bsSeeker.Mouse.GRCm38.wig
   
Select the following coordinates from the output 55844491 55847491

Select chromosomal region corresponding to the above by getting the script at : 

.. code-block:: none

   https://github.com/Multiscale-Genomics/mg-misc-scripts/blob/master/ChIPSeq_Scripts/extractChromosomalRegion.py

And run it using 

.. code-block:: none

   python extractChromosomalRegion.py /path/to/Mouse.CM001012.2.fasta /path/to/output/bsSeeker.Mouse.GRCm38.fasta 55844491 55847491


Re run the pipeline with this fasta file and original fastq files till the alignment step. Take the .bam file and convert it to sam using : 

.. code-block:: none

   samtools view -h -o /path/to/BS_seeker_tests/bsSeeker.Mouse.GRCm38.sam  /path/to/bsSeeker.Mouse.GRCm38.bam
   
Use this sam file to extract the fastq entries from the larger fastq files. 

.. code-block:: none

   python makeFastQFiles.py --samfile /path/to/bsSeeker.Mouse.GRCm38.sam --fastQfile /path/to/SRR892982_1.fastq --pathToOutput /path/to/output/ --fastqOut bsSeeker.Mouse.GRCm38_1.fastq
   
   python makeFastQFiles.py --samfile /path/to/bsSeeker.Mouse.GRCm38.sam --fastQfile /path/to/SRR892982_2.fastq --pathToOutput /path/to/output/ --fastqOut bsSeeker.Mouse.GRCm38_2.fastq

   
The fastq files in the above steps along with the bsSeeker.Mouse.GRCm38.fasta make up the test data for the WGBS pipeline.