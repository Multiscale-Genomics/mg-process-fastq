.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Architectural Design Record
===========================

2017-08-15 - Implementation of pigz
-----------------------------------

It was highlighted that the archival step using the python module tarfile was taking a long time when archiving and compressing large volumes of data. This is an issue in the WGBS pipeline where the FASTQ files are broken down into small sections and are then aligned individually.

The issue was in relation to the compression step. tar, and therefore tarfile, runs as a single process. In linear time to compress a 14GB file takes 18 minutes, which is a long time for the user that have to wait. The archival step is less than 1 min. As a result the compression has moved to using `pigz <https://zlib.net/pigz/>`_ after the archival step to perform the compression in a parallel fashion across all available cores on a machine. On a 4 core machine this allowed the compression of the same file to take only 7 min.

The resultant compressed file remains accessible via gzip. After decompression files have the same structure, so this should be an invisible change.


2018-01-26 - Disable no-self-use for @tasks
-------------------------------------------

The disabling of pylint test for determining if a function should be a static method was made as this affected the behaviour of pycompss preventing it from functioning correctly. As a result all @task functions do not require the no-self-use test to be run.


2018-02-28 - BAM Merge Strategy
-------------------------------

Based on benchmarking of the pipeline the procedure for merging of bam files has been modified to get an optimal balance between running as much as possible in parallel vs the cost of spinning up new jobs to perform the merging. It was found that running each job merging 10 files provided the break even point between the cost of creating a new job and getting the maximum throughput through the system. It also reduces the number of iterative merge procedures which is beneficial when there are alignments that are difficult to merge.


2018-04-26 - BAM Splitting
--------------------------

Added the required functions for tools to be able to split a bam by the number of chromosomes so that the analysis can be done in parallel. As an initial run this has been implemented in the MACS2 tool, where the indexing of the bam file is performed by the tool. The bam and bai files are passed to the jobs so that they can then extract the chromosome that they are required to analyse.

In the future the creation of the bai file could be done at alignment time, but in the case of MACS2 there is a filtering step on the aligned bam file, so a new index would be required.


2018-05-01 - Compression of FASTQ
---------------------------------

Added compression of the split FASTQ files to reduce the amount of space required when processing the data. There is also the removal of the tmp directory after the completion of the splitter to save on space.


2018-05-09 - Handling aligner index decompression
-------------------------------------------------

The code has been modified so that there is a single decompression of the BWA and Bowtie2 common indexes. The index files are then explicitly handed to the alignment task rather than handing over the compressed index. The decompression is performed as a @task so that the index files are already in the COMPSs system. This means that handing the index files to the alignment tasks creates a single symlink in the sandbox temporary file directory rather than duplicating the whole of the index structure for each job.


2018-05-22 - GEM Naming
-----------------------

Update so that the gem files are name <genome-file>.gem.gz inline with requests from WP7 partners so that the name of the index is picked up correctly


2018-05-22 - TrimGalore
-----------------------

To try and improve the quality of the reads that are used for numerous pipelines, TrimGalore has been included as a pipeline to aid in the clipping and removal of low quality regions of reads. The pipeline can be run on single or paired end FASTQ files. A report of the trimmed data is also returned for the user to identify what changes were made.


2018-05-31 - Public genomes and indexes
---------------------------------------

The VRE is making it possible to use centrally stored genomic sequences and the pre-built indexes. This relies on a second key-value pair for the genome and indexes. Pipelines that interact with the public files need to be able to test for the present of "genome_public" and "index_public", if they are not present then they should be able to safely fall back on "genome" and "index" keys-value pairs. With the separation of the tools and the pipelines this means that the changes should only need to happen in the pipelines as this can perform the translation between the "genome_public" to "genome" ready for the tool.


2018-06-01 - Separated WGBS Code Testing
----------------------------------------

To bring down the run time for the TravisCI, the WGBS has been moved to a separate track. This has the benefit of getting the testing started earlier and allowing the other tests to finish sooner.


2018-06-01 - Travis Caching
---------------------------

Travis CI is now able to cache the lib directory so that the build phase is reduced to improve test speeds. If there are changes to the lib then these need to be refreshed in the TravisCI settings to ensure that the new libraries are included, or flushed if there are changes to the versions of packages in the lib.

There is also caching of the pip directory to reduce the load time.


2018-06-04 - Split the WGBS test scripts
----------------------------------------

Split the testing of the WGBS pipeline and tool chains so that they 2 sets can run in parallel. Both take too long when run in series.


2018-06-05 - Use of the logger PROGRESS
---------------------------------------

Added in the use of the logger.progress to indicate the progression of a process.


2018-06-14 - Paired end alignment
---------------------------------

The aligner pipelines has been modified the pass through all the input and metadata to the aligner tools, this simplifies the the passing of a second fastq file and also make using these pipelines for alignment of paired end data possible.


2018-06-18 - Branch tidying during alignment
--------------------------------------------

Modified the way that the alignment pipelines manage the temporary files. These are now deleted once the pipeline has finished using them. The purpose of this is to save space on the file system and prevent large jobs taking up too much space.

There have also been changes to the handling of paired end files for the alignment pipelines improving the clarity of what is happening and simplifying the passing of parameters. There are also changes to the tests to allow for the removal of temporary files and there are tests to make sure that the output bam files are single or paired end.

Other changes include:
- Simplification of the untarring functions
- Modifications to the Bowtie2 index file for consistency with the BWA index file
- Refactored the BWA ALN sai file generation to reduce redundancy to allow for multi-processing when there is paired-end data
- Improved the handling of the suffixes for FASTQ and FASTA files so that it can handle variants


2018-06-27 - Remove reads marked as duplicate by BioBamBam
----------------------------------------------------------

BioBamBam only marks reads as duplicate, but does not remove them after. The Tool has been updated to remove the flagged duplicates using samtools with the parameter `-F 1024`. This matches the pipeline used within the `Blueprints project <http://dcc.blueprint-epigenome.eu/#/md/chip_seq_grch37>`_.

Also performed some tidying of the code to annotate issues that had been highlighted by pylint.


2018-07-05 - Refactoring of repo to avoid naming collisions
-----------------------------------------------------------

The original repo had the tools and tests in a directory that was in the root of the repo. This is problematic when there is sharing of the code as there are collisions in the name space. This has meant that all of the code needs to be moved into a new subdirectory (mg_process_fastq).

There is also the movement of the workflow classes out of the pipeline scripts and into their own module in mg_process_fastq. This means that tworkflows that have been written by others can be more easily imported.

This change will probably necessitate incrementing the major release number to 1.0.0 as this is not a backwards compatible change and will require changes in other repos that rely on mg-process-fastq. As a result this will need to be part of a planned release with other developers.


2018-07-11 - Changes FASTQ splitter file management
---------------------------------------------------

The previous splitter would split the FASTQ files into separate changes, then create the tar file and then the Gzip file. This results in a large amount of wasted tmp space, which is a limited resource. The changes implemented incrementally add the sub-FASTQ files to the archive file, deleting them once they have been added. The whole archive file is then compressed. This has a large advantage when handling larger human datasets.

There has also been some refactoring of the handling of the archiving and compression steps to reduce the duplication of code within the repository.


2018-07-16 - Modified handling of file locations
------------------------------------------------

Updated the handling of file locations to use os.path.join and os.path.split to allow for compatibility between different operating systems for the pipelines and tools.


2018-08-02 - Added in Paired End BAM file handling for MACS2
------------------------------------------------------------

MACS2 is able to automatically handle the files that are handed to it except for paired-end BAM and BED files (BAMPE and BEDPE respectively). The MACS2 tool only accepts BAM files so a check was implemented to determine if the BAM file contained paired-end reads.

There has also been a major rewrite of the MACS2 tool to remove code duplication.


2018-07-16 - Modified handling of file locations
------------------------------------------------

Updated the handling of file locations to use os.path.join and os.path.split to allow for compatibility between different operating systems for the pipelines and tools.


2018-08-07 - Storing tool parameters as part of the metadata
------------------------------------------------------------

To improve the amount of information that is stored about the run of a tool, the parameters that were used are now being included as part of the metadata.


2018-08-07 - Extra output files from MACS2
------------------------------------------

MACS2 is able to generate a plot of the results as well as a bedGraph. These have now been integrated as part of the output files from teh tool.


2018-08-13 - Normalised the use of OSError
------------------------------------------

IOError was depricated in favour of OSError when moving to py3, but to maintian backwards compatibility IOError also needs to be supported. There were places in the code where this was not true and other places that relied on just OSError. Instances of just IOError have been converted to testing for both IOError and OSError and visa versa.


2018-08-15 - Use the config.json execution path
-----------------------------------------------

Using the directory of the input file for building the location of the working directory with outside of a task is not a viable option as it can write data to the wrong directory. The execution path provided in the config.json arguments is the right place. This location is also the location for output files. This issue occurred as the FASTQ splitter was generating a tar file that the aligners were downloading to the wrong location. Even though this was tidied up this was still not the right place to put this file.


2018-08-16 - Prevent further duplicate filtering by MACS2
---------------------------------------------------------

In the process_chipseq.py pipeline the duplicates have already been filtered by BioBamBam2 and samtools so there is no need for further filtering to be done by MACS2.


2018-09-04 - Adding functionality to bam_utils and MACS2
--------------------------------------------------------

MACS2 was previously set to work with the BAMPE option for the -f/--format parameter. Additional functionality has been added to bam_utils and macs2 mode to incorporate the BEDPE option. This has been done for the Atac Seq pipeline to incorporate the processing of bed file rather than bam files if the user would need changes to the result files generated.


2018-09-17 - Updates to tool and pipline run()
----------------------------------------------

Changes to the pipelines so that the run() function matches the definitions within the Tool API. There have also been a number of changes so that the pipeline and tool code is python 3 compatible


2018-08-22 - Improvement of tadbit tools wrappers
-------------------------------------------------

A json with the matrix was included in the outputs of tadbit bin
New normalization method OneD in tadbit normalize
Code update to use last features of the development branch of tadbit tools api
The wrapper of tadbit model was rebuilt to allow the modelling of full genomes, mainly for yeast
General reshape of all the code according to pylint
Inclusion of tests for the wrappers and tools of the tadbit pipelines


2018-09-25 - Converting the Kallisto TSV file to BED
----------------------------------------------------

To display the scores on the genome browser the abundance tsv is used to generate a bed file where the score matches the transcripts per million column from the abundance.tsv output from Kallisto. This module requires the presence of the ensembl gff3 file for the matching assembly. This should be passed by the VRE when passing the FASTA file for the transcripts.


2018-10-18 - Multi File handling for the DamID-seq Pipeline
-----------------------------------------------------------

Ability to handle multiple input single/paired end data and background data files and process them in an orderly fashion. Ability to handle the resultant multiples of generated bam files in the idear tool and its matching individual pipeline.


2018-10-25 - WGBS Pipeline Create BigWig files as standard
----------------------------------------------------------

The output wig files from BS Seeker2 are now converted to BigWig files by default rather than returning wig files. This is so that they are easier to visualise on the JBrowse interface.


2018-10-24 - ChIP-seq Pipeline to use BWA MEM
---------------------------------------------

Changed the default aligner for the ChIP-seq pipeline from BWA ALN to BWA MEM for speed improvements and better handling of short read data.
