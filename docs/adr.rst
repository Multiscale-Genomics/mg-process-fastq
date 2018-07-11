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

BioBamBam only marks reads as duplicate, but does not remove the after. The Tool has been updated to remove the flagged duplicates using samtools with the parameter `-F 1024`. This matches the pipeline used within the `Blueprints project <http://dcc.blueprint-epigenome.eu/#/md/chip_seq_grch37>`_.

Also performed some tidying of the code to annotate issues that had been highlighted by pylint.


2018-07-11 - Changes FASTQ splitter file management
---------------------------------------------------

The previous splitter would split the FASTQ files into separate changes, then create the tar file and then the Gzip file. This results in a large amount of wasted tmp space, which is a limited resource. The changes implemented incrementally add the sub-FASTQ files to the archive file, deleting them once they have been added. The whole archive file is then compressed. This has a large advantage when handling larger human datasets.
