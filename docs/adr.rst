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


2018-05-11 - Sleuth gene differential analysis pipeline
-------------------------------------------------------

This allows for the comparison of multiple RNA-seq experiments to determine if there are any genes that are differentially expressed. This has required changes to the output of the kallisto_quant tool so that it generates only a single tar file containing the abundance and run_info files. There is also the introduction of the bootstrap-sample parameter as part of the quantification to determine the accuracy of the counts.

The first tool uses Sleuth to generate an R object of all the processed tracks. Separate tools are written for each visualisation to allow for a certain amount of parallelisation with the results being saved to an archive file.


2018-05-22 - GEM Naming
-----------------------

Update so that the gem files are name <genome-file>.gem.gz inline with requests from WP7 partners so that the name of the index is picked up correctly


2018-05-22 - TrimGalore
-----------------------

To try and improve the quality of the reads that are used for numerous pipelines, TrimGalore has been included as a pipeline to aid in the clipping and removal of low quality regions of reads. The pipeline can be run on single or paired end FASTQ files. A report of the trimmed data is also returned for the user to identify what changes were made.
