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


2018-05-11 - Sleuth gene differential analysis pipeline
-------------------------------------------------------

This allows for the comparison of multiple RNA-seq experiments to determine if there are any genes that are differentially expressed. This has required changes to the output of the kallisto_quant tool so that it generates only a single tar file containing the abundance and run_info files. There is also the introduction of the bootstrap-sample parameter as part of the quantification to determine the accuracy of the counts.

The first tool uses Sleuth to generate an R object of all the processed tracks. Separate tools are written for each visualisation to allow for a certain amount of parallelisation with the results being saved to an archive file.