.. Copyright 2017 EMBL-European Bioinformatics Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

FASTQ Functions
===============

The following functions are ones that are used for the manipulation of FASTQ
files.

.. automodule:: tool

   Reading
   -------

   The following functions are to provide easy access for iterating through entries
   within a FASTQ file(s) both single and paired.

   .. autoclass:: fastqreader.fastqreader
     :members:

   Splitting
   ---------

   This tool has been created to aid in splitting FASTQ files into manageable
   chunks for parallel processing. It is able to work on single and paired end
   files.

   .. autoclass:: tool.fastq_splitter
     :members:
