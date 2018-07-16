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

WGBS Test Data
==============

Test Data
---------

Dataset
^^^^^^^

+-----------+----------------------------------------------------------------+
| Stable ID | SRR892982                                                      |
+-----------+----------------------------------------------------------------+
| Citation  | `Sun et al 2014 <http://europepmc.org/abstract/MED/24792119>`_ |
+-----------+----------------------------------------------------------------+

Genome
^^^^^^

+------------+---------------+
| Assembly   | GRChm8, chr19 |
+------------+---------------+
| Chromosome | 19            |
+------------+---------------+
| Start      | 3000000       |
+------------+---------------+
| End        | 4020000       |
+------------+---------------+

Method
------
The full dataset was downloaded from ENA aligned to the genome using GEM.

.. code-block:: none
   :linenos:

   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR892/SRR892982/SRR892982_1.fastq.gz
   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR892/SRR892982/SRR892982_2.fastq.gz
   gunzip SRR892982_1.fastq.gz
   gunzip SRR892982_2.fastq.gz

   wget "http://www.ebi.ac.uk/ena/data/view/CM001012&display=fasta&download=fasta&filename=entry.fasta" -O mouse.GRCm38.19.fasta

   sed -n 40022200,41046060p GCA_000001635.6.fasta > mouse.GRCm38.19.fasta

   bowtie2-build mouse.GRCm38.19.fasta mouse.GRCm38.19.idx
   bowtie2 -x genome/mouse.GRCm38.19.idx -U data/SRR892982_1.fastq > SRR892982_1.sam
   samtools view -h -o SRR892982_1.sam SRR892982_1.bam

   awk 'BEGIN { FS = "[[:space:]]+" } ; $3 ~ /CM001012/ {print $1}' SRR892982_1.sam > SRR892982_1.chr19.rows

   python scripts/ExtractRowsFromFASTQs.py --input_1 tests/data/SRR892982_1.fastq --input_2 tests/data/SRR892982_2.fastq --rows SRR892982_1.chr19.rows --output_tag profile

   mv data/tmp/SRR892982_1.profile_0.fastq data/SRR892982_1.chr19.fastq
   mv data/tmp/SRR892982_2.profile_0.fastq data/SRR892982_2.chr19.fastq

   bowtie2 -x genome/mouse.GRCm38.19.idx -1 data/SRR892982_1.chr19.fastq -2 data/SRR892982_2.chr19.fastq > SRR892982_1.chr19.sam

There are a range of positive and negative peaks between 3000000 and 4020000. Subselect the genome and matching FASTQ reads for the required subsection:

.. code-block:: none
   :linenos:

   head -n 1 mouse.GRCm38.19.fasta > 3M/mouse.GRCm38.19.3M.fasta
   sed -n 50001,67001p mouse.GRCm38.19.fasta >> 3M/mouse.GRCm38.19.3M.fasta

   bowtie2-build 3M/mouse.GRCm38.19.3M.fasta 3M/mouse.GRCm38.19.3M.idx
   bowtie2 -p 4 -x 3M/mouse.GRCm38.19.3M.idx -1 SRR892982_1.chr19.fastq -2 SRR892982_2.chr19.fastq > 3M/SRR892982.chr19.3M.sam
   samtools view -h -o 3M/SRR892982.chr19.3M.sam 3M/SRR892982.chr19.3M.bam

   python scripts/ExtractRowsFromFASTQs.py --input_1 SRR892982_1.chr19.fastq --input_2 SRR892982_2.chr19.fastq --rows 3M/SRR892982.chr19.3M.rows --output_tag subset

   mv tmp/SRR892982_1.chr19.subset_0.fastq 3M/SRR892982_1.chr19.3M.fastq
   mv tmp/SRR892982_2.chr19.subset_0.fastq 3M/SRR892982_2.chr19.3M.fastq

This has enough load to test the system, while also generating results in the output files.

These files are then saved in the tests/data directory as:

.. code-block:: none
   :linenos:

   bsSeeker.Mouse.GRCm38.fasta
   bsSeeker.Mouse.SRR892982_1.fastq.gz
   bsSeeker.Mouse.SRR892982_2.fastq.gz

These files were too large for use within the Travis tst environment, so the number of entries was reduced by taking every other read:

.. code-block:: none
   :linenos:

   sed -n -e '0~9{N;N;N;N;p}' tests/data/bsSeeker.Mouse.SRR892982_1.fastq > tests/data/test_bsSeeker.Mouse.SRR892982_1.fastq
   sed -n -e '0~9{N;N;N;N;p}' tests/data/bsSeeker.Mouse.SRR892982_2.fastq > tests/data/test_bsSeeker.Mouse.SRR892982_2.fastq

   mv tests/data/test_bsSeeker.Mouse.SRR892982_1.fastq tests/data/bsSeeker.Mouse.SRR892982_1.fastq
   mv tests/data/test_bsSeeker.Mouse.SRR892982_2.fastq tests/data/bsSeeker.Mouse.SRR892982_2.fastq

   gzip tests/data/bsSeeker.Mouse.SRR892982_1.fastq
   gzip tests/data/bsSeeker.Mouse.SRR892982_2.fastq

Test Scripts
------------

The following are the tests for checking that the tools in the WGBS pipeline are
functioning correctly.

The tests should be run in this order so that the required input files are
generated at the correct stage.

.. code-block:: none
   :linenos:

   pytest -m wgbs tests/test_fastqc_validation.py
   pytest -m wgbs tests/test_bs_seeker_filter.py
   pytest -m wgbs tests/test_bs_seeker_indexer.py
   pytest -m wgbs tests/test_bs_seeker_aligner.py
   pytest -m wgbs tests/test_bs_seeker_methylation_caller.py

These can be called as part of a single tool chain with:

.. code-block:: none
   :linenos:

   python tests/test_toolchains.py --pipeline wgbs