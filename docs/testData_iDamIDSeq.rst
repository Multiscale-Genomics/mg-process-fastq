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

Hi-C Test Data
==============

Test Data
---------

Dataset
^^^^^^^

+-----------+----------------+
| Stable ID | SRR1658573     |
+-----------+----------------+
| Citation  | Rao et al 2014 |
+-----------+----------------+

Genome
^^^^^^

+------------+----------+
| Assembly   | GRCh38   |
+------------+----------+
| Chromosome | 22       |
+------------+----------+
| Start      | 15000000 |
+------------+----------+
| End        | 19999999 |
+------------+----------+

Method
------
The full dataset was downloaded from ENA aligned to the genome using GEM.

.. code-block:: none
   :linenos:

   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR371/005/SRR3714775/SRR3714775.fastq.gz
   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR371/005/SRR3714775/SRR3714776.fastq.gz
   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR371/005/SRR3714775/SRR3714777.fastq.gz
   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR371/005/SRR3714775/SRR3714778.fastq.gz

   wget "http://www.ebi.ac.uk/ena/data/view/CM000684.2&display=fasta&download=fasta&filename=entry.fasta" -O GCA_000001405.22.chr22.fasta

   bwa index GCA_000001405.22.chr22.fasta

   bwa aln -q 5 -f SRR3714775.fastq.sai
   bwa aln -q 5 -f SRR3714776.fastq.sai
   bwa aln -q 5 -f SRR3714777.fastq.sai
   bwa aln -q 5 -f SRR3714778.fastq.sai

   bwa samse -f SRR3714775.fastq.sam GCA_000001405.22.chr22.fasta SRR3714775.fastq.sai SRR3714775.fastq
   bwa samse -f SRR3714776.fastq.sam GCA_000001405.22.chr22.fasta SRR3714776.fastq.sai SRR3714776.fastq
   bwa samse -f SRR3714777.fastq.sam GCA_000001405.22.chr22.fasta SRR3714777.fastq.sai SRR3714777.fastq
   bwa samse -f SRR3714778.fastq.sam GCA_000001405.22.chr22.fasta SRR3714778.fastq.sai SRR3714778.fastq

   samtools view -h -o SRR3714775.sam SRR3714775.bam
   samtools view -h -o SRR3714776.sam SRR3714776.bam
   samtools view -h -o SRR3714777.sam SRR3714777.bam
   samtools view -h -o SRR3714778.sam SRR3714778.bam

   awk 'BEGIN { FS = "[[:space:]]+" } ; $3 ~ /CM000684/ {print $1}' SRR3714775.sam > SRR3714775.chr22.sam
   awk 'BEGIN { FS = "[[:space:]]+" } ; $3 ~ /CM000684/ {print $1}' SRR3714776.sam > SRR3714776.chr22.sam
   awk 'BEGIN { FS = "[[:space:]]+" } ; $3 ~ /CM000684/ {print $1}' SRR3714777.sam > SRR3714777.chr22.sam
   awk 'BEGIN { FS = "[[:space:]]+" } ; $3 ~ /CM000684/ {print $1}' SRR3714778.sam > SRR3714778.chr22.sam

   python scripts/ExtractRowsFromFASTQs.py --input_1 tests/data/SRR3714775.fastq --rows tests/data/SRR3714775.chr22.sam --output_tag profile
   python scripts/ExtractRowsFromFASTQs.py --input_1 tests/data/SRR3714776.fastq --rows tests/data/SRR3714776.chr22.sam --output_tag profile
   python scripts/ExtractRowsFromFASTQs.py --input_1 tests/data/SRR3714777.fastq --rows tests/data/SRR3714777.chr22.sam --output_tag profile
   python scripts/ExtractRowsFromFASTQs.py --input_1 tests/data/SRR3714778.fastq --rows tests/data/SRR3714778.chr22.sam --output_tag profile

The alignments were then filtered with BioBamBam and peak calling was performed with iDEAR and a suitable region with a number of peaks was identified. The chromosomal region was extract and the matching reads to this region were identified. To reduce the number of reads that matched so that it could be used as a test set for the code base every other read was selected so that a reasonable number of reads we present. This mean that there are results when running the test, but generating the alignments does not take too long to compute.

.. code-block:: none
   :linenos:
   sed -n 539916,556583p GCA_000001405.chr22.fasta > GCA_000001405.chr22ss.fasta
   sed -n 5002,11669p idear.Human.GCA_000001405.22.fasta >> idear.Human.GCA_000001405.22.subset.fasta


Test Scripts
------------

The following are the tests for checking that the tools in the Hi-C pipeline are
functioning correctly.

The tests should be run in this order so that the required input files are
generated at the correct stage.

.. code-block:: none
   :linenos:

   pytest -m idamidseq tests/test_bwa_indexer.py
   pytest -m idamidseq tests/test_bwa_aligner.py
   pytest -m idamidseq tests/test_biobambam.py
   pytest -m idamidseq tests/test_bsgenome.py
   pytest -m idamidseq tests/test_idear.py

These can be called as part of a single tool chain with:

.. code-block:: none
   :linenos:

   python tests/test_toolchains.py --pipeline idamidseq