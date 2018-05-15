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

Sleuth Test Data
================

Test Data
---------

Dataset
^^^^^^^

+------------+--------------------------------------------------------------+
| Stable IDs | ERR030856, ERR030857, ERR030858, ERR030872, ERR030903        |
+------------+--------------------------------------------------------------+
| Project    | `PRJEB2445 <https://www.ebi.ac.uk/ena/data/view/PRJEB2445>`_ |
+------------+--------------------------------------------------------------+

Genome
^^^^^^

CDNA was downloaded from `ensembl 92 <http://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz>`_

+-------------+--------+
| Assembly    | GRCh38 |
+-------------+--------+
| Transcripts | 1000   |
+-------------+--------+

Method
------
The full dataset was downloaded from ENA aligned to the cDNA using kallisto producing a pseudo alignment bam file. Sleuth was used to calculate the most significant hits of which the top 1000 were picked. These were used to select the matching FASTQ reads from the pseudo alignment files


hiseq_info.txt file:

.. code-block:: none
   :linenos:

   ERR030856\tcontrol
   ERR030857\tcontrol
   ERR030858\tcontrol
   ERR030872\tthyroid
   ERR030903\tthyroid

.. code-block:: R
   :linenos:

   library("sleuth")
   sample_id <- dir(file.path("data", "results"))
   kal_dirs <- file.path("data", "results", sample_id, "kallisto")

   s2c <- read.table(file.path("data", "hiseq_info.txt"), header = TRUE, stringsAsFactors=FALSE)
   s2c <- dplyr::select(s2c, sample, condition)
   s2c <- dplyr::mutate(s2c, path = kal_dirs)

   so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, num_cores = 1)
   so <- sleuth_fit(so, ~condition, 'full')
   so <- sleuth_fit(so, ~1, 'reduced')
   so <- sleuth_lrt(so, 'reduced', 'full')

   sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
   sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

   # Generate a set of transcripts to use for code testing
   sample(sleuth_significant$target_id, 1000)


.. code-block:: none
   :linenos:

   # Kallisto Quantification
   kallisto quant -i GRCh38.cdna.fasta.idx -o ERR030856 --pseudobam --single -l 100 -s 0.01 ERR030856/ERR030856.fastq > ERR030856/ERR030856.sam
   kallisto quant -i GRCh38.cdna.fasta.idx -o ERR030857 --pseudobam --single -l 100 -s 0.01 ERR030857/ERR030857.fastq > ERR030857/ERR030857.sam
   kallisto quant -i GRCh38.cdna.fasta.idx -o ERR030858 --pseudobam --single -l 100 -s 0.01 ERR030858/ERR030858.fastq > ERR030858/ERR030858.sam
   kallisto quant -i GRCh38.cdna.fasta.idx -o ERR030903 --pseudobam --single -l 75 -s 0.0133333 ERR030903/ERR030903.fastq > ERR030903/ERR030903.sam


.. code-block:: none
   :linenos:

   # Extract the FASTQ read IDs for the selected transcripts
   grep -f sleuth_sample_transcripts.txt ERR030872/ERR030872.sam | tr "\t" "~" | cut -d"~" -f1 | grep -v @ > ERR030872/ERR030872.reads
   grep -f sleuth_sample_transcripts.txt ERR030903/ERR030903.sam | tr "\t" "~" | cut -d"~" -f1 | grep -v @ > ERR030903/ERR030903.reads
   grep -f sleuth_sample_transcripts.txt ERR030856/ERR030856.sam | tr "\t" "~" | cut -d"~" -f1 | grep -v @ > ERR030856/ERR030856.reads
   grep -f sleuth_sample_transcripts.txt ERR030857/ERR030857.sam | tr "\t" "~" | cut -d"~" -f1 | grep -v @ > ERR030857/ERR030857.reads
   grep -f sleuth_sample_transcripts.txt ERR030858/ERR030858.sam | tr "\t" "~" | cut -d"~" -f1 | grep -v @ > ERR030858/ERR030858.reads


.. code-block:: none
   :linenos:

   # Extract the original reads from teh FASTQ files
   python scripts/ExtractRowsFromFASTQs.py --input_1 ERR030856/ERR030856.fastq --rows ERR030856/ERR030856.reads --prop 0.1 --output_tag subset
   python scripts/ExtractRowsFromFASTQs.py --input_1 ERR030857/ERR030857.fastq --rows ERR030857/ERR030857.reads --prop 0.1 --output_tag subset
   python scripts/ExtractRowsFromFASTQs.py --input_1 ERR030858/ERR030858.fastq --rows ERR030858/ERR030858.reads --prop 0.1 --output_tag subset
   python scripts/ExtractRowsFromFASTQs.py --input_1 ERR030872/ERR030872_1.fastq --input_2 ERR030872/ERR030872_2.fastq --rows ERR030872/ERR030872.reads --prop 0.1 --output_tag subset
   python scripts/ExtractRowsFromFASTQs.py --input_1 ERR030903/ERR030903.fastq --rows ERR030903/ERR030903.reads --prop 0.1 --output_tag subset

Due to the number of reads that match to the transcripts, only 1% have been kept for code testing
