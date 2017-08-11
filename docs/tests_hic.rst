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
| Chromosome | 21       |
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

   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR165/003/SRR1658573/SRR1658573_1.fastq.gz
   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR165/003/SRR1658573/SRR1658573_2.fastq.gz

   wget "http://www.ebi.ac.uk/ena/data/view/CM000663.2,CM000664.2,CM000665.2,CM000666.2,CM000667.2,CM000668.2,CM000669.2,CM000670.2,CM000671.2,CM000672.2,CM000673.2,CM000674.2,CM000675.2,CM000676.2,CM000677.2,CM000678.2,CM000679.2,CM000680.2,CM000681.2,CM000682.2,CM000683.2,CM000684.2,CM000685.2,CM000686.2,J01415.2&display=fasta&download=fasta&filename=entry.fasta" -O GCA_000001405.22.fasta

.. code-block:: python
   :linenos:

   from tool.common import common

   genome_file = "GCA_000001405.22.fasta"
   new_genome_file = "GCA_000001405.22_gem.fasta"

   common_handle = common()
   common_handle.replaceENAHeader(genome_file, new_genome_file)

   idx_loc = common_handle.gem_index_genome(new_genome_file, new_genome_file)

   from pytadbit.mapping.mapper import full_mapping

   fastq_file_1 = "SRR1658573_1.fastq"
   fastq_file_2 = "SRR1658573_2.fastq"
   output_dir = "mapped_reads"
   windows = None
   enzyme_name = "MboI"

   map_files_1 = full_mapping(
       idx_loc, fastq_file, output_dir,
       r_enz=enzyme_name, windows=windows, frag_map=True, nthreads=32,
       clean=True, temp_dir='/tmp/'
   )
   map_files_2 = full_mapping(
       idx_loc, fastq_file, output_dir,
       r_enz=enzyme_name, windows=windows, frag_map=True, nthreads=32,
       clean=True, temp_dir='/tmp/'
   )


A list of FASTQ ids that had aligned to the genome were then extracted:

.. code-block:: none
   :linenos:

   grep chr21 mapped_reads/SRR1658573_1_full_1-100.map | tr "\t" "~" | cut -d"~" -f1 -f5 | tr ":" "\t" | awk '(NR==1) || (($4>15000000) && ($4<20000000))' | tr "\t" "~" | cut -d "~" -f1 > SRR1658573_1_chr21_1-100.row
   grep chr21 mapped_reads/SRR1658573_2_full_1-100.map | tr "\t" "~" | cut -d"~" -f1 -f5 | tr ":" "\t" | awk '(NR==1) || (($4>15000000) && ($4<20000000))' | tr "\t" "~" | cut -d "~" -f1 > SRR1658573_2_chr21_1-100.row

The IDs were filtered to ensure matching pairs in both files:

.. code-block:: none
   :linenos:

   grep -Fx -f SRR1658573_1_chr21_1-100.row SRR1658573_2_chr21_1-100.row > SRR1658573_chr21.row

The `split_paired_fastq.py` was used to divide the original FASTQ files into
chunks of 1000000 reads. The `ExtractRowsFromFASTQ.py` script was then used to
extract the matching FASTQ pairs from each of the FASTQ files in parallel. All
of the individual FASTQ files were then concatenated together to form the final
2 FASTQ test files. The commands for this were:

.. code-block:: none
   :linenos:

   python split_paired_fastq.py --input_1 SRR1658573_1.fastq --input_2 SRR1658573_1.fastq


Test Scripts
------------

The following are the tests for checking that the tools in the Hi-C pipeline are
functioning correctly.

The tests should be run in this order so that the required input files are
generated at the correct stage.

.. code-block:: none
   :linenos:

   pytest tests/test_gem_indexer.py
   pytest tests/test_tb_full_mapping.py
   pytest tests/test_tb_parse_mapping.py
   pytest tests/test_tb_filter.py
   pytest tests/test_tb_generate_tads.py
   pytest tests/test_tb_save_hdf5_matrix.py

These can be called as part of a single tool chain with:

.. code-block:: none
   :linenos:

   python tests/test_toolchains.py --pipeline hic