Hi-C Test Data
==============

Test Data
---------

Dataset
^^^^^^^

+-----------+----------------+
| Stable ID | SRR1658573     |
| Citation  | Rao et al 2014 |
+-----------+----------------+

Genome
^^^^^^

+------------+----------+
| Assembly   | GRCh38   |
| Chromosome | 21       |
| Start      | 15000000 |
| End        | 19999999 |
+------------+----------+

Method
------
The full dataset was aligned to the genome using GEM. A list of FASTQ ids that
had aligned to the genome were then extracted:

.. code-block:: none
   :linenos:

   grep chr21 SRR1658573_1_full_1-100.map | tr "\t" "~" | cut -d"~" -f1 -f5 | tr ":" "\t" | awk '(NR==1) || (($4>15000000) && ($4<20000000))' | tr "\t" "~" | cut -d "~" -f1 > SRR1658573_1_chr21_1-100.row
   grep chr21 SRR1658573_2_full_1-100.map | tr "\t" "~" | cut -d"~" -f1 -f5 | tr ":" "\t" | awk '(NR==1) || (($4>15000000) && ($4<20000000))' | tr "\t" "~" | cut -d "~" -f1 > SRR1658573_2_chr21_1-100.row

The IDs were filtered to ensure matching pairs in both files:

.. code-block:: none
   :linenos:

   grep -Fx -f SRR1658573_1_chr21_1-100.row SRR1658573_2_chr21_1-100.row > SRR1658573_chr21.row

The `split_paired_fastq.py` was used to divide the original FASTQ files into
chunks of 1000000 reads. The `ExtractRowsFromFASTQ.py` script was then used to
extract the matching FASTQ pairs from each of the FASTQ files in parallel. All
of the individual FASTQ files were then concatenated together to form the final
2 FASTQ test files.


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