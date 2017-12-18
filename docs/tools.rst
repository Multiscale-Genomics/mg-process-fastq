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

Tools for processing FastQ files
================================

.. automodule:: tool

   File Validation
   ===============

   Pipelines and functions assessing the quality of input files.

   FastQC
   ------
   .. autoclass:: tool.validate_fastqc.fastqcTool
      :members:

   Indexers
   ========

   Bowtie 2
   --------
   .. autoclass:: tool.bowtie_indexer.bowtieIndexerTool
      :members:

   BS-Seeker2 Indexer
   ------------------
   .. autoclass:: tool.bs_seeker_indexer.bssIndexerTool
      :members:

   BWA
   ---
   .. autoclass:: tool.bwa_indexer.bwaIndexerTool
      :members:

   GEM
   ---
   .. autoclass:: tool.gem_indexer.gemIndexerTool
      :members:

   Kallisto
   --------
   .. autoclass:: tool.kallisto_indexer.kallistoIndexerTool
      :members:


   Aligners
   ========

   BWA
   ---
   .. autoclass:: tool.bwa_aligner.bwaAlignerTool
      :members:


   BS-Seeker2 Aligner
   ------------------
   .. autoclass:: tool.bs_seeker_aligner.bssAlignerTool
      :members:


   File Manipultation
   ==================
   .. autoclass:: tool.fastq_splitter.fastq_splitter
      :members:


   Filters
   =======

   BioBamBam Filter
   ----------------
   .. autoclass:: tool.biobambam_filter.biobambam
      :members:

   BS-Seeker2 Filter
   -----------------
   .. autoclass:: tool.bs_seeker_filter.filterReadsTool
      :members:


   Peak Calling
   ============

   BS-Seeker2 Methylation Caller
   -----------------------------
   .. autoclass:: tool.bs_seeker_methylation_caller.bssMethylationCallerTool
      :members:

   iNPS
   ----
   .. autoclass:: tool.inps.inps
      :members:

   Kallisto Quantification
   -----------------------
   .. autoclass:: tool.kallisto_quant.kallistoQuantificationTool
      :members:

   MACS2
   -----
   .. autoclass:: tool.macs2.macs2
      :members:


   Hi-C Parsing
   ============

   The following tools are a split out of the Hi-C pipelines generated to use
   the TADbit library.

   FASTQ mapping
   -------------
   .. autoclass:: tool.tb_full_mapping.tbFullMappingTool
      :members:

   Map Parsing
   -----------
   .. autoclass:: tool.tb_parse_mapping.tbParseMappingTool
      :members:

   Filter Aligned Reads
   --------------------
   .. autoclass:: tool.tb_filter.tbFilterTool
      :members:

   Save Matrix to HDF5 File
   ------------------------
   .. autoclass:: tool.tb_save_hdf5_matrix.tbSaveAdjacencyHDF5Tool
      :members:

   Generate TAD Predictions
   ------------------------
   .. autoclass:: tool.tb_generate_tads.tbGenerateTADsTool
      :members: