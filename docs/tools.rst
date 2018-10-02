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

.. automodule:: mg_process_fastq.tool

   File Validation
   ===============

   Pipelines and functions assessing the quality of input files.

   FastQC
   ------
   .. autoclass:: mg_process_fastq.tool.validate_fastqc.fastqcTool
      :members:

   TrimGalore
   ----------
   .. autoclass:: mg_process_fastq.tool.trimgalore.trimgalore
      :members:

   Indexers
   ========

   Bowtie 2
   --------
   .. autoclass:: mg_process_fastq.tool.bowtie_indexer.bowtieIndexerTool
      :members:

   BSgenome Index
   --------------
   .. autoclass:: mg_process_fastq.tool.forge_bsgenome.bsgenomeTool
      :members:

   BS-Seeker2 Indexer
   ------------------
   .. autoclass:: mg_process_fastq.tool.bs_seeker_indexer.bssIndexerTool
      :members:

   BWA
   ---
   .. autoclass:: mg_process_fastq.tool.bwa_indexer.bwaIndexerTool
      :members:

   GEM
   ---
   .. autoclass:: mg_process_fastq.tool.gem_indexer.gemIndexerTool
      :members:

   Kallisto
   --------
   .. autoclass:: mg_process_fastq.tool.kallisto_indexer.kallistoIndexerTool
      :members:


   Aligners
   ========

   Bowtie2
   -------
   .. autoclass:: mg_process_fastq.tool.bowtie_aligner.bowtie2AlignerTool
      :members:

   BWA - ALN
   ---------
   .. autoclass:: mg_process_fastq.tool.bwa_aligner.bwaAlignerTool
      :members:

   BWA - MEM
   ---------
   .. autoclass:: mg_process_fastq.tool.bwa_mem_aligner.bwaAlignerMEMTool
      :members:

   BS-Seeker2 Aligner
   ------------------
   .. autoclass:: mg_process_fastq.tool.bs_seeker_aligner.bssAlignerTool
      :members:


   Filters
   =======

   BioBamBam Filter
   ----------------
   .. autoclass:: mg_process_fastq.tool.biobambam_filter.biobambam
      :members:

   BS-Seeker2 Filter
   -----------------
   .. autoclass:: mg_process_fastq.tool.bs_seeker_filter.filterReadsTool
      :members:

   Trim Galore
   -----------
   .. autoclass:: mg_process_fastq.tool.trimgalore.trimgalore
      :members:


   Peak Calling
   ============

   BS-Seeker2 Methylation Caller
   -----------------------------
   .. autoclass:: mg_process_fastq.tool.bs_seeker_methylation_caller.bssMethylationCallerTool
      :members:

   iDEAR
   -----
   .. autoclass:: mg_process_fastq.tool.idear.idearTool
      :members:

   iNPS
   ----
   .. autoclass:: mg_process_fastq.tool.inps.inps
      :members:

   Kallisto Quantification
   -----------------------
   .. autoclass:: mg_process_fastq.tool.kallisto_quant.kallistoQuantificationTool
      :members:

   MACS2
   -----
   .. autoclass:: mg_process_fastq.tool.macs2.macs2
      :members:


   Hi-C Parsing
   ============

   The following tools are a split out of the Hi-C pipelines generated to use
   the TADbit library.

   FASTQ mapping
   -------------
   .. autoclass:: mg_process_fastq.tool.tb_full_mapping.tbFullMappingTool
      :members:

   Map Parsing
   -----------
   .. autoclass:: mg_process_fastq.tool.tb_parse_mapping.tbParseMappingTool
      :members:

   Filter Aligned Reads
   --------------------
   .. autoclass:: mg_process_fastq.tool.tb_filter.tbFilterTool
      :members:

   Identify TADs and Compartments
   ------------------------------
   .. autoclass:: mg_process_fastq.tool.tb_segment.tbSegmentTool
      :members:

   Normalize paired end reads file
   -------------------------------
   .. autoclass:: mg_process_fastq.tool.tb_normalize.tbNormalizeTool
      :members:

   Extract binned matrix from paired end reads file
   ------------------------------------------------
   .. autoclass:: mg_process_fastq.tool.tb_bin.tbBinTool
      :members:

   Save Matrix to HDF5 File
   ------------------------
   .. autoclass:: mg_process_fastq.tool.tb_save_hdf5_matrix.tbSaveAdjacencyHDF5Tool
      :members:

   Generate TAD Predictions
   ------------------------
   .. autoclass:: mg_process_fastq.tool.tb_generate_tads.tbGenerateTADsTool
      :members: