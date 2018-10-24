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

Workflows
=========

Workflows are the scripts that define the pipelines and launc each of the tools

.. automodule:: mg_process_fastq.workflow

   Multitool Workflows
   ===================

   ChIP-seq
   --------
   Used for the analysis of FastQ files. Handles paired or single end FastQ inputs.
   It uses the BWA MEM alginer, filters out repeats and technical artifacts using
   BioBamBam2 and the peak calling uses MACS2.

   Methods
   ^^^^^^^

   .. autoclass:: mg_process_fastq.workflow.chipseq.process_chipseq
      :members:

   Genome Indexer
   --------------
   Used to generate assembly indexes for Bowtie2, BWA and GEM. This does not include
   the index for the WGBS Bowtie2 indexes, these are calculated as they are required.

   Methods
   ^^^^^^^

   .. autoclass:: mg_process_fastq.workflow.genome_indexer.process_genome
      :members:


   MNase-Seq
   ---------
   This is a workflow for the identification of nucleosome binding sites. It uses
   BWA ALN for aligning single or paired end FASTQ data. BioBamBam2 is used for the
   filtering of duplicate data and technical artifacts. iNPS is used for the
   peak calling of binding sites.

   Methods
   ^^^^^^^

   .. autoclass:: mg_process_fastq.workflow.mnaseseq.process_mnaseseq
      :members:

   RNA-seq
   -------
   This workflow uses Kallisto to index and quantify the level of gene expression.

   Methods
   ^^^^^^^

   .. autoclass:: mg_process_fastq.workflow.rnaseq.process_rnaseq
      :members:


   Whole Genome Bisulfate Sequencing (WGBS)
   ----------------------------------------
   Use for the analysis of Whole genome Bisulfate sequencing data to identify
   protein-DNA bringing sites. The workflow uses the BS Seeker2 suit of tools for
   filtering of the FASTQ data, indexing and alignment to the assembly and then
   the peak calling.

   Methods
   ^^^^^^^

   .. autoclass:: mg_process_fastq.workflow.wgbs.process_wgbs
      :members:

   Single Tool Workflows
   =====================

   Bowtie2
   -------

   Methods
   ^^^^^^^

   .. autoclass:: mg_process_fastq.workflow.bowtie2.process_bowtie
      :members:


   BWA ALN
   -------

   Methods
   ^^^^^^^

   .. autoclass:: mg_process_fastq.workflow.bwa_aln.process_bwa_aln
      :members:


   BWA MEM
   -------

   Methods
   ^^^^^^^

   .. autoclass:: mg_process_fastq.workflow.bwa_mem.process_bwa_mem
      :members:


