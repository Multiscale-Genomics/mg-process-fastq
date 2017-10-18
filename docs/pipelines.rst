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

Pipelines
=========

Download and index genome files
-------------------------------

.. automodule:: process_genome

   This pipeline is for the indexing of genomes once they have been loaded into
   the VRE. It indexes each new genome with Bowtie2, BWA and GEM. These indexes
   can then be used by the other pipelines.

   Running from the command line
   =============================

   Parameters
   ----------
   taxon_id : int
      Species taxonomic ID
   assembly : str
      Genomic assembly ID
   genome : str
      Location of the genomes FASTA file

   Returns
   -------
   Bowtie2 index files
   BWA index files
   GEM index file

   Example
   -------
   When running the pipeline on a local machine

   .. code-block:: none
      :linenos:

      python process_genome.py                              \
          --taxon_id 9606                                   \
          --assembly GRCh38                                 \
          --genome /<dataset_dir>/Homo_sapiens.GRCh38.fasta

   When using a local version of the [COMPS virtual machine](https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/):

   .. code-block:: none
      :linenos:

      runcompss                               \
          --lang=python                       \
          --library_path=${HOME}/bin          \
          --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \
          --log_level=debug process_genome.py \
          --taxon_id 9606                     \
          --genome /<dataset_dir>/tb.Human.GCA_000001405.22.fasta \
          --assembly GRCh38


   Methods
   =======

   .. autoclass:: process_genome.process_genome
      :members:


Hi-C Analysis
-------------

.. automodule:: process_hic

   This pipeline can process paired end FASTQ files to identify structural
   interactions that occur so that the genome can fold into the confines of the
   nucleus

   Running from the command line
   =============================

   Parameters
   ----------
   genome : str
      Location of the genomes FASTA file
   genome_gem : str
      Location of the genome GEM index file
   taxon_id : int
      Species taxonomic ID
   assembly : str
      Genomic assembly ID
   file1 : str
      Location of FASTQ file 1
   file2 : str
      Location of FASTQ file 2
   resolutions : str
      Comma separated list of resolutions to calculate the matrix for.
      [DEFAULT : 1000000,10000000]
   enzyme_name : str
      Name of the enzyme used to digest the genome (example 'MboI')
   window_type : str
      iter | frag. Analysis windowing type to use
   windows1 : str
      FASTQ sampling window sizes to use for the first paired end FASTQ file,
      the default is to use `[[1,25], [1,50], [1,75], [1,100]]`. This would be
      represented as `1,25,50,75,100` as input for this variable
   windows2 : str
      FASTQ sampling window sizes to use for the second paired end FASTQ file,
      the default is to use `[[1,25], [1,50], [1,75], [1,100]]`. This would be
      represented as `1,25,50,75,100` as input for this variable
   normalized : int
      1 | 0. Determines whether the counts of alignments
      should be normalized
   tag : str
      Name for the experiment output files to use

   Returns
   -------
   Adjacency List : file
   HDF5 Adjacency Array : file

   Example
   -------
   REQUIREMENT - Needs the indexing step to be run first

   When running the pipeline on a local machine:

   .. code-block:: none
      :linenos:

      python process_hic.py                                   \
         --genome /<dataset_dir>/Homo_sapiens.GRCh38.fasta    \
         --genome_gem /<dataset_dir>/Homo_sapiens.GRCh38.gem  \
         --assembly GCA_000001405.25                          \
         --taxon_id 9606                                      \
         --file1 /<dataset_dir>/<file_name>_1.fastq           \
         --file2 /<dataset_dir>/<file_name>_2.fastq           \
         --resolutions 1000000,10000000                       \
         --enzyme_name MboI                                   \
         --windows1 1,100                                     \
         --windows2 1,100                                     \
         --normalized 1                                       \
         --tag Human.SRR1658573                            \
         --window_type frag

   When using a local version of the [COMPS virtual machine](https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/):

   .. code-block:: none
      :linenos:

      runcompss                        \
         --lang=python                 \
         --library_path=${HOME}/bin    \
         --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \
         --log_level=debug             \
         process_hic.py                \
            --taxon_id 9606            \
            --genome /<dataset_dir>/.Human.GCA_000001405.22_gem.fasta \
            --assembly GRCh38          \
            --file1 /<dataset_dir>/Human.SRR1658573_1.fastq \
            --file2 /<dataset_dir>/Human.SRR1658573_2.fastq \
            --genome_gem /<dataset_dir>/Human.GCA_000001405.22_gem.fasta.gem \
            --enzyme_name MboI         \
            --resolutions 10000,100000 \
            --windows1 1,100           \
            --windows2 1,100           \
            --normalized 1             \
            --tag Human.SRR1658573     \
            --window_type frag

   Methods
   =======
   .. autoclass:: process_hic.process_hic
      :members:


ChIP-Seq Analysis
-----------------
.. automodule:: process_chipseq

   This pipeline can process FASTQ to identify protein-DNA binding sites.

   Running from the command line
   =============================

   Parameters
   ----------
   genome : str
      Genome accession (e.g. GCA_000001405.22)
   taxon_id : int
      Species (e.g. 9606)
   file : str
      Location of FASTQ input file
   bgd_file : str
      Location of FASTQ background file

   Returns
   -------
   bed : file
      Bed files with the locations of transcription factor binding sites
      within the genome

   Example
   -------
   REQUIREMENT - Needs the indexing step to be run first

   When running the pipeline on a local machine:

   .. code-block:: none
      :linenos:

      python process_chipseq.py                            \
         --genome /<dataset_dir>/Homo_sapiens.GRCh38.fasta \
         --assembly GCA_000001405.25                       \
         --taxon_id 9606                                   \
         --file /<dataset_dir>/<file_name>.fastq           \
         --bgd_file /<dataset_dir>/<bgd_file_name>.fastq

   When using a local version of the [COMPS virtual machine](https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/):

   .. code-block:: none
      :linenos:

      runcompss                     \
         --lang=python              \
         --library_path=${HOME}/bin \
         --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \
         --log_level=debug          \
         process_chipseq.py         \
            --taxon_id 9606         \
            --genome /<dataset_dir>/Human.GCA_000001405.22.fasta \
            --assembly GRCh38       \
            --file /<dataset_dir>/DRR000150.22.fastq


   Methods
   =======
   .. autoclass:: process_chipseq.process_chipseq
      :members:


Mnase-Seq Analysis
------------------
.. automodule:: process_mnaseseq

   This pipeline can process FASTQ to identify nucleosome binding sites.

   Running from the command line
   =============================

   Parameters
   ----------
   assembly : str
      Genome assembly ID (e.g. GCA_000001635.2)
   taxon_id : int
      Taxon ID (e.g. 10090)
   genome : str
      Location of genome assembly FASTA file
   species : str
      Species (e.g. mus_musculus)
   file : str
      Location of FASTQ input file

   Returns
   -------
   bed : file
      Bed files with the locations of nucleosome binding sites within the genome

   Example
   -------
   REQUIREMENT - Needs the indexing step to be run first

   When running the pipeline on a local machine:

   .. code-block:: none
      :linenos:

      python process_mnaseseq.py                           \
         --genome /<dataset_dir>/Homo_sapiens.GRCh38.fasta \
         --assembly GCA_000001405.25                       \
         --taxon_id 9606                                   \
         --file /<dataset_dir>/<file_name>.fastq

   When using a local version of the [COMPS virtual machine](https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/):

   .. code-block:: none
      :linenos:

      runcompss                                        \
         --lang=python                                 \
         --library_path=${HOME}/bin                    \
         --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \
         --log_level=debug                             \
         process_mnaseseq.py                           \
            --taxon_id 10090                           \
            --genome /<dataset_dir>/Mouse.GRCm38.fasta \
            --assembly GRCm38                          \
            --file /<dataset_dir>/DRR000386.fastq


   Methods
   =======
   .. autoclass:: process_mnaseseq.process_mnaseseq
      :members:


RNA-Seq Analysis
----------------
.. automodule:: process_rnaseq

   This pipeline can process FASTQ to quantify the level of expression of cDNAs.

   Running from the command line
   =============================

   Parameters
   ----------
   assembly : str
      Genome assembly ID (e.g. GCA_000001405.22)
   taxon_id : int
      Taxon ID (e.g. 9606)
   genome : str
      Location of genome assembly FASTA file
   file : str
      Location of FASTQ input file
   file2 : str
      [OPTIONAL] Location if FASTQ file if the data is paired end

   Returns
   -------
   bed : file
      WIG file with the levels of expression for genes

   Example
   -------
   When running the pipeline on a local machine:

   .. code-block:: none
      :linenos:

      python process_rnaseq.py                                       \
         --genome /<dataset_dir>/Homo_sapiens.GRCh38.cdna.all.fasta  \
         --assembly GCA_000001405.25                                 \
         --taxon_id 9606                                             \
         --file /<dataset_dir>/expt_1.fastq                          \
         --file2 /<dataset_dir>/expt_2.fastq


   When using a local version of the [COMPS virtual machine](https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/):

   .. code-block:: none
      :linenos:

      runcompss                                        \
         --lang=python                                 \
         --library_path=${HOME}/bin                    \
         --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \
         --log_level=debug                             \
         process_rnaseq.py                             \
            --taxon_id 9606                            \
            --genome /<dataset_dir>/Human.GRCh38.fasta \
            --assembly GRCh38                          \
            --file /<dataset_dir>/ERR030872_1.fastq    \
            --file2 /<dataset_dir>/ERR030872_2.fastq

   Methods
   =======
   .. autoclass:: process_rnaseq.process_rnaseq
      :members:


Whole Genome BiSulphate Sequencing Analysis
-------------------------------------------
.. automodule:: process_wgbs

   This pipeline can process FASTQ to identify methylation sites.

   Running from the command line
   =============================

   Parameters
   ----------
   assembly : str
      Genome assembly ID (e.g. GCA_000001405.22)
   taxon_id : int
      Taxon ID (e.g. 9606)
   genome : str
      Location of genome assembly FASTA file
   fastq1 : str
      Location of FASTQ input file
   fastq2 : str
      Location of FASTQ input file
   aligner : str
      The aligner for BS-Seeker to use (bowtie2)
   aligner_path
      Location of the bowtie2 executable
   bss_path : str
      Location of the BS-Seeker2 scripts

   Returns
   -------
   wig_file : file
      Location of a wig file containing the methylation peak calling results
   cgmap_file : file
      Location of the CGmap file from BS-Seeker2
   atcgmap_file : file
      Location of the ATCGmap file from BS-Seeker2

   A full description of the BS-Seeker2 files can be found at
   https://github.com/BSSeeker/BSseeker2

   Example
   -------
   When running the pipeline on a local machine:

   .. code-block:: none
      :linenos:

      python process_rnaseq.py                   \
         --assembly GCA_000001405.22             \
         --taxon_id 9606                         \
         --genome <data_dir>/GCA_000001405.22.fa \
         --fastq1 <data_dir>/expt1_a.fastq       \
         --fastq2 <data_dir>/expt1_b.fastq       \
         --aligner bowtie2                       \
         --aligner_path /home/compss/lib         \
         --bss_path <script_dir>/BS-Seeker2

   When using a local version of the [COMPS virtual machine](https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/):

   .. code-block:: none
      :linenos:

      runcompss                                  \
         --lang=python                           \
         --library_path=${HOME}/bin              \
         --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \
         --log_level=debug                       \
         process_wgbs.py                         \
            --taxon_id 10090                     \
            --genome /<dataset_dir>/Mouse.GRCm38.fasta \
            --assembly GRCm38                    \
            --fastq1 /<dataset_dir>/expt_1.fastq \
            --fastq2 /<dataset_dir>/expt_2.fastq \
            --aligner bowtie2                    \
            --aligner_path ${HOME}/lib/bowtie2-2.3.2 \
            --bss_path ${HOME}/lib/BSseeker2


   Methods
   =======
   .. autoclass:: process_wgbs.process_wgbs
      :members:

