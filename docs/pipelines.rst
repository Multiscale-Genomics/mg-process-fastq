Pipelines
=========

Download and index genome files
-------------------------------
.. automodule:: process_genome
   
   Methods
   =======
   .. autoclass:: process_genome.process_genome
      :members:


Hi-C Anslysis
-------------
.. automodule:: process_hic
   
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
   When running the pipeline on a local machine:

   .. code-block:: none
      :linenos:

      python process_chipseq.py                                      \\
         --genome /<dataset_dir>/Homo_sapiens.GRCh38.fasta           \\
         --assembly GCA_000001405.25                                 \\
         --taxon_id 9606                                             \\
         --file /<dataset_dir>/<file_name>.fastq                     \\
         --bgd_file /<dataset_dir>/<bgd_file_name>.fastq

   When using a local verion of the [COMPS virtual machine](http://www.bsc.es/computer-sciences/grid-computing/comp-superscalar/downloads-and-documentation):
   
   .. code-block:: none
      :linenos:
      
      runcompss
         --comm=integratedtoolkit.gat.master.GATAdaptor                   \\
         --log_level=debug                                                \\
         --lang=python /home/compss/mg-process-fastq/process_chipseq.py   \\
            --taxon_id 9606                                               \\
            --assembly GCA_000001405.25                                   \\
            --file /<dataset_dir>/<file_name>.fastq                       \\
            --bgd_file /<dataset_dir>/<bgd_file_name>.fastq
            
   
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
   When running the pipeline on a local machine:

   .. code-block:: none
      :linenos:

      python process_mnaseseq.py                                     \\
         --genome /<dataset_dir>/Homo_sapiens.GRCh38.fasta           \\
         --assembly GCA_000001405.25                                 \\
         --taxon_id 9606                                             \\
         --file /<dataset_dir>/<file_name>.fastq

   When using a local verion of the [COMPS virtual machine](http://www.bsc.es/computer-sciences/grid-computing/comp-superscalar/downloads-and-documentation):
   
   .. code-block:: none
      :linenos:
      
      runcompss
         --comm=integratedtoolkit.gat.master.GATAdaptor                  \\
         --log_level=debug                                               \\
         --lang=python /home/compss/mg-process-fastq/process_mnaseseq.py \\
            --assembly GCA_000001635.2                                   \\
            --taxon_id 10090                                             \\
            --genome <data_dir>/GCA_000001635.2.fa                       \\
            --file /<dataset_dir>/<file_name>.fastq
   
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

      python process_rnaseq.py                                       \\
         --genome /<dataset_dir>/Homo_sapiens.GRCh38.cdna.all.fasta  \\
         --assembly GCA_000001405.25                                 \\
         --taxon_id 9606                                             \\
         --file /<dataset_dir>/expt1.fastq


   When using a local verion of the [COMPS virtual machine](http://www.bsc.es/computer-sciences/grid-computing/comp-superscalar/downloads-and-documentation):
   
   .. code-block:: none
      :linenos:
      
      runcompss
         --comm=integratedtoolkit.gat.master.GATAdaptor                   \\
         --log_level=debug                                                \\
         --lang=python /home/compss/mg-process-fastq/process_rnaseq.py    \\
            --assembly GCA_000001405.22                                   \\
            --taxon_id 9606                                               \\
            --genome <data_dir>/GCA_000001405.25.fa                       \\
            --file <data_dir>/expt1.fastq
   
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
   species : str
      Species (e.g. homo_sapien)
   fastq1 : str
      Location of FASTQ input file
   fastq2 : str
      Location of FASTQ input file
   aligner : str
      The aligner for BS-Seeker to use (bowtie2)
   aligner_path
      Location of the bowtie2 executable

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
   When using a local verion of the [COMPS virtual machine](http://www.bsc.es/computer-sciences/grid-computing/comp-superscalar/downloads-and-documentation):
   
   .. code-block:: none
      :linenos:
      
      runcompss
         --comm=integratedtoolkit.gat.master.GATAdaptor              \\
         --log_level=debug                                           \\
         --lang=python /home/compss/mg-process-fastq/process_wgbs.py \\
            --assembly GCA_000001405.22                              \\
            --taxon_id 9606                                          \\
            --genome <data_dir>/GCA_000001405.22.fa                  \\
            --fastq1 <data_dir>/expt1_a.fastq                        \\
            --fastq2 <data_dir>/expt1_b.fastq                        \\
            --aligner bowtie2                                        \\
            --aligner_path /home/compss/lib


   Methods
   =======
   .. autoclass:: process_wgbs.process_wgbs
      :members:

