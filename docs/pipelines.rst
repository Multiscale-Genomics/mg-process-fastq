Pipelines
=========

Download and index genome files
-------------------------------
.. automodule:: process_genome
   
   Process Methods
   ---------------
   .. autoclass:: process_genome.process_genome
      :members:


Hi-C Anslysis
-------------
.. automodule:: process_hic
   
   Process Methods
   ---------------
   .. autoclass:: process_hic.process_hic
      :members:


ChIP-Seq Analysis
-----------------
.. automodule:: process_chipseq
   
   Running the pipeline from the command line:
   
   Parameters
   ----------
   genome : str
       Genome accession (e.g. GCA_000001405.22)
   species : str
       Species (e.g. homo_sapiens)
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
   When using a local verion of the [COMPS virtual machine](http://www.bsc.es/computer-sciences/grid-computing/comp-superscalar/downloads-and-documentation):
   
   :: console
      runcompss --lang=python /home/compss/mg-process-fastq/process_hic.py --genome GCA_000001405.22 --dataset GSE63525 --expt_name rao2014 --expt_list /home/compss/mg-process-fastq/exptList.tsv --tmp_dir /home/compss/tmp/ --data_dir /home/compss/data/
   
   Process Methods
   ---------------
   .. autoclass:: process_chipseq.process_chipseq
      :members:


Mnase-Seq Analysis
------------------
.. automodule:: process_mnaseseq
   
   Process Methods
   ---------------
   .. autoclass:: process_mnaseseq.process_mnaseseq
      :members:


RNA-Seq Analysis
----------------
.. automodule:: process_rnaseq
   
   Process Methods
   ---------------
   .. autoclass:: process_rnaseq.process_rnaseq
      :members:


Whole Genome BiSulphate Sequencing Analysis
-------------------------------------------
#.. automodule:: process_wgbs
#   
#   Process Methods
#   ---------------
#   .. autoclass:: process_wgbs.process_wgbs
#      :members:

