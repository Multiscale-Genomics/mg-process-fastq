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
   
   This pipeline can process FASTQ to identify protein-DNA binding sites.
   
   Running from the command line
   =============================
   
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
    .. code-block:: none
       :linenos:
       
       {
         "user" : "user_name",
         "submission_name" : "test_01",
         "expts" : [
           {
             "project_id": "ena_project_id",
             "group_name": "user_defined_name_1",
             "local": '',
             "run_ids" : ["ena_run_accn_1", "ena_run_accn_2", ...]
             "bgd_ids" : ["ena_run_accn_3", "ena_run_accn_4", ...]
           },
           {
             "project_id": "ena_project_id",
             "group_name": "user_defined_name_2",
             "run_ids" : ["ena_run_accn_5", "ena_run_accn_6", ...]
             "bgd_ids" : ["ena_run_accn_3", "ena_run_accn_4", ...]
           }
         ]
       }
      
   When using a local verion of the [COMPS virtual machine](http://www.bsc.es/computer-sciences/grid-computing/comp-superscalar/downloads-and-documentation):
   
   .. code-block:: none
      :linenos:
      
      runcompss --comm=integratedtoolkit.gat.master.GATAdaptor --log_level=debug --lang=python /home/compss/mg-process-fastq/process_chipseq.py --species homo_sapiens --assembly GRCh38 --project_id PRJDA34559 --run_ids PRJEB2445.json --data_dir /home/compss/VMShare
   
   Process Methods
   ---------------
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
   genome : str
       Genome accession (e.g. GCA_000001405.22)
   species : str
       Species (e.g. homo_sapiens)
   file : str
       Location of FASTQ input file
   
   Returns
   -------
   bed : file
       Bed files with the locations of nucleosome binding sites within the genome
   
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

