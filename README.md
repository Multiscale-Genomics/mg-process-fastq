# mg-process-fastq
Scripts required for the processing of FASTQ files (eg generating adjacency lists for Hi-C data)

# Requirements
- Python 2.7.12+
- Python 3
- Python Modules:
  - numpy
  - h5py
  - scipy
  - matplotlib
  - TADbit
  - pysam
  - MACS2 - can be installed with pip, but runs on command line
  - rpy2
- GEMtools
- HDF5
- Bowtie2
- BWA
- SAMtools
- BS-Seeker2
- libmaus2
- BioBamBam2
- imp (for 3D modelling with TADbit)
- mcl
- R (2.9.1)
- iNPS


# Whole Genome Bisulfite Sequencing (WGBS) Data Processing
A set of scripts that can get run on the COMPS infrastructure to convert the paired FastQ data for WGBS into the matching wig, ATCGmap and CGmap files.

## Process
* FASTQ and assemblies are downloaded
* BSseeker2 filters the reads to remove low quality reads
* FASTQ files are split
* BSseeker2 uses BOWTIE2 to align the reads to the genomes
* Bam files are sorted and merged.
* BSseeker2 then calls the methylation regions

## Code
### Filter the FastQ files

### Split the list
Reads through the paired input FastQ files, removed pairs that are only in one side then splits the original 2 FastQ files into pairs of 1000000 ready for alignment.

```
python split_paired_fastq.py --input_1=<file_1_paired.fastq> --input_2=<file_2_paired.fastq> --output_tag=<e.g. matching>
```

## COMPS Code
`process_wgbs.py`

### Parameters:
* --genome

   Genome accession (e.g. GCA_000001405.22)
* --dataset

   Dataset ID from the reference paper (e.g. SRR1536575)
* --data_dir \<data_dir\>/

   This is where the initial FastQ files will be downloaded to and the output files will get saved.
* --tmp_dir \<temp_dir\>/

   This will be the location for the intermediary files. This gets shared between process, but they should not collide.

### Example
When using a local verion of the [COMPS virtual machine](http://www.bsc.es/computer-sciences/grid-computing/comp-superscalar/downloads-and-documentation):
```
runcompss --lang=python /home/compss/mg-process-fastq/process_wgbs.py --genome GCA_000001405.22 --dataset SRR1536575 --tmp_dir /home/compss/tmp/ --data_dir /home/compss/data/
```

# Hi-C Data Processing
Processing of paired end FastQ files from Hi-C experiments. Generates adjacency matrixes, computes TADs and generates the matching HDF5 files for using the REST API (mg-rest-hdf5). The mojority of the code has been wrapped up in the COMPS script, but this can be extracted and run locally. For the moment you need to comment out the `@task(...)` and `@constraint(...)` flags.

## Process
* FastQ and assembly is downloaded
* TADbit is used to parse the genome and FASTQ files
* TADbit generates the adjacency matrices
* TAD calls are made
* HDF5 files are generated

## Code - Generate Hi-C Adjacency Matrix
Takes an SRA file as the input and generates an adjacency matrix for the interacting sections of the genome.
```
python generate_adjacency_matrix.py --genome=<genome_accession> --dataset=<dataset_id> --sra_id=<sra_id> --library=<library_id> --enzyme_name=<enzyme_name> --resolution=<bin size> --tmp_dir=<tmp_dir> --data_dir=<download_directory>
```

## COMPS Code
`process_hic.py`

### Parameters:
* --genome

   Genome accession (e.g. GCA_000001405.22)
* --dataset

   Dataset ID from the reference paper (e.g. GSE63525)
* --expt_name

   A manual name to provide extra separation if there are multiple Hi-C sets within the same dataset
* --expt_list /home/compss/MuG_pipelines/mg-process-fastq/exptList.tsv

   A list of the experiments in a tab separated format including the following columns

   1. SRR ID
   2. Library name
   3. Retrictions enzyme (eg MboI)
* --tmp_dir \<temp_dir\>/

   This will be the location for the intermediary files. This gets shared between process, but they should not collide.
* --data_dir \<data_dir\>/

   This is where the initial FastQ files will be downloaded to and the output files will get saved.

### Output Files
* Adjacency list saved to an HDF5 file
* TSV file of the TAD regions

### Example
When using a local verion of the [COMPS virtual machine](http://www.bsc.es/computer-sciences/grid-computing/comp-superscalar/downloads-and-documentation):
```
runcompss --lang=python /home/compss/mg-process-fastq/process_hic.py --genome GCA_000001405.22 --dataset GSE63525 --expt_name rao2014 --expt_list /home/compss/mg-process-fastq/exptList.tsv --tmp_dir /home/compss/tmp/ --data_dir /home/compss/data/
```


# RNA-Seq Differential Expression Data Processing
A set of scripts that can get run on the COMPS infrastructure to process the paired FastQ data for RNA-Seq differenctial expression data.

## Process
* cDNAs and FastQ are downloaded
* Kallisto is used to map the reads to the cDNAs
* Kallisto then generates the quantification files

## COMPS Code
`process_rnaseq.py`

### Parameters:
* --species

   Species (e.g. homo_sapiens)
* --assembly

   Assembly accession (e.g. GRCh38)
* --project_id

   Project ID from the ENA website (e.g. PRJEB2445)
* --run_id

   Experiment Run Accession ID from the ENA website (e.g. ERR030872)
* --data_dir \<data_dir\>/

   This is where the initial FastQ files will be downloaded to and the output files will get saved.

### Output Files
* bed
* wig

### Example
When using a local verion of the [COMPS virtual machine](http://www.bsc.es/computer-sciences/grid-computing/comp-superscalar/downloads-and-documentation):
```
runcompss --lang=python /home/compss/mg-process-fastq/process_rnaseq.py --species homo_sapiens --assembly GRCh38 --project_id PRJEB2445 --run_id ERR030872 --data_dir /home/compss/data/
```


# CHiP-Seq Analysis
A set of scripts that can get run on the COMPS infrastructure to process the CHiP-Seq data and produce peak calls.

## Process
* BWA for the alignment
* BioBamBam2 for the filtering
* MACS2 for the peak calling

## COMPS Code
`process_chipseq.py`

### Parameters
* --species

   Species (e.g. homo_sapiens)
* --assembly

   Assembly accession (e.g. GRCh38)
* --project_id

   Project ID from the ENA website (e.g. PRJDA34559)
* --run_id

   Experiment Run Accession ID from the ENA website (e.g. DRR000150)
* --data_dir \<data_dir\>/

   This is where the initial FastQ files will be downloaded to and the output files will get saved.

### Example JSON file:
```
{
  "user" : "user_name",
  "submission_name" : "test_01",
  "expts" : [
    {
      "project_id": "ena_project_id",
      "group_name": "user_defined_name_1",
      "run_ids" : ["ena_run_accn_1", "ena_run_accn_2", ...]
      "bgd_ids" : ["ena_run_accn_3", "ena_run_accn_4", ...]
    },
    {
      "project_id": "ena_project_id",
      "group_name": "user_defined_name_2",
      "run_ids" : ["ena_run_accn_5", "ena_run_accn_6", ...]
      "bgd_ids" : ["ena_run_accn_3", "ena_run_accn_4", ...]
    },
    { ... }
  ]
}
```

### Output Files
* bed
* wig

## Example
When using a local verion of the [COMPS virtual machine](http://www.bsc.es/computer-sciences/grid-computing/comp-superscalar/downloads-and-documentation):
```
runcompss --comm=integratedtoolkit.gat.master.GATAdaptor --log_level=debug --lang=python /home/compss/mg-process-fastq/process_chipseq.py --species homo_sapiens --assembly GRCh38 --project_id PRJDA34559 --run_ids PRJEB2445.json --data_dir /home/compss/VMShare
```


# MNase-Seq Analysis
A set of scripts that can get run on the COMPS infrastructure to process the MNase-Seq data and produce peak calls. This uses the iNPS package for teh peak calling.

## COMPS Code
`process_mnaseseq.py`

### Parameters
* --species

   Species (e.g. homo_sapiens)
* --assembly

   Assembly accession (e.g. GRCh38)
* --project_id

   Project ID from the ENA website (e.g. PRJDA34559)
* --run_id

   Experiment Run Accession ID from the ENA website (e.g. DRR000150)
* --data_dir \<data_dir\>/

   This is where the initial FastQ files will be downloaded to and the output files will get saved.

### Example JSON file:
```
{
  "user" : "user_name",
  "submission_name" : "test_01",
  "expts" : [
    {
      # Paired-end dataset
      "project_id": "ena_project_id",
      "group_name": "user_defined_name_1",
      "run_ids" : ["ena_run_accn_1", "ena_run_accn_2"]
    },
    {
      # Single ended dataset
      "project_id": "ena_project_id",
      "group_name": "user_defined_name_2",
      "run_ids" : ["ena_run_accn_5"]
    },
    { ... }
  ]
}
```

## Example
When using a local verion of the [COMPS virtual machine](http://www.bsc.es/computer-sciences/grid-computing/comp-superscalar/downloads-and-documentation):
```
runcompss --lang=python /home/compss/mg-process-fastq/process_mnaseseq.py --species homo_sapiens --assembly GRCh38 --project_id PRJEB2445 --run_id ERR030872 --data_dir /home/compss/data/
```


