# mg-process-fastq
Scripts required for the processing of FASTQ files (eg generating adjacency lists for Hi-C data)

# Requirements
- Python 2.7.12 (required for MACS2 in ChIP-Seq pipeline)
- Python 3.5.2 (required for iNPS in MNase-Seq pipeline)
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




