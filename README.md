# mg-process-fastq
Scripts required for the processing of FASTQ files (eg generating adjacency lists for Hi-C data)

# Requirements
- Python 2.7+
- Python Modules:
  - numpy
  - h5py
  - scipy
  - matplotlib
  - TADbit
  - pysam
- BS-Seeker2
  - Bowtie2
- imp (for 3D modelling with TADbit)
- mcl


# Whole Genome Bisulfite Sequencing (WGBS) Data Processing
## Code
### Filter the FastQ files

### Split the list
Reads through the paired input FastQ files, removed pairs that are only in one side then splits the original 2 FastQ files into pairs of 1000000 ready for alignment.

```
python split_paired_fastq.py --input_1=<file_1_paired.fastq> --input_2=<file_2_paired.fastq> --output_tag=<e.g. matching>
```

## COMPS Code
`process_wgbs.py` is a script that can get run on the COMPS infrastructure to convert the paired FastQ data for WGBS into the matching wig, ATCGmap and CGmap files.


# Hi-C Data Processing
Processing of paired end FastQ files from Hi-C experiments. Generates adjacency matrixes, computes TADs and generates the matching HDF5 files for using the REST API (mg-rest-hdf5). The mojority of the code has been wrapped up in the COMPS script, but this can be extracted and run locally. For the moment you need to comment out the `@task(...)` and `@constraint(...)` flags.

## Code - Generate Hi-C Adjacency Matrix
Takes an SRA file as the input and generates an adjacency matrix for the interacting sections of the genome.
```
python generate_adjacency_matrix.py --genome=<genome_accession> --dataset=<dataset_id> --sra_id=<sra_id> --library=<library_id> --enzyme_name=<enzyme_name> --resolution=<bin size> --tmp_dir=<tmp_dir> --data_dir=<download_directory>
```

## COMPS Code
`process_hic.py`

### Parameters:
* --genome
....Genome accession (e.g. GCA_000001405.22)
* --dataset
....Dataset ID from the reference paper (e.g. GSE63525)
* --expt_name
....A manual name to provide extra separation if there are multiple Hi-C sets within the same dataset
* --expt_list /home/compss/MuG_pipelines/mg-process-fastq/exptList.tsv
....A list of the experiments in a tab separated format including the following columns
....1. SRR ID
....2. Library name
....3. Retrictions enzyme (eg MboI)
* --tmp_dir \<temp_dir\>/
....This will be the location for the intermediary files. This gets shared between process, but they should not collide.
* --data_dir \<data_dir\>/
....This is where the initial FastQ files will be downloaded to and the output files will get saved.
