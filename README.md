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


# Whole Genome Bisulfite Sequencing (WGBS)
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
## Code - Generate Hi-C Adjacency Matrix
Takes an SRA file as the input and generates an adjacency matrix for the interacting sections of the genome.
```
python generate_adjacency_matrix.py --genome=<genome_accession> --dataset=<dataset_id> --sra_id=<sra_id> --library=<library_id> --enzyme_name=<enzyme_name> --resolution=<bin size> --tmp_dir=<tmp_dir> --data_dir=<download_directory>
```

## Code - Generate TAD predictions

## COMPS Code
This script will act to proceess the SRA file all the way through to generating the adjacency matrix

