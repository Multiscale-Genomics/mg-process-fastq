import os, os.path

import pytadbit
from pytadbit.mapping               import get_intersection
from pytadbit.mapping.filter        import apply_filter
from pytadbit.mapping.filter        import filter_reads
from pytadbit.mapping.mapper        import full_mapping
from pytadbit.parsers.map_parser    import parse_map
from pytadbit.parsers.hic_parser    import load_hic_data_from_reads
from pytadbit.parsers.genome_parser import parse_fasta
from pytadbit.utils.file_handling   import mkdir
import numpy as np


class fastq2adjacency:
    """
    These are the parts of teh TADbit library that are required for processing
    FASTQ data. They have been packaged into chunks that can be easily handled
    by the COMPS infrastructure.
    
    At the moment this assumes 1 SRA generates a single adjacency matrix. Often
    there are multiple SRA files that get merged for a final result. This needs
    to be integrated into this pipeline.
    """
    
    def __init__(self):
        """
        Initialise the module and 
        """
        self.genome_accession = 'GCA_000001405.22' # GRChg38, current
        self.dataset     = 'GSE63525'
        self.sra_id      = '' # 'SRR1658632'
        self.library     = '' # 'HiC036'
        self.enzyme_name = '' # 'NcoI'
        self.resolution  = 100000

        self.temp_root = '/' # '/<tmp_area>/'
        self.data_root = '/' # '/<data_dir>/'
        
        self.gem_file     = ''
        self.genome_file  = ''
        self.fastq_file_1  = ''
        self.fastq_file_2  = ''
        self.map_dir      = ''
        self.tmp_dir      = ''
        
        self.windows1 = ((1,25), (1,50), (1,75),(1,100))
        self.windows2 = ((101,125), (101,150), (101,175),(101,200))
        
        self.mapped_r1 = None
        self.mapped_r2 = None
        
        self.genome_seq = None
        
        self.hic_data = None

    def set_params(self, genome_accession, dataset, sra_id, library, enzyme_name, resolution, tmp_dir, data_dir, same_fastq=True, windows1=None, windows2=None):
        self.genome_accession = genome_accession
        self.dataset     = dataset
        self.sra_id      = sra_id
        self.library     = library
        self.enzyme_name = enzyme_name
        self.resolution  = resolution

        self.temp_root = tmp_dir
        self.data_root = data_dir
    
        self.gem_file    = self.data_root + self.genome_accession + '/' + self.genome_accession + '.gem'
        self.genome_file = self.data_root + self.genome_accession + '/chromos/' + self.genome_accession + '.fa'
        
        self.fastq_file_1  = self.data_root + self.dataset + '/' + self.library + '/'
        self.fastq_file_2  = self.data_root + self.dataset + '/' + self.library + '/'
        if same_fastq == True:
            self.fastq_file_1  = self.fastq_file_1 + self.sra_id + '.fastq'
            self.fastq_file_2  = self.fastq_file_1 + self.sra_id + '.fastq'
        else:
            self.fastq_file_1  = self.fastq_file_1 + self.sra_id + '_1.fastq'
            self.fastq_file_2  = self.fastq_file_1 + self.sra_id + '_2.fastq'
            
        self.map_dir     = self.data_root + self.dataset + '/' + self.library + '/01_it-mapped_read'
        self.tmp_dir     = self.temp_root + self.dataset + '/' + self.library
        self.parsed_reads_dir = tmp_dir + '/parsed_reads'
        
        mkdir(self.parsed_reads_dir)
        
        if windows1 != None:
            self.windows1 = windows1
        if windows2 != None:
            self.windows2 = windows2
    
    def mapWindows(self, side=1, same_file=True):
        """
        Map the reads to the genome
        """
        
        if side == 1:
            mapped_r1 = full_mapping(self.gem_file, self.fastq_file_1, self.map_dir + str(1), windows=self.windows1, frag_map=False, nthreads=8, clean=True, temp_dir=self.tmp_dir)
        elif side == 2:
            mapped_r2 = full_mapping(self.gem_file, self.fastq_file_2, self.map_dir + str(2), windows=self.windows2, frag_map=False, nthreads=8, clean=True, temp_dir=self.tmp_dir)
    
    def getMappedWindows(self):
        """
        Populate the mapped_rN values so that it is not reliant on a single
        process
        """
        
        mapped_r1 = []
        mapped_r2 = []

        r1_dir = self.map_dir + str(1)
        r2_dir = self.map_dir + str(2)

        for mapped in os.listdir(r1_dir):
            mapped_r1.append(os.path.join(r1_dir, mapped))
        for mapped in os.listdir(r2_dir):
            mapped_r2.append(os.path.join(r2_dir, mapped))
        
        return {'mapped_r1': mapped_r1, 'mapped_r2': mapped_r2}
    
    def parseGenomeSeq(self):
        """
        Loads the genome
        """
        self.genome_seq = parse_fasta(self.genome_file)
    
    def parseMaps(self):
        """
        Merge the 2 read maps together 
        """
        # new file with info of each "read1" and its placement with respect to RE sites
        reads1 = self.parsed_reads_dir + '/read1.tsv'
        # new file with info of each "read2" and its placement with respect to RE sites
        reads2 = self.parsed_reads_dir + '/read2.tsv'
        
        mapped_rN = getMappedWindows()

        print 'Parse MAP files...'
        parse_map(mapped_rN["mapped_r1"], mapped_rN["mapped_r2"], out_file1=reads1, out_file2=reads2, genome_seq=self.genome_seq, re_name=self.enzyme_name, verbose=True, ncpus=8)
    
    def mergeMaps(self):
        """
        Merging mapped "read1" and "read2"
        """
        # Output file
        reads  = parsed_reads_dir + '/both_map.tsv'
        # new file with info of each "read1" and its placement with respect to RE sites
        reads1 = parsed_reads_dir + '/read1.tsv'
        # new file with info of each "read2" and its placement with respect to RE sites
        reads2 = self.parsed_reads_dir + '/read2.tsv'
        get_intersection(reads1, reads2, reads, verbose=True)
    
    def filterReads(self, conservative = True):
        """
        Filter the reads to remove duplicates and experimental abnormalities
        Requires 4 CPU
        """
        
        reads      = self.parsed_reads_dir + '/both_map.tsv'
        filt_reads = self.parsed_reads_dir + '/filtered_map.tsv'
        
        masked = filter_reads(reads, max_molecule_length=610, min_dist_to_re=915, over_represented=0.005, max_frag_size=100000, min_frag_size=100, re_proximity=4)

        if conservative == True:
            # Ignore filter 5 (based on docs) as not very helpful
            apply_filter(reads, filt_reads, masked, filters=[1,2,3,4,6,7,8,9,10])
        else:
            # Less conservative option
            apply_filter(reads, filt_reads, masked, filters=[1,2,3,9,10])
    
    def load_hic_data(self):
        """
        Load the interactions into the HiC-Data data type
        """
        filt_reads = self.parsed_reads_dir + '/filtered_map.tsv'
        self.hic_data = load_hic_data_from_reads(filt_reads, resolution=resolution)
    
    def normalise_hic_data(self, iterations=0):
        """
        Normalise the Hi-C data
        Example has the iterations set to 9, but setting to 0 to match that
        done by Rao et al 2014
        """
        self.hic_data.normalize_hic(iterations=iterations, max_dev=0.1)
    
    def save_hic_data(self, normalized=False):
        """
        Save the hic_data object to a file. This is saved as an NxN array with
        the values for all positions being set.
        """
        if normalized == False:
            # Dump the data pre-normalized
            adj_list = self.parsed_reads_dir + '/adjlist_map.tsv'
            self.hic_data.write_matrix(adj_list, normalized=False)
        else:
            adj_list = self.parsed_reads_dir + '/adjlist_map_norm.tsv'
            self.hic_data.write_matrix(adj_list, normalized=True)

