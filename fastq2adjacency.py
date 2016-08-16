import os, os.path, urllib2

import pytadbit
from pytadbit.mapping               import get_intersection
from pytadbit.mapping.filter        import apply_filter
from pytadbit.mapping.filter        import filter_reads
from pytadbit.mapping.mapper        import full_mapping
from pytadbit.parsers.map_parser    import parse_map
from pytadbit.parsers.hic_parser    import load_hic_data_from_reads
from pytadbit.parsers.hic_parser    import read_matrix
from pytadbit.parsers.genome_parser import parse_fasta
from pytadbit.utils.file_handling   import mkdir
import numpy as np


class fastq2adjacency:
    """
    These are the parts of the TADbit library that are required for processing
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
        self.genome_accession = '' # GCA_000001405.22 - GRChg38, current
        self.dataset     = '' # 'GSE63525'
        self.sra_id      = '' # 'SRR1658632'
        self.library     = '' # 'HiC036'
        self.enzyme_name = '' # 'NcoI'
        self.resolution  = 1000000

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

    def set_params(self, genome_accession, dataset, sra_id, library, enzyme_name, resolution, tmp_dir, data_dir, expt_name = None, same_fastq=True, windows1=None, windows2=None):
        self.genome_accession = genome_accession
        self.dataset     = dataset
        self.sra_id      = sra_id
        self.library     = library
        self.enzyme_name = enzyme_name
        self.resolution  = resolution

        self.temp_root = tmp_dir
        self.data_root = data_dir
        
        self.gem_file = self.data_root + self.genome_accession + "/" + self.genome_accession + ".gem"
        self.genome_file = self.data_root + self.genome_accession + "/chroms/" + self.genome_accession + ".fa"
        
        if expt_name != None:
            self.expt_name = expt_name + '/'
        else:
            self.expt_name = ''
        
        self.library_dir = self.data_root + self.expt_name + self.dataset + '/' + self.library + '/'
        if same_fastq == True:
            self.fastq_file_1  = self.library_dir + self.sra_id + '.fastq'
            self.fastq_file_2  = self.library_dir + self.sra_id + '.fastq'
        else:
            self.fastq_file_1  = self.library_dir + self.sra_id + '_1.fastq'
            self.fastq_file_2  = self.library_dir + self.sra_id + '_2.fastq'
        
        if expt_name != None:
            self.expt_name = expt_name + '/'
        else:
            self.expt_name = ''
        
        
        self.map_dir     = self.library_dir + '01_it-mapped_read'
        self.tmp_dir     = self.temp_root + self.expt_name + self.dataset + '/' + self.library
        self.parsed_reads_dir = self.tmp_dir + '/parsed_reads'
        
        try:
            os.makedirs(self.map_dir)
            os.makedirs(self.tmp_dir)
            os.makedirs(self.parsed_reads_dir)
        except:
            pass
        
        if windows1 != None:
            self.windows1 = windows1
        if windows2 != None:
            self.windows2 = windows2
    
    def getFastqData(self):
        f_index = urllib2.urlopen(
        'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=' + str(self.sra_id) + '&result=read_run&fields=study_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp&download=txt')
        data = f_index.read()
        rows = data.split("\n")
        row_count = 0
        for row in rows:
            if row_count == 0:
                row_count += 1
                continue
            
            row = row.rstrip()
            row = row.split("\t")
            
            if len(row) == 0:
                continue
            
            project = row[0]
            srr_ir = row[1]
            fastq_files = row[6].split(';')
            row_count += 1
            
            for fastq_file in fastq_files:
                file_name = fastq_file.split("/")
                print fastq_file
                print self.data_root + self.expt_name + self.dataset + '/' + self.library + '/' + file_name[-1]
                print file_name[-1]
                dl_file = open(self.data_root + self.expt_name + self.dataset + '/' + self.library + '/' + file_name[-1], "w")
                dl = urllib2.urlopen("ftp://" + fastq_file)
                dl_file.write(dl.read())
                dl_file.close()
        
    
    def mapWindows(self, side=1):
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
        Requires 8 CPU
        """
        # new file with info of each "read1" and its placement with respect to RE sites
        reads1 = self.parsed_reads_dir + '/read1.tsv'
        # new file with info of each "read2" and its placement with respect to RE sites
        reads2 = self.parsed_reads_dir + '/read2.tsv'
        
        mapped_rN = self.getMappedWindows()

        print 'Parse MAP files...'
        parse_map(mapped_rN["mapped_r1"], mapped_rN["mapped_r2"], out_file1=reads1, out_file2=reads2, genome_seq=self.genome_seq, re_name=self.enzyme_name, verbose=True, ncpus=8)
    
    def mergeMaps(self):
        """
        Merging mapped "read1" and "read2"
        """
        # Output file
        reads  = self.parsed_reads_dir + '/both_map.tsv'
        # new file with info of each "read1" and its placement with respect to RE sites
        reads1 = self.parsed_reads_dir + '/read1.tsv'
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
    
    def merge_adjacency_data(self, adjacency_matrixes):
        """
        The recommended route is to merge the normalised data is to sum the
        normalised values from previous steps. This should be the final step of
        the initial phase for main function should have finished by normalising
        each of the individual adjacency files.
        """
        merged_matrix = Chromosome(name=self.expt_name, centromere_search=True)
        merged_matrix.add_experiment(self.expt_name, resolution=self.resolution)
        merged_exp = merged_matrix.experiment[self.expt_name]
        print adjacency_data
        
        for m in adjacency_matrixes:
            new_chrom = Chromosome(name=self.expt_name, centromere_search=True)
            new_chrom.add_experiment(self.expt_name, hic_data=m, resolution=self.resolution)
            merged_exp = merged_exp + new_chrom.experiment[self.expt_name]
    
    def generate_tads(self):
        """
        Uses TADbit to generate the TAD borders based on the computed hic_data
        """
        my_chrom = Chromosome(name=self.expt_name, centromere_search=True)
        my_chrom.add_experiment(self.expt_name, hic_data=self.hic_data, resolution=self.resolution)
        
        # Run core TADbit function to find TADs on each expt.
        # For the current dataset required 61GB of RAM
        my_chrom.find_tad(self.expt_name, n_cpus=15)
        
        exp = my_chrom.experiments[self.expt_name]
        tad_file = self.library_dir + expt_name + '_tads_' + str(resolution) + '.tsv'
        tad_out = open(tad_file)
        tad_out.write(exp.write_tad_borders())
        tad_out.close()
    
    def load_hic_read_data(self):
        """
        Load the interactions into the HiC-Data data type
        """
        filter_reads = self.parsed_reads_dir + '/filtered_map.tsv'
        self.hic_data = load_hic_data_from_reads(filter_reads, resolution=self.resolution)
    
    def load_hic_matrix_data(self, norm=True):
        """
        Load the interactions from Hi-C adjacency matrix into the HiC-Data data
        type
        """
        if norm == True:
            # Dump the data pre-normalized
            adj_list = self.parsed_reads_dir + '/adjlist_map.tsv'
        else:
            adj_list = self.parsed_reads_dir + '/adjlist_map_norm.tsv'
        
        self.hic_data = read_matrix(adj_list, resolution=self.resolution)
    
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

