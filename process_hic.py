#!/usr/bin/python

"""
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
"""

# -*- coding: utf-8 -*-
"""process Hi-C paired end FastQ files"""
import argparse, urllib2, gzip, shutil, shlex, subprocess, os, sys, json, time

from basic_modules import Tool, Workflow, Metadata
from dmp import dmp

from functools import wraps

import tool

try :
    from pycompss.api.parameter import *
    from pycompss.api.task import task
except ImportError :
    print "[Warning] Cannot import \"pycompss\" API packages."
    print "          Using mock decorators."
    
    from dummy_pycompss import *

from tool import tb_full_mapping
from tool import tb_parse_mapping
from tool import tb_filter

# ------------------------------------------------------------------------------
class process_hic(Workflow):
    """
    Functions for downloading and processing Mnase-seq FastQ files. Files are
    downloaded from the European Nucleotide Archive (ENA), then aligned,
    filtered and analysed for peak calling
    """
    
    def run(self, file_ids, metadata):
        """
        Main run function for processing MNase-Seq FastQ data. Pipeline aligns
        the FASTQ files to the genome using BWA. iNPS is then used for peak
        calling to identify nucleosome position sites within the genome.
                
        Parameters
        ----------
        files_ids : list
            List of file locations
        metadata : list
            Required meta data
        
        Returns
        -------
        outputfiles : list
            List of locations for the output bam, bed and tsv files
        """
        
        from fastq2adjacency import fastq2adjacency

        genome_file  = files_ids[0]
        genome_gem   = files_ids[1]
        fastq_file_1 = files_ids[2]
        fastq_file_2 = files_ids[3]
        enzyme_name = metadata['enzyme_name']
        resolutions = metadata['resolutions']
        windows1    = metadata['windows1']
        windows2    = metadata['windows2']

        tmp_name = fastq_file_1.split('/')
        tmp_dir = '/'.join(tmp_name[0:-1])
        try:
            os.makedirs(tmp_dir)
        except:
            pass

        tfm1 = tb_full_mapping.tbFullMappingTool()
        tfm1_files, tfm1_meta = tfm1.run([genome_gem, fastq_file_1], {'windows' : windows1})

        tfm2 = tb_full_mapping.tbFullMappingTool()
        tfm2_files, tfm2_meta = tfm2.run([genome_gem, fastq_file_2], {'windows' : windows2})

        tpm = tb_parse_mapping.tb_parse_mapping()
        files = [genome_file] + tfm1_files + tfm2_files
        metadata = {'enzyme_name' : enzyme_name}
        tpm_files, tpm_meta = tpm.run(files, metadata)

        tf = tb_filter.tbFilterTool()
        tf_files, tf_meta = tf.run(tpm_files, {'conservative' : True})

        #adjlist_loc = f2a.save_hic_data()

        # List of files to get saved
        return (adjlist_filtered_loc, adjlist_loc)


# ------------------------------------------------------------------------------


if __name__ == "__main__":
    start = time.time()
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="Generate adjacency files")
    parser.add_argument("--genome", help="Genome assembly FASTA file") #             default="GCA_000001405.22")
    parser.add_argument("--genome_gem", help="Genome assembly GEM file")
    parser.add_argument("--taxon_id", help="Species (9606)")
    parser.add_argument("--assembly", help="Assembly (GRCh38)")
    parser.add_argument("--file1", help="Location of FASTQ file 1")
    parser.add_argument("--file2", help="Location of FASTQ file 2")
    parser.add_argument("--resolutions", help="CSV string of the resolutions to be computed for the models")
    parser.add_argument("--enzyme_name", help="Enzyme used to digest the DNA")
    parser.add_argument("--windows1", help="FASTQ windowing - start locations", default="1,25,50,75,100")
    parser.add_argument("--windows2", help="FASTQ windowing - paired end locations", default="1,25,50,75,100")
    #parser.add_argument("--file_out")
    #parser.add_argument("--tmp_dir", help="Temporary data dir")
    
    # Get the matching parameters from the command line
    args = parser.parse_args()

    # Assumes that there are 2 fastq files for the paired ends
    #windows1 = ((1,25), (1,50), (1,75),(1,100))
    #windows2 = ((1,25), (1,50), (1,75),(1,100))
    #windows2 = ((101,125), (101,150), (101,175),(101,200))
    
    genome_fa     = args.genome
    genome_gem    = args.genome_gem
    taxon_id      = args.taxon_id
    assembly      = args.assembly
    fastq_01_file = args.file1
    fastq_02_file = args.file2
    tmp_dir       = args.tmp_dir
    enzyme_name   = args.enzyme_name
    resolutions   = args.resolutions
    windows1arg   = args.windows1
    windows2arg   = args.windows1

    if windows1arg is not None:
        w1 = [int(i) for i in windows1arg.split(",")]
        windows1 = ((w1[0], j) for j in w1[1:])
    if windows2arg is not None:
        w2 = [int(i) for i in windows2arg.split(",")]
        windows2 = ((w1[0], j) for j in w2[1:])

    if resolutions is None:
         #resolutions = [1000, 2500, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000, 10000000]
        resolutions = [1000000, 10000000]
    else:
        resolutions = resolutions.split(',')

    da = dmp(test=True)
    
    print da.get_files_by_user("test")
    
    genome_file = da.set_file("test", genome_fa, "fasta", "Assembly", taxon_id, meta_data={'assembly' : assembly})
    genome_idx  = da.set_file("test", genome_gem, "gem", "Assembly Index", taxon_id, meta_data={'assembly' : assembly})
    
    metadata = {
        'user_id'     : 'test',
        'assembly'    : assembly,
        'resolutions' : resolutions,
        'enzyme_name' : enzyme_name,
        'windows1'    : windows1,
        'windows2'    : windows2,
    }

    fastq_01_file_in = da.set_file("test", fastq_01_file, "fastq", "Hi-C", taxon_id, meta_data=metadata)
    fastq_02_file_in = da.set_file("test", fastq_02_file, "fastq", "Hi-C", taxon_id, meta_data=metadata)
    
    print da.get_files_by_user("test")
    
    # 3. Instantiate and launch the App
    #from basic_modules import WorkflowApp
    #app = WorkflowApp()
    
    files = [
        genome_fa,
        genome_gem,
        fastq_01_file_in,
        fastq_02_file_in
    ]

    #results = app.launch(process_hic, files, metadata)

    ph = process_hic()
    results = ph.run(files, metadata)
    
    print results
