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

import argparse, urllib2, gzip, shutil, shlex, subprocess, os.path, json

from basic_modules.tool import Tool
from basic_modules.workflow import Workflow
from basic_modules.metadata import Metadata


from functools import wraps

from dmp import dmp

from tool import bwa_aligner
from tool import biobambam_filter
from tool import macs2

import os

try:
    from pycompss.api.parameter import FILE_IN, FILE_OUT
    from pycompss.api.task import task
    from pycompss.api.constraint import constraint
    from pycompss.api.api import compss_wait_on
except ImportError :
    print "[Warning] Cannot import \"pycompss\" API packages."
    print "          Using mock decorators."
    
    from dummy_pycompss import *

try :
    import pysam
except ImportError :
    print "[Error] Cannot import \"pysam\" package. Have you installed it?"
    exit(-1)

# ------------------------------------------------------------------------------

class process_chipseq(Workflow):
    """
    Functions for processing Chip-Seq FastQ files. Files are the aligned,
    filtered and analysed for peak calling
    """
    
    configuration = {}
    
    def __init__(self, configuration={}):
        """
        Initialise the tool with its configuration.


        Parameters
        ----------
        configuration : dict
            a dictionary containing parameters that define how the operation
            should be carried out, which are specific to each Tool.
        """
        self.configuration.update(configuration)

    
    def run(self, file_ids, metadata, output_files):
        """
        Main run function for processing ChIP-seq FastQ data. Pipeline aligns
        the FASTQ files to the genome using BWA. MACS 2 is then used for peak
        calling to identify transcription factor binding sites within the
        genome.
                
        Parameters
        ----------
        files_ids : list
            List of file locations
        metadata : list
            Required meta data
        output_files : list
            List of output file locations
        
        Returns
        -------
        outputfiles : list
            List of locations for the output bam, bed and tsv files
        """
        
        # TODO - Handle multiple file and background files
        genome_fa = file_ids[0]
        bwa_amb = file_ids[1]
        bwa_ann = file_ids[2]
        bwa_bwt = file_ids[3]
        bwa_pac = file_ids[4]
        bwa_sa  = file_ids[5]
        file_loc = file_ids[6]
        file_bgd_loc = file_ids[7]
        

        
        out_bam = file_loc.replace(".fastq", '.bam')
        
        bwa = bwa_aligner.bwaAlignerTool(self.configuration)
        out_bam = file_loc.replace(".fastq", '.bam')
        bwa_results = bwa.run(
            [genome_fa, file_loc, bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_sa],
            {},
            [out_bam]
        )
        
        #bwa_results = compss_wait_on(bwa_results)
        
        if file_bgd_loc != None:
            out_bgd_bam = file_bgd_loc.replace(".fastq", '.bam')
            bwa_results_bgd = bwa.run(
                [genome_fa, file_bgd_loc, bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_sa],
                {},
                [out_bgd_bam]
            )
            #bwa_results_bgd = compss_wait_on(bwa_results_bgd)
            
        
        # TODO - Multiple files need merging into a single bam file
       
        b3f_file_bgd_out = ''

        # Filter the bams
        b3f = biobambam_filter.biobambam(self.configuration)
        b3f_file_out = file_loc.replace('.fastq', '.filtered.bam')
        b3f_results = b3f.run([out_bam], {}, [b3f_file_out])
        
        #b3f_results = compss_wait_on(b3f_results)
        
        if file_bgd_loc != None:
            b3f_bgd_file_out = file_bgd_loc.replace(".fastq", '.filtered.bam')
            b3f_results_bgd = b3f.run([out_bgd_bam, b3f_bgd_file_out], {})
            #b3f_results_bgd = compss_wait_on(b3f_results_bgd)
        else:
            b3f_bgd_file_out = None
        
        # MACS2 to call peaks
        m = macs2.macs2(self.configuration)
        
        mac_root_name = b3f_file_out.split("/")
        mac_root_name[-1] = mac_root_name[-1].replace('.bam', '')
        
        name = root_name[-1]
        
        summits_bed = '/'.join(mac_root_name) + "_summits.bed"
        narrowPeak  = '/'.join(mac_root_name) + "_narrowPeak"
        broadPeak   = '/'.join(mac_root_name) + "_broadPeak"
        gappedPeak  = '/'.join(mac_root_name) + "_gappedPeak"
        
        output_files = [
            summits_bed,
            narrowPeak,
            broadPeak,
            gappedPeak
        ]
        
        if file_bgd_loc != None:
            m_results = m.run([b3f_file_out,  b3f_bgd_file_out], {}, output_files)
        else:
            m_results = m.run([b3f_file_out], {}, output_files)
        
        #m_results = compss_wait_on(m_results)
        
        return ([b3f_file_out, b3f_file_bgd_out] + m_results[0], [b3f_results[1],m_results[1]])


# -----------------------------------------------------------------------------

def main(inputFiles, inputMetadata, outputFiles):
    """
    Main function
    -------------

    This function launches the app.
    """

    # import pprint  # Pretty print - module for dictionary fancy printing

    # 1. Instantiate and launch the App
    print "1. Instantiate and launch the App"
    from apps.workflowapp import WorkflowApp
    app = WorkflowApp()
    result = app.launch(process_chipseq, inputFiles, inputMetadata,
                        outputFiles, {})

    # 2. The App has finished
    print "2. Execution finished"
    print result
    return result

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    import os
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="ChIP-seq peak calling")
    parser.add_argument("--taxon_id", help="Taxon_ID (9606)")
    parser.add_argument("--genome", help="Genome FASTA file")
    parser.add_argument("--assembly", help="Genome assembly ID (GCA_000001405.25)")
    parser.add_argument("--file", help="Location of FASTQ input file")
    parser.add_argument("--bgd_file", help="Location of FASTQ background file", default=None)
    
    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    taxon_id    = args.taxon_id
    genome_fa   = args.genome
    assembly    = args.assembly
    file_loc    = args.file
    file_bg_loc = args.bgd_file
    
    #
    # MuG Tool Steps
    # --------------
    # 
    # 1. Create data files
    
    # Get the assembly
    
    #2. Register the data with the DMP
    da = dmp(test=True)
    
    print da.get_files_by_user("test")
    
    genome_file = da.set_file("test", genome_fa, "fasta", "Assembly",
        taxon_id, None, [], meta_data={"assembly" : assembly})
    genome_file_idx1 = da.set_file("test", genome_fa + ".amb", "amb", "Assembly",
        taxon_id, None, [genome_file], meta_data={'assembly' : assembly})
    genome_file_idx2 = da.set_file("test", genome_fa + ".ann", "ann", "Assembly",
        taxon_id, None, [genome_file], meta_data={'assembly' : assembly})
    genome_file_idx3 = da.set_file("test", genome_fa + ".bwt", "bwt", "Assembly",
        taxon_id, None, [genome_file], meta_data={'assembly' : assembly})
    genome_file_idx4 = da.set_file("test", genome_fa + ".pac", "pac", "Assembly",
        taxon_id, None, [genome_file], meta_data={'assembly' : assembly})
    genome_file_idx5 = da.set_file("test", genome_fa + ".sa", "sa", "Assembly",
        taxon_id, None, [genome_file], meta_data={'assembly' : assembly})

    # Read metadata file and build a dictionary with the metadata:
    from basic_modules.metadata import Metadata
    # Maybe it is necessary to prepare a metadata parser from json file
    # when building the Metadata objects.
    metadata = [
        Metadata("fasta", "Assembly"),
        Metadata("index", "Assembly"),
        Metadata("index", "Assembly"),
        Metadata("index", "Assembly"),
        Metadata("index", "Assembly"),
        Metadata("index", "Assembly"),
        Metadata("fasta", "ChIP-seq")
    ]

    file_in = da.set_file("test", file_loc, "fasta", "ChIP-seq", taxon_id,
        None, [], meta_data={'assembly' : assembly})
    if file_bg_loc:
        file_bg_in = da.set_file("test", file_bg_loc, "fasta", "ChIP-seq",
            taxon_id, None, [], meta_data={'assembly' : assembly})
        metadata.append(Metadata("fasta", "ChIP-seq"))
    
    print da.get_files_by_user("test")

    files = [
        genome_fa,
        genome_fa + ".amb",
        genome_fa + ".ann",
        genome_fa + ".bwt",
        genome_fa + ".pac",
        genome_fa + ".sa",
        file_loc,
        file_bg_loc
    ]

    out_bam = file_loc.replace(".fastq", '.filtered.bam')
    out_peaks_narrow = file_loc.replace(".fastq", '_peaks.narrowPeak')
    out_peaks_xls = file_loc.replace(".fastq", '_peaks.xls')
    out_peaks_broad = file_loc.replace(".fastq", '_summits.bed')
    out_peaks_gapped = file_loc.replace(".fastq", '_summits.bed')
    out_summits = file_loc.replace(".fastq", '_summits.bed')


    files_out = [
        out_bam,
        out_peaks_narrow,
        out_peaks_xls,
        out_peaks_broad,
        out_peaks_gapped,
        out_summits
    ]

    # 3. Instantiate and launch the App
    #from basic_modules import WorkflowApp
    #app = WorkflowApp()
    #results = app.launch(process_chipseq, [genome_file, file_in, file_bg_in], {'user_id' : 'test'})

    #pc = process_chipseq()
    #results_files, results_meta = pc.run(files, {"user_id" : "test"})

    #results_files, results_meta = main(files, metadata, files_out)
    results = main(files, metadata, files_out)

    #print(results_files)
    #print(results_meta)
    print(results)
    
    print da.get_files_by_user("test")
    
