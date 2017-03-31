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
"""process whole genome bisulfate sequencing FastQ files"""
import argparse, time, urllib2, gzip, shutil

try :
    from pycompss.api.parameter import *
    from pycompss.api.task import task
    from pycompss.api.constraint import *
except ImportError :
    print "[Warning] Cannot import \"pycompss\" API packages."
    print "          Using mock decorators."
    
    from dummy_pycompss import *

from tool.common import common
from tool import bs_seeker_aligner
from tool import bs_seeker_filter
from tool import bs_seeker_indexer
from tool import bs_seeker_methylation_caller

import pysam

# ------------------------------------------------------------------------------

class process_wgbs:
    """
    Functions for downloading and processing whole genome bisulfate sequencings
    (WGBS) files. Files are filtered, aligned and analysed for points of
    methylation
    """

    def Splitter(self, in_file1, in_file2, tag = 'tmp'):
        """
        Function to divide the FastQ files into separte sub files of 1000000
        sequences so that the aligner can run in parallel.
        
        Parameters
        ----------
        in_file1 : str
            Location of first paired end FASTQ file
        in_file2 : str
            Location of second paired end FASTQ file
        tag : str
            DEFAULT = tmp
            Tag used to identify the files. Useful if this is getting run
            manually on a single machine multiple times to prevent collisions of
            file names

        
        Returns
        -------
        Returns: Returns a list of lists of the files that have been generated.
                 Each sub list containing the two paired end files for that
                 subset.
        paired_files : list
            List of lists of pair end files. Each sub list containing the two
            paired end files for that subset.
        """
        
        fqr = fastqreader()
        fqr.openPairedFastQ(in_file1, in_file2)
        fqr.createOutputFiles(tag)

        r1 = fqr.next(1)
        r2 = fqr.next(2)

        count_r1 = 0
        count_r2 = 0
        count_r3 = 0
        
        f1 = fqr.fastq1.split("/")
        f1[-1] = f1[-1].replace(".fastq", "." + str(fqr.output_tag) + "_" + str(fqr.output_file_count) + ".fastq")
        f1.insert(-1, tag)
        
        f2 = fqr.fastq2.split("/")
        f2[-1] = f2[-1].replace(".fastq", "." + str(fqr.output_tag) + "_" + str(fqr.output_file_count) + ".fastq")
        f2.insert(-1, tag)
        files_out = [["/".join(f1), "/".join(f2)]]

        while fqr.eof(1) == False and fqr.eof(2) == False:
            r1_id = r1["id"].split(" ")
            r2_id = r2["id"].split(" ")
            
            if r1_id[0] == r2_id[0]:
                fqr.writeOutput(r1, 1)
                fqr.writeOutput(r2, 2)
                
                r1 = fqr.next(1)
                r2 = fqr.next(2)
                
                count_r1 += 1
                count_r2 += 1
                count_r3 += 1
            elif r1_id[0] < r2_id[0]:
                r1 = fqr.next(1)
                count_r1 += 1
            else:
                r2 = fqr.next(2)
                count_r2 += 1
            
            if count_r3 % 1000000 == 0:
                fqr.incrementOutputFiles()
                f1 = fqr.fastq1.split("/")
                f1[-1] = f1[-1].replace(".fastq", "." + str(fqr.output_tag) + "_" + str(fqr.output_file_count) + ".fastq")
                f1.insert(-1, "tmp")
                
                f2 = fqr.fastq2.split("/")
                f2[-1] = f2[-1].replace(".fastq", "." + str(fqr.output_tag) + "_" + str(fqr.output_file_count) + ".fastq")
                f2.insert(-1, "tmp")
                
                files_out.append(["/".join(f1), "/".join(f2)])

        fqr.closePairedFastQ()
        fqr.closeOutputFiles()
        
        return files_out


    def run(self, file_ids, metadata):
        """
        This pipeline processes paired-end FASTQ files to idententify
        methylated regions within the genome.

        Parameters
        ----------
        file_ids : list
            List of strings for the locations of files. These should include:

            genome_fa : str
                Genome assembly in FASTA
            fastq1 : str
                FASTQ file for the first pair end file
            fastq2 : str
                FASTQ file for the second pair end file

        Returns
        -------
        wig_file : str
            Location of the wig file containing the methylation peak calls
        cgmap_file : str
            Location of the CGmap file generated by BS-Seeker2
        atcgmap_file : str
            Location of the ATCGmap file generated by BS-Seeker2
        """
        genome_fa = file_ids[0]
        fastq1 = file_ids[1]
        fastq2 = file_ids[2]

        output_metadata = {}

        # Filter the FASTQ reads to remove duplicates
        frt = bs_seeker_filter.filterReadsTool()
        fastq1f, filter1_meta = frt.run([fastq1], {})
        fastq2f, filter2_meta = frt.run([fastq2], {})

        output_metadata['fastq1'] = filter1_meta
        output_metadata['fastq2'] = filter2_meta

        # Build the matching WGBS genome index
        builder = bs_seeker_indexer.bssIndexerTool()
        genome_idx, gidx_meta = builder.run(
            [genome_fa],
            metadata
        )
        output_metadata['genome_idx'] = gidx_meta

        # Split the FASTQ files into smaller, easier to align packets
        tmp_fastq = self.Splitter(fastq1f[0], fastq2f[0], 'tmp')
        bam_sort_files = []
        bam_merge_files = []
        fastq_for_alignment = []
        for bams in tmp_fastq:
            bam_root = bams[0] + "_bspe.bam"
            tmp = bams
            tmp.append(aligner)
            tmp.append(aligner_path)
            tmp.append(genome_fa)
            tmp.append(bam_root)

            fastq_for_alignment.append(tmp)
            bam_sort_files.append([bam_root, bam_root + ".sorted.bam"])
            bam_merge_files.append(bam_root + ".sorted.bam")
        
        # Run the bs_seeker2-align.py steps on the split up fastq files
        for ffa in fastq_for_alignment:
            bss_aligner = bs_seeker_aligner.bssAlignerTool()
            bam, bam_meta = bss_aligner.run(
                [ffa[0], ffa[1], ffa[2], ffa[3][ffa[2]], ffa[4], ffa[5]],
                {}
            )
            if 'alignment' in output_metadata:
                output_metadata['alignment'] = bam_meta
            else:
                output_metadata['alignment'].append(bam_meta)

        # Sort and merge the aligned bam files
        # Pre-sort the original input bam files
        for bfs in bam_sort_files:
            pysam.sort("-o", bfs[1], bfs[0])
        
        f_bam = fastq1.split("/")
        f_bam[-1] = f_bam[-1].replace(".fastq", ".sorted.bam")
        out_bam_file = "/".join(f_bam)

        pysam.merge(out_bam_file, *bam_merge_files)
        
        pysam.sort("-o", out_bam_file, "-T", out_bam_file + "_sort", out_bam_file)
        
        pysam.index(out_bam_file)

        # Methylation peak caller
        wig_file     = out_bam_file.replace('.bam', '.wig')
        cgmap_file   = out_bam_file.replace('.bam', '.cgmap')
        atcgmap_file = out_bam_file.replace('.bam', '.atcgmap')
        mc = bs_seeker_methylation_caller.bssMethylationCallerTool()
        peak_calls, peak_meta = mc.run(
            [
                metadata['aligner_path'], out_bam_file, wig_file,
                cgmap_file, atcgmap_file, genome_idx
            ],
            {}
        )
        output_metadata['peak_calling'] = peak_meta




        return (wig_file, cgmap_file, atcgmap_file)


# ------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    import os
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="Parse WGBS data")
    parser.add_argument("--fastq1", help="Location of first paired end FASTQ")
    parser.add_argument("--fastq2", help="Location of second paired end FASTQ")
    parser.add_argument("--taxon_id", help="Taxon_ID (10090)")
    parser.add_argument("--assembly", help="Assembly (GRCm38)")
    parser.add_argument("--genome", help="Genome assembly FASTA file")
    parser.add_argument("--aligner", help="Aligner to use (eg bowtie2)")
    parser.add_argument("--aligner_path", help="Directory for the aligner program")

    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    fastq1 = args.fastq1
    fastq2 = args.fastq2
    genome   = args.genome
    taxon_id   = args.taxon_id
    assembly    = args.assembly
    aligner  = args.aligner
    aligner_path = args.aligner_path
    
    metadata = {
        'user_id' : 'test',
        'aligner' : aligner,
        'aligner_path' : aligner_path
    }


    da = dmp(test=True)
    
    print da.get_files_by_user("test")
    
    genome_file = da.set_file("test", genome_fa, "fasta", "Assembly", taxon_id, meta_data={'assembly' : assembly})
    fastq_file_1 = da.set_file("test", fastq1, "fastq", "wgbs", taxon_id, meta_data={'assembly' : assembly})
    fastq_file_2 = da.set_file("test", fastq2, "fastq", "wgbs", taxon_id, meta_data={'assembly' : assembly})
    
    print da.get_files_by_user("test")

    files = [
        genome_fa,
        fastq1,
        fastq2
    ]
    
    # 3. Instantiate and launch the App
    #from basic_modules import WorkflowApp
    #app = WorkflowApp()
    #results = app.launch(process_wgbs, [genome_fa, fastq_file_1, fastq_file_2], metadata)
    
    pw = process_wgbs()
    peak_file = pw.run(files, metadata)
    
    print peak_files