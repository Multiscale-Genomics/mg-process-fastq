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

import os, shutil, shlex, subprocess

try :
    from pycompss.api.parameter import FILE_IN, FILE_OUT, FILE_INOUT, IN
    from pycompss.api.task import task
except ImportError :
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")
    
    from dummy_pycompss import *

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

from . import common

# ------------------------------------------------------------------------------

class macs2(Tool):
    """
    Tool for peak calling for ChIP-seq data
    """
    
    def __init__(self):
        """
        Init function
        """
        print ("MACS2 Peak Caller")
    
    @task(
            name = IN,
            bam_file = FILE_IN,
            bam_file_bgd = FILE_IN,
            narrowPeak = FILE_OUT,
            summit_bed = FILE_OUT,
            broadPeak = FILE_OUT,
            gappedPeak = FILE_OUT)
    def macs2_peak_calling(self,
        name, bam_file, bam_file_bgd = None,
        narrowPeak = None, summits_bed = None, broadPeak = None, gappedPeak = None):
        """
        Function to run MACS2 for peak calling on aligned sequence files and
        normalised against a provided background set of alignments.
        
        Parameters
        ----------
        name : str
            Name to be used to identify the files
        bam_file : str
            Location of the aligned FASTQ files as a bam file
        bam_file_bgd : str
            Location of the aligned FASTQ files as a bam file representing
            background values for the cell
        
        
        Returns
        -------
        narrowPeak : file
            BED6+4 file - ideal for transcription factor binding site
            identification
        summitPeak : file
            BED4+1 file - Contains the peak summit locations for everypeak
        broadPeak : file
            BED6+3 file - ideal for histone binding site identification
        gappedPeak : file
            BED12+3 file - Contains a merged set of the broad and narrow peak
            files
        
        Definitions defined for each of these files have come from the MACS2
        documentation described in the docs at https://github.com/taoliu/MACS
        """
        od = bam_file.split("/")
        output_dir = "/".join(od[0:-1])

        bgd_command = ''
        if bam_file_bgd is not None:
            bgd_command = ' -c ' + bam_file_bgd
        
        command_line = 'macs2 callpeak -t ' + bam_file + ' -n ' + name + bgd_command + ' --outdir ' + output_dir
        
        if name == 'fastQForSelRegion':
            # This is for when running the test data
            command_line = command_line + ' --nomodel'
        
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()
        
        return True
    
    @task(
            name = IN,
            bam_file = FILE_IN,
            narrowPeak = FILE_OUT,
            summit_bed = FILE_OUT,
            broadPeak = FILE_OUT,
            gappedPeak = FILE_OUT)
    def macs2_peak_calling_nobgd(self,
        name, bam_file,
        narrowPeak = None, summits_bed = None, broadPeak = None, gappedPeak = None):
        """
        Function to run MACS2 for peak calling on aligned sequence files without
        a background dataset for normalisation.
        
        Parameters
        ----------
        name : str
            Name to be used to identify the files
        bam_file : str
            Location of the aligned FASTQ files as a bam file
        
        Returns
        -------
        narrowPeak : file
            BED6+4 file - ideal for transcription factor binding site
            identification
        summitPeak : file
            BED4+1 file - Contains the peak summit locations for everypeak
        broadPeak : file
            BED6+3 file - ideal for histone binding site identification
        gappedPeak : file
            BED12+3 file - Contains a merged set of the broad and narrow peak
            files
        
        Definitions defined for each of these files have come from the MACS2
        documentation described in the docs at https://github.com/taoliu/MACS
        """
        od = bam_file.split("/")
        output_dir = "/".join(od[0:-1])

        command_line = 'macs2 callpeak -t ' + bam_file + ' -n ' + name + ' --outdir ' + output_dir
        
        if name == 'fastQForSelRegion':
            # This is for when running the test data
            command_line = command_line + ' --nomodel'
        
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()
        
        return True
    
    
    def run(self, input_files, metadata):
        """
        The main function to run MACS 2 for peak calling over a given BAM file
        and matching background BAM file.
        
        Parameters
        ----------
        input_files : list
            List of input bam file locations where 0 is the bam data file and 1
            is the matching background bam file
        metadata : dict
        
        
        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects
        
        """
        
        bam_file = input_files[0]
        bam_file_bgd = None
        if len(input_files) == 2 and input_files[1] is not None:
            bam_file_bgd = input_files[1]
        
        root_name = bam_file.split("/")
        root_name[-1] = root_name[-1].replace('.bam', '')
        
        name = root_name[-1]
        
        summits_bed = '/'.join(root_name) + "_summits.bed"
        narrowPeak  = '/'.join(root_name) + "_narrowPeak"
        broadPeak   = '/'.join(root_name) + "_broadPeak"
        gappedPeak  = '/'.join(root_name) + "_gappedPeak"
        
        # input and output share most metadata
        output_metadata = {"bed_types" : ["bed4+1", "bed6+4", "bed6+3", "bed12+3"]}
        
        peak_bed = None
        # handle error
        if bam_file_bgd is None:
            results = self.macs2_peak_calling_nobgd(name, bam_file,
                narrowPeak, summits_bed, broadPeak, gappedPeak)
            #if not self.macs2_peak_calling_nobgd(name, bam_file,
            #    narrowPeak, summits_bed, broadPeak, gappedPeak):
            #    output_metadata['error'] = Exception(
            #            "macs2_peak_calling_nobgd: Could not process files {}, {}."#.format(*input_files)
            #    )
            #    peak_bed    = None
            #    summits_bed = None
            #    narrowPeak  = None
            #    broadPeak   = None
            #    gappedPeak  = None
        else:
            results = self.macs2_peak_calling(name, bam_file, bam_file_bgd,
                narrowPeak, summits_bed, broadPeak, gappedPeak)
            #if not self.macs2_peak_calling(name, bam_file, bam_file_bgd,
            #    narrowPeak, summits_bed, broadPeak, gappedPeak):
            #    output_metadata['error'] = Exception(
            #            "macs2_peak_calling: Could not process files {}, {}.".format(*input_files)
            #    )
            #    peak_bed    = None
            #    summits_bed = None
            #    narrowPeak  = None
            #    broadPeak   = None
            #    gappedPeak  = None
        return ([peak_bed, summits_bed, narrowPeak, broadPeak, gappedPeak], output_metadata)

# ------------------------------------------------------------------------------
