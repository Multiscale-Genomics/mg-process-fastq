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

import os

from pycompss.api.parameter import FILE_IN, FILE_OUT
from pycompss.api.task import task

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

from common import common

# ------------------------------------------------------------------------------

class macs2(Tool):
    """
    Tool for peak calling for ChIP-seq data
    """
    
    def __init__(self):
        """
        Init function
        """
        print "MACS2 Peak Caller"
    
    @task(name = IN, bam_file = FILE_IN, bam_file_bg = FILE_IN, peak_bed = FILE_OUT, summit_bed = FILE_OUT, narrowPeak = FILE_OUT, broadPeak = FILE_OUT, gappedPeak = FILE_OUT)
    def macs2_peak_calling(self, name, bam_file, bam_file_bg, peak_bed, summits_bed, narrowPeak, broadPeak, gappedPeak):
        """
        Function to run MACS2 for peak calling on aligned sequence files and
        normalised against a provided background set of alignments.
        
        background might need to be optional.
        
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
        Definitions defined for each of these files have come from the MACS2
        documentation
        
        peak_bed : file
            
        summit_bed : file
            
        narrowPeak : file
            BED6+4 file - ideal for transcription factor binding site
            identification
        broardPeak : file
            BED6+3 file - ideal for histone binding site identification
        gappedPeak : file
            BED12+3 file - Contains a merged set of the broad and narrow peak
            files
        
        All these files are described in the docs at https://github.com/taoliu/MACS
        """
        command_line = 'macs2 callpeak -t ' + bam_file + ' -n ' + name + ' -c ' + bam_file_bg + ' --outdir ' + data_dir + project_id
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()
        
        peak_bed    = name + "_peaks.bed"
        summits_bed = name + "_summits.bed"
        narrowPeak  = name + "_narrowPeak"
        broadPeak   = name + "_broadPeak"
        gappedPeak  = name + "_gappedPeak"
        
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
        
        root_name = bam_file.split("/")
        root_name[-1].replace('.fa', '')
        
        name = root_name[-1]
        
        peak_bed    = '/'.join(root_name) + "_peaks.bed"
        summits_bed = '/'.join(root_name) + "_summits.bed"
        narrowPeak  = '/'.join(root_name) + "_narrowPeak"
        broadPeak   = '/'.join(root_name) + "_broadPeak"
        gappedPeak  = '/'.join(root_name) + "_gappedPeak"
        
        # input and output share most metadata
        output_metadata = dict(
            data_type=metadata[0]["data_type"],
            file_type=metadata[0]["file_type"],
            meta_data=metadata[0]["meta_data"])
        
        # handle error
        if not self.macs2_peak_calling(name, input_files[0], input_files[1], peak_bed, summits_bed, narrowPeak, broadPeak, gappedPeak):
            output_metadata.set_exception(
                Exception(
                    "macs2_peak_calling: Could not process files {}, {}.".format(*input_files)))
            peak_bed    = None
            summits_bed = None
            narrowPeak  = None
            broadPeak   = None
            gappedPeak  = None
        return ([peak_bed, summits_bed, narrowPeak, broadPeak, gappedPeak], [output_metadata])

# ------------------------------------------------------------------------------
