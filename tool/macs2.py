"""
.. Copyright 2016 EMBL-European Bioinformatics Institute
 
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

from .. import common

# ------------------------------------------------------------------------------

class macs2(Tool):
    """
    Tool for peak calling for ChIP-seq data
    """
    
    @task(name = IN, bam_file = FILE_IN, bam_file_bg = FILE_IN, peak_bed = FILE_OUT, summit_bed = FILE_OUT, narrowPeak = FILE_OUT, broadPeak = FILE_OUT, gappedPeak = FILE_OUT)
    def macs2_peak_calling(self, name, bam_file, bam_file_bg):
        """
        Function to run MACS2 for peak calling
        
        background might need to be optional.
        
        Parameters
        ----------
        name : str
        bam_file : str
        bam_file_bgd : str
        
        Returns
        -------
        peak_bed : file
        summit_bed : file
        narrowPeak : file
        broardPeak : file
        gappedPeak : file
        
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
        Standard function to call a task
        """
        
        # handle error
        if not self.macs2_peak_calling("TODO_Name", input_files[0], input_files[1]):
            output_metadata.set_exception(
                Exception(
                    "macs2_peak_calling: Could not process files {}, {}.".format(*input_files)))
output_file = None
        return ([output_file], [output_metadata])

# ------------------------------------------------------------------------------
