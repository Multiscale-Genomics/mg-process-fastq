"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.

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
from __future__ import print_function

import sys
import glob, os
from subprocess import CalledProcessError, PIPE, Popen

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, FILE_INOUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
    # from pycompss.api.constraint import constraint
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT, FILE_INOUT, IN
    from dummy_pycompss import task
    from dummy_pycompss import compss_wait_on
    #from dummy_pycompss import constraint

from basic_modules.tool import Tool

# ------------------------------------------------------------------------------

class tbSegmentTool(Tool):
    """
    Tool for finding tads and compartments in an adjacency matrix
    """

    def __init__(self):
        """
        Init function
        """
        print("TADbit - Normalize")
        Tool.__init__(self)

    @task(bamin=FILE_IN, biases=FILE_IN, resolution=IN, workdir=IN)
    # @constraint(ProcessorCoreCount=16)
    def tb_segment(self, bamin, biases, resolution, callers, chromosomes, workdir, fasta=None ,ncpus="1"):
        """
        Function to find tads and compartments in the Hi-C
        matrix

        Parameters
        ----------
        bamin : str
            Location of the tadbit bam paired reads
        biases : str
            Location of the pickle hic biases
        resolution : int
            Resolution of the Hi-C
        callers: str
            1 for ta calling, 2 for compartment calling
        workdir : str
            Location of working directory
        ncpus : int
            Number of cpus to use
                
        Returns
        -------
        compartments : str
            Location of tsv file with compartment definition
        tads : str
            Location of tsv file with tad definition
        filtered_bins : str
            Location of filtered_bins png

        """
        print("TB SEGMENT:",bamin, resolution, workdir)

        _cmd = [
                'tadbit', 'segment', 
            '--nosql', '--mreads', bamin,
            '--workdir', workdir,
            '--resolution', resolution,
            '--cpu', str(ncpus)
            ]
    
        if '2' not in callers:
            _cmd.append('--only_tads')
        if '1' not in callers:
            _cmd.append('--only_compartments')
        
        if chromosomes:
            _cmd.append('--chromosomes')
            _cmd.append(chromosomes)
        if fasta:
            _cmd.append('--fasta')
            _cmd.append(fasta)
        
        if biases:
            _cmd.append('--biases')
            _cmd.append(biases)
            
        output_metadata = {}
        output_files = []
        
        out, err = Popen(_cmd, stdout=PIPE, stderr=PIPE).communicate()
        print(out)
        print(err)

        if '1' in callers:
            tad_dir = os.path.join(workdir, '06_segmentation',
                                 'tads_%s' % (nice(int(resolution))))
            output_files.append(tad_dir) 
        if '2' in callers:
            cmprt_dir = os.path.join(workdir, '06_segmentation',
                                  'compartments_%s' % (nice(int(resolution))))
            output_files.append(cmprt_dir)
            
        return (output_files, output_metadata)

    def run(self, input_files, output_files, metadata=None):
        """
        The main function to the predict TAD sites and compartments for a given resolution from
        the Hi-C matrix

        Parameters
        ----------
        input_files : list
            bamin : str
                Location of the tadbit bam paired reads
            biases : str
                Location of the pickle hic biases
        metadata : dict
            resolution : int
                Resolution of the Hi-C
            workdir : str
                Location of working directory
            ncpus : int
                Number of cpus to use



        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects

        """

        bamin = input_files[0]
        
        if not os.path.isfile(bamin.replace('.bam','.bam.bai')):
            print('Creating bam index')
            _cmd = ['samtools', 'index', bamin]
            out, err = Popen(_cmd, stdout=PIPE, stderr=PIPE).communicate()
            print(out)
            print(err)
        
        resolution = '1000000'
        if 'resolution' in metadata:
            resolution = metadata['resolution']
            
        ncpus=1
        if 'ncpus' in metadata:
            ncpus = metadata['ncpus']
        biases = chromosomes = fasta = None
        if len(input_files) > 1:
            biases = input_files[1]
        if "chromosomes" in metadata:
            chromosomes = metadata['chromosomes']
        if "fasta" in metadata:
            fasta = metadata['fasta']
        callers = "1"
        if "callers" in metadata:
            callers = metadata['callers']
            
        root_name = bamin.split("/")
        if 'workdir' in metadata:
            root_name = metadata['workdir']
        
        # input and output share most metadata
        
        
        output_files, output_metadata = self.tb_segment(bamin, biases, resolution, callers, chromosomes, root_name, fasta, ncpus)

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------

def nice(reso):
    if reso >= 1000000:
        return '%dMb' % (reso / 1000000)
    return '%dkb' % (reso / 1000)
