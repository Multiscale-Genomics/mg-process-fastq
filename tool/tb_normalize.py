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

class tbNormalizeTool(Tool):
    """
    Tool for normalizing an adjacency matrix
    """

    def __init__(self):
        """
        Init function
        """
        print("TADbit - Normalize")
        Tool.__init__(self)

    @task(bamin=FILE_IN, resolution=IN, min_perc=IN, max_perc=IN, workdir=IN)
    # @constraint(ProcessorCoreCount=16)
    def tb_normalize(self, bamin, resolution, min_perc, max_perc, workdir):
        """
        Function to the predict TAD sites for a given resolution from the Hi-C
        matrix

        Parameters
        ----------
        expt_name : str
                Location of the adjacency list
        matrix_file : str
            Location of the HDF5 output matrix file
        resolution : int
            Resolution to read the Hi-C adjacency list at
        tad_file : str
            Location of the output TAD file

        Returns
        -------
        tad_file : str
            Location of the output TAD file

        """
        #chr_hic_data = read_matrix(matrix_file, resolution=int(resolution))

        print("TB NORMALIZATION:",bamin, resolution, min_perc, max_perc, workdir)

        _cmd = [
                'tadbit', 'normalize', 
            '--bam', bamin,
            '--workdir', workdir,
            '--resolution', resolution,
            '--cpus', '32'
            ]
    
        if min_perc:
            _cmd.append('--min_perc')
            _cmd.append(min_perc)
        if max_perc:
            _cmd.append('--max_perc')
            _cmd.append(max_perc)
        
        try:
            out, err = Popen(_cmd, stdout=PIPE, stderr=PIPE).communicate()
        except CalledProcessError as e:
            print(out)
            print(err)
            raise Exception(e.output)

        return True

    def run(self, input_files, output_files, metadata=None):
        """
        The main function to the predict TAD sites for a given resolution from
        the Hi-C matrix

        Parameters
        ----------
        input_files : list
            adj_list : str
                Location of the adjacency list
        metadata : dict
            resolutions : list
                Levels of resolution for the adjacency list to be daved at
            assembly : str
                Assembly of the aligned sequences



        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects

        """

        bamin = input_files[0]

        resolution = '1000000'
        if 'resolution' in metadata:
            resolution = metadata['resolution']

        min_perc = max_perc = None 
        
        if 'min_perc' in metadata:
            min_perc = metadata['min_perc']
        if 'max_perc' in metadata:
            max_perc = metadata['max_perc']
            
        root_name = bamin.split("/")
        if 'workdir' in metadata:
            root_name = metadata['workdir']
        
        # input and output share most metadata
        output_metadata = {}
        
        hic_biases = 'hic_biases' 
        self.tb_normalize(bamin, resolution, min_perc, max_perc, root_name)

        return ([hic_biases], output_metadata)

# ------------------------------------------------------------------------------
