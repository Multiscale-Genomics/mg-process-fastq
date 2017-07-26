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

from __future__ import print_function
import sys

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    from pycompss.api.constraint import constraint
    from pycompss.api.api import compss_wait_on
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT, IN
    from dummy_pycompss import task
    from dummy_pycompss import constraint
    from dummy_pycompss import compss_wait_on

from basic_modules.tool import Tool

import pytadbit
from pytadbit.mapping.mapper import full_mapping

# ------------------------------------------------------------------------------

class tbFullMappingTool(Tool):
    """
    Tool for mapping fastq paired end files to the GEM index files
    """

    def __init__(self):
        """
        Init function
        """
        print("TADbit full_mapping")
        Tool.__init__(self)

    @task(
        gem_file=FILE_IN, fastq_file=FILE_IN, windows=IN, window1=FILE_OUT,
        window2=FILE_OUT, window3=FILE_OUT, window4=FILE_OUT)
    @constraint(ProcessorCoreCount=32)
    def tb_full_mapping_iter(
            self, gem_file, fastq_file, windows,
            window1, window2, window3, window4):
        """
        Function to map the FASTQ files to the GEM file over different window
        sizes ready for alignment

        Parameters
        ----------
        gem_file : str
            Location of the genome GEM index file
        fastq_file_bgd : str
            Location of the FASTQ file
        windows : list
            List of lists with the window sizes to be computed
        window1 : str
            Location of the first window index file
        window2 : str
            Location of the second window index file
        window3 : str
            Location of the third window index file
        window4 : str
            Location of the fourth window index file


        Returns
        -------
        window1 : str
            Location of the first window index file
        window2 : str
            Location of the second window index file
        window3 : str
            Location of the third window index file
        window4 : str
            Location of the fourth window index file

        """
        print("tb_full_mapping_iter")
        od_loc = fastq_file.split("/")
        output_dir = "/".join(od_loc[0:-1])

        map_files = full_mapping(
            gem_file, fastq_file, output_dir,
            windows=windows, frag_map=False, nthreads=32, clean=True,
            temp_dir='/tmp/'
        )

        return True


    @task(
        gem_file=FILE_IN, fastq_file=FILE_IN, enzyme_name=IN, windows=IN,
        full_file=FILE_OUT, frag_file=FILE_OUT)
    @constraint(ProcessorCoreCount=16)
    def tb_full_mapping_frag(
            self, gem_file, fastq_file, enzyme_name, windows,
            full_file, frag_file):
        """
        Function to map the FASTQ files to the GEM file based on fragments
        derived from the restriction enzyme that was used.

        Parameters
        ----------
        gem_file : str
            Location of the genome GEM index file
        fastq_file_bgd : str
            Location of the FASTQ file
        enzyme_name : str
            Restriction enzyme name (MboI)
        windows : list
            List of lists with the window sizes to be computed
        window_file : str
            Location of the first window index file


        Returns
        -------
        window_file : str
            Location of the window index file

        """
        print("tb_full_mapping_frag")
        od_loc = fastq_file.split("/")
        output_dir = "/".join(od_loc[0:-1])

        map_files = full_mapping(
            gem_file, fastq_file, output_dir,
            r_enz=enzyme_name, windows=windows, frag_map=True, nthreads=32,
            clean=True, temp_dir='/tmp/'
        )

        return True


    def run(self, input_files, output_files, metadata=None):
        """
        The main function to map the FASTQ files to the GEM file over different
        window sizes ready for alignment

        Parameters
        ----------
        input_files : list
            gem_file : str
                Location of the genome GEM index file
            fastq_file_bgd : str
                Location of the FASTQ file
        metadata : dict
            windows : list
                List of lists with the window sizes to be computed
            enzyme_name : str
                Restriction enzyme used [OPTIONAL]


        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects

        """

        gem_file = input_files[0]
        fastq_file = input_files[1]
        windows = metadata['windows']

        root_name = fastq_file.split("/")
        root_name[-1] = root_name[-1].replace('.fasta', '')
        root_name[-1] = root_name[-1].replace('.fa', '')

        #name = root_name[-1]

        # input and output share most metadata
        output_metadata = {}

        root_path = '/'.join(root_name)

        if 'enzyme_name' in metadata:
            full_file = root_path + "_full_" + str(windows[0][0]) + "-" + str(windows[0][1]) + ".map"
            frag_file = root_path + "_frag_" + str(windows[0][0]) + "-" + str(windows[0][1]) + ".map"

            results = self.tb_full_mapping_frag(
                gem_file, fastq_file, metadata['enzyme_name'], windows,
                full_file, frag_file
            )
            results = compss_wait_on(results)

            output_metadata['func'] = 'frag'
            return ([full_file, frag_file], output_metadata)

        window1 = root_path + "_full_" + str(windows[0][0]) + "-" + str(windows[0][1]) + ".map"
        window2 = root_path + "_full_" + str(windows[1][0]) + "-" + str(windows[1][1]) + ".map"
        window3 = root_path + "_full_" + str(windows[2][0]) + "-" + str(windows[2][1]) + ".map"
        window4 = root_path + "_full_" + str(windows[3][0]) + "-" + str(windows[3][1]) + ".map"

        results = self.tb_full_mapping_iter(
            gem_file, fastq_file, windows,
            window1, window2, window3, window4
        )
        results = compss_wait_on(results)

        output_metadata['func'] = 'iter'
        return ([window1, window2, window3, window4], output_metadata)

# ------------------------------------------------------------------------------
