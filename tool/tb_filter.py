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
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT, IN
    from dummy_pycompss import task

from basic_modules.tool import Tool

from pytadbit.mapping.filter import apply_filter
from pytadbit.mapping.filter import filter_reads

# ------------------------------------------------------------------------------

class tbFilterTool(Tool):
    """
    Tool for filtering out experimetnal artifacts from the aligned data
    """

    def __init__(self):
        """
        Init function
        """
        print("TADbit filter aligned reads")
        Tool.__init__(self)

    @task(reads=FILE_IN, filter_reads=FILE_OUT, conservative=IN)
    def tb_filter(self, reads, filter_reads_file, conservative):
        """
        Function to filter out expoerimental artifacts

        Parameters
        ----------
        reads : str
            Location of the reads thats that has a matching location at both
            ends of the paired reads
        filtered_reads_file : str
            Location of the filtered reads
        conservative : bool
            Level of filtering [DEFAULT : True]


        Returns
        -------
        filtered_reads : str
            Location of the filtered reads

        """

        masked = filter_reads(
            reads,
            max_molecule_length=610,
            min_dist_to_re=915,
            over_represented=0.005,
            max_frag_size=100000,
            min_frag_size=100,
            re_proximity=4)

        if conservative is True:
            # Ignore filter 5 (based on docs) as not very helpful
            apply_filter(reads, filter_reads_file, masked, filters=[1, 2, 3, 4, 6, 7, 8, 9, 10])
        else:
            # Less conservative option
            apply_filter(reads, filter_reads_file, masked, filters=[1, 2, 3, 9, 10])

        return True


    def run(self, input_files, metadata):
        """
        The main function to filter the reads to remove experimental artifacts

        Parameters
        ----------
        input_files : list
            reads : str
                Location of the reads thats that has a matching location at both
                ends of the paired reads
        metadata : dict
            conservative : bool
                Level of filtering to apply [DEFAULT : True]


        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects

        """

        reads = input_files[0]
        conservative = True
        if 'conservative_filtering' in metadata:
            conservative = metadata['conservative_filtering']

        root_name = reads.split("/")
        filtered_reads = "/".join(root_name[0:-1]) + '/filtered_map.tsv'

        # input and output share most metadata
        output_metadata = {}

        # handle error
        results = self.tb_filter(reads, filter_reads, conservative)
        return ([filtered_reads], output_metadata)

# ------------------------------------------------------------------------------
