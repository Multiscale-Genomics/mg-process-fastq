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

    @task(
        reads=FILE_IN, filter_reads_file=FILE_OUT, conservative=IN,
        output_de=FILE_OUT, output_d=FILE_OUT, output_e=FILE_OUT,
        output_ed=FILE_OUT, output_or=FILE_OUT, output_rb=FILE_OUT,
        output_sc=FILE_OUT, output_tc=FILE_OUT, output_tl=FILE_OUT,
        output_ts=FILE_OUT)
    def tb_filter(
            self, reads, filter_reads_file, conservative, output_de, output_d,
            output_e, output_ed, output_or, output_rb, output_sc, output_tc,
            output_tl, output_ts):
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


    def run(self, input_files, output_files, metadata=None):
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
        filtered_reads_file = "/".join(root_name[0:-1]) + '/' + metadata['expt_name'] + '_filtered_map.tsv'

        output_de = filtered_reads_file + '_dangling-end.tsv'
        output_d = filtered_reads_file + '_duplicated.tsv'
        output_e = filtered_reads_file + '_error.tsv'
        output_ed = filtered_reads_file + '_extra_dangling-end.tsv'
        output_or = filtered_reads_file + '_over-represented.tsv'
        output_rb = filtered_reads_file + '_random_breaks.tsv'
        output_sc = filtered_reads_file + '_self-circle.tsv'
        output_tc = filtered_reads_file + '_too_close_from_RES.tsv'
        output_tl = filtered_reads_file + '_too_large.tsv'
        output_ts = filtered_reads_file + '_too_short.tsv'

        # input and output share most metadata
        output_metadata = {}

        # handle error
        results = self.tb_filter(
            reads, filtered_reads_file, conservative,
            output_de, output_d, output_e, output_ed, output_or, output_rb,
            output_sc, output_tc, output_tl, output_ts)
        return ([filtered_reads_file], output_metadata)

# ------------------------------------------------------------------------------
