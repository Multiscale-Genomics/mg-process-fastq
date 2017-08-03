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
import os.path

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT, IN
    from dummy_pycompss import task
    from dummy_pycompss import compss_wait_on

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
        
        reads_tmp = reads.replace(".tsv", '')
        
        with open(reads_tmp + "_tmp.tsv", "wb") as f_out:
            with open(reads, "rb") as f_in:
                f_out.write(f_in.read())
        
        masked = filter_reads(
            reads_tmp + "_tmp.tsv",
            max_molecule_length=610,
            min_dist_to_re=915,
            over_represented=0.005,
            max_frag_size=100000,
            min_frag_size=100,
            re_proximity=4)
        
        filter_reads_file_tmp = filter_reads_file.replace(".tsv", '')
        
        if conservative is True:
            # Ignore filter 5 (based on docs) as not very helpful
            apply_filter(reads_tmp + "_tmp.tsv", filter_reads_file_tmp + "_tmp.tsv", masked, filters=[1, 2, 3, 4, 6, 7, 8, 9, 10])
        else:
            # Less conservative option
            apply_filter(reads_tmp + "_tmp.tsv", filter_reads_file_tmp + "_tmp.tsv", masked, filters=[1, 2, 3, 9, 10])
        
        with open(filter_reads_file, "wb") as f_out:
            with open(filter_reads_file_tmp + "_tmp.tsv", "rb") as f_in:
                f_out.write(f_in.read())
        
        filters_suffixes = ['dangling-end', 'duplicated', 'error', 'extra_dangling-end', 'over-represented', 'random_breaks', 'self-circle', 'too_close_from_RES', 'too_large', 'too_short']
        for i in filters_suffixes:
            report_file_loc = reads_tmp + '_tmp.tsv_' + i + '.tsv'
            print(report_file_loc)
            if os.path.isfile(report_file_loc) is True:
                print("- Present", os.path.getsize(report_file_loc))
                with open(report_file_loc, "rb") as f_in:
                    if i == 'dangling-end':
                        print("- Saving to:", output_de)
                        with open(output_de, "wb") as f_out:
                            f_out.write(f_in.read())
                    elif i == 'duplicated':
                        print("- Saving to:", output_d)
                        with open(output_d, "wb") as f_out:
                            f_out.write(f_in.read())
                    elif i == 'error':
                        print("- Saving to:", output_e)
                        with open(output_e, "wb") as f_out:
                            f_out.write(f_in.read())
                    elif i == 'extra_dangling-end':
                        print("- Saving to:", output_ed)
                        with open(output_ed, "wb") as f_out:
                            f_out.write(f_in.read())
                    elif i == 'over-represented':
                        print("- Saving to:", output_or)
                        with open(output_or, "wb") as f_out:
                            f_out.write(f_in.read())
                    elif i == 'random_breaks':
                        print("- Saving to:", output_rb)
                        with open(output_rb, "wb") as f_out:
                            f_out.write(f_in.read())
                    elif i == 'self-circle':
                        print("- Saving to:", output_sc)
                        with open(output_sc, "wb") as f_out:
                            f_out.write(f_in.read())
                    elif i == 'too_close_from_RES':
                        print("- Saving to:", output_tc)
                        with open(output_tc, "wb") as f_out:
                            f_out.write(f_in.read())
                    elif i == 'too_large':
                        print("- Saving to:", output_tl)
                        with open(output_tl, "wb") as f_out:
                            f_out.write(f_in.read())
                    elif i == 'too_short':
                        print("- Saving to:", output_ts)
                        with open(output_ts, "wb") as f_out:
                            f_out.write(f_in.read())

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
        results = compss_wait_on(results)

        return ([filtered_reads_file], output_metadata)

# ------------------------------------------------------------------------------
