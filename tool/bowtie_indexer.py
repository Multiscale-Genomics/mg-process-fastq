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

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    print ("[Warning] Cannot import \"pycompss\" API packages.")
    print ("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT, IN
    from dummy_pycompss import task
    from dummy_pycompss import compss_wait_on

from basic_modules.tool import Tool

from tool.common import common

# ------------------------------------------------------------------------------

class bowtieIndexerTool(Tool):
    """
    Tool for running indexers over a genome FASTA file
    """

    def __init__(self):
        """
        Init function
        """
        print("Bowtie2 Indexer")
        Tool.__init__(self)

    @task(
        file_loc=FILE_IN, bt_file1=FILE_OUT, bt_file2=FILE_OUT,
        bt_file3=FILE_OUT, bt_file4=FILE_OUT, bt_filer1=FILE_OUT, bt_filer2=FILE_OUT)
    def bowtie2_indexer(
            self, file_loc,
            bt_file1, bt_file2, bt_file3, bt_file4, bt_filer1, bt_filer2): # pylint: disable=unused-argument
        """
        Bowtie2 Indexer

        Parameters
        ----------
        file_loc : str
            Location of the genome assembly FASTA file
        idx_loc : str
            Location of the output index file
        """

        file_name = file_loc.split('/')
        file_name[-1] = file_name[-1].replace('.fasta', '')
        file_name[-1].replace('.fa', '')
        file_name = "/".join(file_name)

        common_handle = common()
        common_handle.bowtie_index_genome(file_loc, file_name + ".tmp")

        with open(bt_file1, "wb") as f_out:
            with open(file_name + ".tmp.1.bt2", "rb") as f_in:
                f_out.write(f_in.read())

        with open(bt_file2, "wb") as f_out:
            with open(file_name + ".tmp.2.bt2", "rb") as f_in:
                f_out.write(f_in.read())

        with open(bt_file3, "wb") as f_out:
            with open(file_name + ".tmp.3.bt2", "rb") as f_in:
                f_out.write(f_in.read())

        with open(bt_file4, "wb") as f_out:
            with open(file_name + ".tmp.4.bt2", "rb") as f_in:
                f_out.write(f_in.read())

        with open(bt_filer1, "wb") as f_out:
            with open(file_name + ".tmp.rev.1.bt2", "rb") as f_in:
                f_out.write(f_in.read())

        with open(bt_filer2, "wb") as f_out:
            with open(file_name + ".tmp.rev.2.bt2", "rb") as f_in:
                f_out.write(f_in.read())

        return True

    def run(self, input_files, output_files, metadata=None):
        """
        Tool for generating assembly aligner index files for use with the
        Bowtie 2 aligner

        Parameters
        ----------
        input_files : list
            List with a single str element with the location of the genome
            assembly FASTA file
        metadata : list

        Returns
        -------
        array : list
            First element is a list of the index files. Second element is a
            list of the matching metadata
        """

        file_name = input_files[0].split('/')
        file_name[-1] = file_name[-1].replace('.fasta', '')
        file_name[-1].replace('.fa', '')

        # input and output share most metadata
        output_metadata = {}

        files_out = [
            "/".join(file_name) + '.1.bt2',
            "/".join(file_name) + '.2.bt2',
            "/".join(file_name) + '.3.bt2',
            "/".join(file_name) + '.4.bt2',
            "/".join(file_name) + '.rev.1.bt2',
            "/".join(file_name) + '.rev.2.bt2',
        ]

        print("BWA INDEXER - files_out:", files_out)

        results = self.bowtie2_indexer(
            input_files[0],
            files_out[0],
            files_out[1],
            files_out[2],
            files_out[3],
            files_out[4],
            files_out[5]
        )

        results = compss_wait_on(results)

        return (files_out, output_metadata)

# ------------------------------------------------------------------------------
