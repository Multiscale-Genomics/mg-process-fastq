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
    #from pycompss.api.constraint import constraint
    from pycompss.api.api import compss_wait_on
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN
    from utils.dummy_pycompss import task
    # from utils.dummy_pycompss import constraint
    from utils.dummy_pycompss import compss_wait_on

from basic_modules.tool import Tool

from pytadbit.parsers.genome_parser import parse_fasta
from pytadbit.parsers.map_parser import parse_map
from pytadbit.mapping import get_intersection

# ------------------------------------------------------------------------------

class tbParseMappingTool(Tool):
    """
    Tool for parsing the mapped reads and generating the list of paired ends
    that have a match at both ends.
    """

    def __init__(self):
        """
        Init function
        """
        print("TADbit parse mapping")
        Tool.__init__(self)

    @task(
        genome_seq=IN, enzyme_name=IN,
        window1_1=FILE_IN, window1_2=FILE_IN, window1_3=FILE_IN, window1_4=FILE_IN,
        window2_1=FILE_IN, window2_2=FILE_IN, window2_3=FILE_IN, window2_4=FILE_IN,
        reads=FILE_OUT)
    # @constraint(ProcessorCoreCount=32)
    def tb_parse_mapping_iter(
            self, genome_seq, enzyme_name,
            window1_1, window1_2, window1_3, window1_4,
            window2_1, window2_2, window2_3, window2_4,
            reads):
        """
        Function to map the aligned reads and return the matching pairs

        Parameters
        ----------
        genome_seq : dict
            Object containing the sequence of each of the chromosomes
        enzyme_name : str
            Name of the enzyme used to digest the genome
        window1_1 : str
            Location of the first window index file
        window1_2 : str
            Location of the second window index file
        window1_3 : str
            Location of the third window index file
        window1_4 : str
            Location of the fourth window index file
        window2_1 : str
            Location of the first window index file
        window2_2 : str
            Location of the second window index file
        window2_3 : str
            Location of the third window index file
        window2_4 : str
            Location of the fourth window index file
        reads : str
            Location of the reads thats that has a matching location at both
            ends of the paired reads


        Returns
        -------
        reads : str
            Location of the intersection of mapped reads that have matching
            reads in both pair end files

        """

        reads1 = reads + '_reads_1.tsv'
        reads2 = reads + '_reads_2.tsv'
        reads_both = reads + '_reads_both.tsv'

        parse_map(
            [window1_1, window1_2, window1_3, window1_4],
            [window2_1, window2_2, window2_3, window2_4],
            out_file1=reads1,
            out_file2=reads2,
            genome_seq=genome_seq,
            re_name=enzyme_name,
            verbose=True,
            #ncpus=32
        )

        intersections = get_intersection(reads1, reads2, reads_both, verbose=True)

        with open(reads, "wb") as f_out:
            with open(reads_both, "rb") as f_in:
                f_out.write(f_in.read())

        return True

    @task(
        genome_seq=IN, enzyme_name=IN,
        window1_full=FILE_IN, window1_frag=FILE_IN,
        window2_full=FILE_IN, window2_frag=FILE_IN,
        reads=FILE_OUT)
    # @constraint(ProcessorCoreCount=32)
    def tb_parse_mapping_frag(
            self, genome_seq, enzyme_name,
            window1_full, window1_frag, window2_full, window2_frag,
            reads):
        """
        Function to map the aligned reads and return the matching pairs

        Parameters
        ----------
        genome_seq : dict
            Object containing the sequence of each of the chromosomes
        enzyme_name : str
            Name of the enzyme used to digest the genome
        window1_full : str
            Location of the first window index file
        window1_frag : str
            Location of the second window index file
        window2_full : str
            Location of the first window index file
        window2_frag : str
            Location of the second window index file
        reads : str
            Location of the reads thats that has a matching location at both
            ends of the paired reads


        Returns
        -------
        reads : str
            Location of the intersection of mapped reads that have matching
            reads in both pair end files

        """

        print("TB WINDOWS - full 1", window1_full)
        print("TB WINDOWS - frag 1", window1_frag)
        print("TB WINDOWS - full 2", window2_full)
        print("TB WINDOWS - frag 2", window2_frag)

        #root_name = reads.split("/")

        #reads1 = "/".join(root_name) + '/reads_1.tsv'
        #reads2 = "/".join(root_name) + '/reads_2.tsv'
        reads1 = reads + '_reads_1.tsv'
        reads2 = reads + '_reads_2.tsv'
        reads_both = reads + '_reads_both.tsv'

        parse_map(
            [window1_frag, window1_full],
            [window2_frag, window2_full],
            out_file1=reads1,
            out_file2=reads2,
            genome_seq=genome_seq,
            re_name=enzyme_name,
            verbose=True
        )

        intersections = get_intersection(reads1, reads2, reads_both, verbose=True)

        with open(reads, "wb") as f_out:
            with open(reads_both, "rb") as f_in:
                f_out.write(f_in.read())

        return True

    def run(self, input_files, output_files, metadata=None):
        """
        The main function to map the aligned reads and return the matching
        pairs. Parsing of the mappings can be either iterative of fragment
        based. If it is to be iteractive then the locations of 4 output file
        windows for each end of the paired end window need to be provided. If
        it is fragment based, then only 2 window locations need to be provided
        along within an enzyme name.

        Parameters
        ----------
        input_files : list
            genome_file : str
                Location of the genome FASTA file
            window1_1 : str
                Location of the first window index file
            window1_2 : str
                Location of the second window index file
            window1_3 : str
                [OPTIONAL] Location of the third window index file
            window1_4 : str
                [OPTIONAL] Location of the fourth window index file
            window2_1 : str
                Location of the first window index file
            window2_2 : str
                Location of the second window index file
            window2_3 : str
                [OPTIONAL] Location of the third window index file
            window2_4 : str
                [OPTIONAL] Location of the fourth window index file
        metadata : dict
            windows : list
                List of lists with the window sizes to be computed
            enzyme_name : str
                Restricture enzyme name
            mapping : list
                The mapping function used. The options are iter or frag.


        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : dict
            Dict of matching metadata dict objects

        Example
        -------

        Iterative:

        .. code-block:: python

            from tool import tb_parse_mapping

            genome_file = 'genome.fasta'

            root_name_1 = "/tmp/data/expt_source_1".split
            root_name_2 = "/tmp/data/expt_source_2".split
            windows = [[1,25], [1,50], [1,75], [1,100]]

            windows1 = []
            windows2 = []

            for w in windows:
                tail = "_full_" + w[0] + "-" + w[1] + ".map"
                windows1.append('/'.join(root_name_1) + tail)
                windows2.append('/'.join(root_name_2) + tail)

            files = [genome_file] + windows1 + windows2

            tpm = tb_parse_mapping.tb_parse_mapping()
            metadata = {'enzyme_name' : 'MboI', 'mapping' : ['iter', 'iter'], 'expt_name' = 'test'}
            tpm_files, tpm_meta = tpm.run(files, metadata)


        Fragment based mapping:

        .. code-block:: python

            from tool import tb_parse_mapping

            genome_file = 'genome.fasta'

            root_name_1 = "/tmp/data/expt_source_1".split
            root_name_2 = "/tmp/data/expt_source_2".split
            windows = [[1,100]]

            start = windows[0][0]
            end   = windows[0][1]

            window1_1 = '/'.join(root_name_1) + "_full_" + start + "-" + end + ".map"
            window1_2 = '/'.join(root_name_1) + "_frag_" + start + "-" + end + ".map"

            window2_1 = '/'.join(root_name_2) + "_full_" + start + "-" + end + ".map"
            window2_2 = '/'.join(root_name_2) + "_frag_" + start + "-" + end + ".map"

            files = [
                genome_file,
                window1_1, window1_2,
                window2_1, window2_2,
            ]

            tpm = tb_parse_mapping.tb_parse_mapping()
            metadata = {'enzyme_name' : 'MboI', 'mapping' : ['frag', 'frag'], 'expt_name' = 'test'}
            tpm_files, tpm_meta = tpm.run(files, metadata)

        """

        genome_file = input_files[0]

        enzyme_name = metadata['enzyme_name']
        mapping_list = metadata['mapping']
        expt_name = metadata['expt_name']

        root_name = input_files[1].split("/")

        reads = "/".join(root_name[0:-1]) + '/'

        genome_seq = parse_fasta(genome_file)

        chromosome_meta = []
        for k in genome_seq:
            chromosome_meta.append([k, len(genome_seq[k])])

        # input and output share most metadata
        output_metadata = {
            'chromosomes' : chromosome_meta
        }

        if mapping_list[0] == mapping_list[1]:
            if mapping_list[0] == 'iter':
                window1_1 = input_files[1]
                window1_2 = input_files[2]
                window1_3 = input_files[3]
                window1_4 = input_files[4]

                window2_1 = input_files[5]
                window2_2 = input_files[6]
                window2_3 = input_files[7]
                window2_4 = input_files[8]

                read_iter = reads + expt_name + '_iter.tsv'

                results = self.tb_parse_mapping_iter(
                    genome_seq, enzyme_name,
                    window1_1, window1_2, window1_3, window1_4,
                    window2_1, window2_2, window2_3, window2_4,
                    read_iter)
                results = compss_wait_on(results)
                return ([read_iter], output_metadata)

            elif mapping_list[0] == 'frag':
                window1_full = input_files[1]
                window1_frag = input_files[2]

                window2_full = input_files[3]
                window2_frag = input_files[4]

                read_frag = reads + expt_name + '_frag.tsv'

                results = self.tb_parse_mapping_frag(
                    genome_seq, enzyme_name,
                    window1_full, window1_frag,
                    window2_full, window2_frag,
                    read_frag)

                results = compss_wait_on(results)
                return ([read_frag], output_metadata)

            reads = None
            return ([reads], output_metadata)

# ------------------------------------------------------------------------------
