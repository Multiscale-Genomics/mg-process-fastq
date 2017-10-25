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

import os
import shlex
import shutil
import subprocess
import sys
import tarfile

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
from basic_modules.metadata import Metadata

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

    @task(file_loc=FILE_IN, index_loc=FILE_OUT)
    def bowtie2_indexer(
            self, file_loc, index_loc): # pylint: disable=unused-argument
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
        common_handle.bowtie_index_genome(file_loc, file_name)

        # tar.gz the index
        print("BT - index_loc", index_loc, index_loc.replace('.tar.gz', ''))
        idx_out_pregz = index_loc.replace('.tar.gz', '.tar')

        index_dir = index_loc.replace('.tar.gz', '')
        os.mkdir(index_dir)

        idx_split = index_dir.split("/")

        shutil.move(file_name + ".1.bt2", index_dir)
        shutil.move(file_name + ".2.bt2", index_dir)
        shutil.move(file_name + ".3.bt2", index_dir)
        shutil.move(file_name + ".4.bt2", index_dir)
        shutil.move(file_name + ".rev.1.bt2", index_dir)
        shutil.move(file_name + ".rev.2.bt2", index_dir)

        index_folder = idx_split[-1]

        tar = tarfile.open(idx_out_pregz, "w")
        tar.add(index_dir, arcname=index_folder)
        tar.close()

        command_line = 'pigz ' + idx_out_pregz
        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        return True

    def run(self, input_files, metadata, output_files):
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

        results = self.bowtie2_indexer(
            input_files["genome"],
            output_files["index"]
        )

        output_metadata = {
            "index": Metadata(
                data_type="sequence_mapping_index_bowtie",
                file_type="TAR",
                file_path=output_files["index"],
                sources=[metadata["genome"].file_path],
                taxon_id=metadata["genome"].taxon_id,
                meta_data={
                    "assembly": metadata["genome"].meta_data["assembly"],
                    "tool": "bowtie_indexer"
                }
            )
        }

        results = compss_wait_on(results)

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
