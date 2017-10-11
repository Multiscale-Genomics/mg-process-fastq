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
import sys
import tarfile

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT
    from dummy_pycompss import task
    from dummy_pycompss import compss_wait_on

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

from tool.common import common

# ------------------------------------------------------------------------------


class bwaAlignerTool(Tool):
    """
    Tool for aligning sequence reads to a genome using BWA
    """

    def __init__(self, configuration=None):
        """
        Init function
        """
        print("BWA Aligner")
        Tool.__init__(self)

    @task(returns=bool, genome_file_loc=FILE_IN, read_file_loc=FILE_IN,
          bam_loc=FILE_OUT, genome_idx=FILE_IN, isModifier=False)
    def bwa_aligner(  # pylint: disable=too-many-arguments
            self, genome_file_loc, read_file_loc, bam_loc, genome_idx):  # pylint: disable=unused-argument
        """
        BWA Aligner

        Parameters
        ----------
        genome_file_loc : str
            Location of the genomic fasta
        read_file_loc : str
            Location of the FASTQ file

        Returns
        -------
        bam_loc : str
            Location of the output file
        """
        g_dir = genome_idx.split("/")
        g_dir = "/".join(g_dir[:-1])

        print("genome_idx:", genome_idx)
        tar = tarfile.open(genome_idx)
        tar.extractall(path=g_dir)
        tar.close()

        gfl = genome_file_loc.split("/")
        genome_fa_ln = genome_idx.replace('.tar.gz', '/') + gfl[-1]
        os.symlink(genome_file_loc, genome_fa_ln)

        out_bam = read_file_loc + '.out.bam'
        common_handle = common()
        common_handle.bwa_align_reads(genome_fa_ln, read_file_loc, out_bam)

        with open(bam_loc, "wb") as f_out:
            with open(out_bam, "rb") as f_in:
                f_out.write(f_in.read())

        return True

    def run(self, input_files, metadata, output_files):
        """
        The main function to align bam files to a genome using BWA

        Parameters
        ----------
        input_files : dict
            File 0 is the genome file location, file 1 is the FASTQ file
        metadata : dict
        output_files : dict

        Returns
        -------
        output_files : dict
            First element is a list of output_bam_files, second element is the
            matching meta data
        output_metadata : dict
        """

        print("BWA ALIGNER (before):", type(output_files), output_files["output"])

        results = self.bwa_aligner(
            str(input_files["genome"]), str(input_files["loc"]), str(output_files["output"]),
            str(input_files["index"]))

        results = compss_wait_on(results)

        print("BWA ALIGNER:", os.path.isfile(output_files["output"]))

        # print("BWA ALIGNER - METADATA:", metadata)

        bam_meta = Metadata(
            metadata['loc'].data_type, "bam", output_files["output"],
            [metadata['genome'].file_path, metadata['loc'].file_path],
            {
                'assembly': metadata['genome'].meta_data['assembly'],
                "tool": "bwa_aligner"
            }
        )

        # print("BWA ALIGNER - METADATA:", bam_meta)

        return ({"bam": output_files["output"]}, {"bam": bam_meta})

# ------------------------------------------------------------------------------
