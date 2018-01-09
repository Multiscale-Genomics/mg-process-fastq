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

import shlex
import subprocess
import os.path

class cd(object):
    """
    Context manager for changing the current working directory
    """

    def __init__(self, newpath):
        self.newpath = os.path.expanduser(newpath)

    def __enter__(self):
        self.savedpath = os.getcwd()
        os.chdir(self.newpath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedpath)

class common(object):
    """
    Functions for downloading and processing *-seq FastQ files. Functions
    provided allow for the downloading andindexing of the genome assemblies.
    """

    def __init__(self):
        """
        Initialise the module
        """
        print("Common functions")

    def replaceENAHeader(self, file_path, file_out):
        """
        The ENA header has pipes in the header as part of teh stable_id. This
        function removes the ENA stable_id and replaces it with the final
        section after splitting the stable ID on the pipe.
        """
        with open(file_out, 'w') as new_file:
            with open(file_path) as old_file:
                for line in old_file:
                    if line[0] == '>':
                        space_line = line.split(" ")
                        new_file.write(">" + space_line[0].split("|")[-1].replace(">", "") + "\n")
                    else:
                        new_file.write(line)

        return True

    def gem_index_genome(self, genome_file, gem_file):
        """
        Create an index of the genome FASTA file with GEM. These are saved
        alongside the assembly file.

        Parameters
        ----------
        genome_file : str
            Location of the assembly file in the file system

        """
        command_line = 'gem-indexer -i ' + genome_file + ' -o ' + gem_file

        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        return True

    def bowtie_index_genome(self, genome_file, index_name=None):
        """
        Create an index of the genome FASTA file with Bowtie2. These are saved
        alongside the assembly file.

        Parameters
        ----------
        genome_file : str
            Location of the assembly file in the file system
        """
        file_name = genome_file.split("/")

        output_file = index_name
        if output_file is None:
            output_file = file_name[-1].replace('.fa', '')

        with cd("/".join(file_name[0:-1])):
            command_line = 'bowtie2-build ' + genome_file + ' ' + output_file
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()

        return True

    def bwa_index_genome(self, genome_file):
        """
        Create an index of the genome FASTA file with BWA. These are saved
        alongside the assembly file. If the index has already been generated
        then the locations of the files are returned

        Parameters
        ----------
        genome_file : str
            Location of the assembly file in the file system

        Returns
        -------
        amb_file : str
            Location of the amb file
        ann_file : str
            Location of the ann file
        bwt_file : str
            Location of the bwt file
        pac_file : str
            Location of the pac file
        sa_file : str
            Location of the sa file

        Example
        -------
        .. code-block:: python
           :linenos:

           from tool.common import common
           cf = common()

           indexes = cf.bwa_index_genome('/<data_dir>/human_GRCh38.fa.gz')
           print(indexes)


        """
        command_line = 'bwa index ' + genome_file

        amb_name = genome_file + '.amb'
        ann_name = genome_file + '.ann'
        bwt_name = genome_file + '.bwt'
        pac_name = genome_file + '.pac'
        sa_name = genome_file + '.sa'

        if os.path.isfile(bwt_name) is False:
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()

        return (amb_name, ann_name, bwt_name, pac_name, sa_name)

    def bwa_aln_align_reads_single(self, genome_file, reads_file, bam_loc):
        """
        Map the reads to the genome using BWA.
        Parameters
        ----------
        genome_file : str
            Location of the assembly file in the file system
        reads_file : str
            Location of the reads file in the file system
        bam_loc : str
            Location of the output file
        """

        intermediate_file = reads_file + '.sai'
        intermediate_sam_file = reads_file + '.sam'
        output_bam_file = bam_loc

        command_lines = [
            'bwa aln -q 5 -f ' + intermediate_file + ' ' + genome_file + ' ' + reads_file,
            'bwa samse -f ' + intermediate_sam_file  + ' ' + genome_file + ' ' +
            intermediate_file + ' ' + reads_file,
            'samtools view -b -o ' + output_bam_file + ' ' + intermediate_sam_file
        ]

        print("BWA COMMAND LINES:", command_lines)

        for command_line in command_lines:
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()

        return output_bam_file

    def bwa_aln_align_reads_paired(self, genome_file, reads_file_1, reads_file_2, bam_loc):
        """
        Map the reads to the genome using BWA.
        Parameters
        ----------
        genome_file : str
            Location of the assembly file in the file system
        reads_file : str
            Location of the reads file in the file system
        bam_loc : str
            Location of the output file
        """

        intermediate_file = reads_file_1 + '.sai'
        intermediate_sam_file = reads_file_2 + '.sam'
        output_bam_file = bam_loc

        command_lines = [
            # 'bwa aln -q 5 -f ' + intermediate_file + ' ' + genome_file + ' ' + reads_file,
            # 'bwa samse -f ' + intermediate_sam_file  + ' ' + genome_file + ' ' +
            # intermediate_file + ' ' + reads_file,
            # 'samtools view -b -o ' + output_bam_file + ' ' + intermediate_sam_file
        ]

        print("BWA COMMAND LINES:", command_lines)

        for command_line in command_lines:
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()

        return output_bam_file

    def bwa_mem_align_reads_single(self, genome_file, reads_file, bam_loc):
        """
        Map the reads to the genome using BWA.

        Parameters
        ----------
        genome_file : str
            Location of the assembly file in the file system
        reads_file : str
            Location of the reads file in the file system
        bam_loc : str
            Location of the output file
        """

        intermediate_sam_file = reads_file + '.sam'
        output_bam_file = bam_loc

        command_line = 'bwa mem ' + genome_file + ' ' + reads_file
        args = shlex.split(command_line)
        with open(intermediate_sam_file, "w") as f_out:
            sub_proc = subprocess.Popen(args, stdout=f_out)
            sub_proc.wait()

        command_line = 'samtools view -b -o ' + output_bam_file + ' ' + intermediate_sam_file
        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        return output_bam_file

    def bwa_mem_align_reads_paired(self, genome_file, reads_file_1, reads_file_2, bam_loc):
        """
        Map the reads to the genome using BWA.

        Parameters
        ----------
        genome_file : str
            Location of the assembly file in the file system
        reads_file : str
            Location of the reads file in the file system
        bam_loc : str
            Location of the output file
        """

        intermediate_sam_file = reads_file_1 + '.sam'
        output_bam_file = bam_loc

        command_line = 'bwa mem ' + genome_file + ' ' + reads_file_1 + ' ' + reads_file_2
        args = shlex.split(command_line)
        with open(intermediate_sam_file, "w") as f_out:
            sub_proc = subprocess.Popen(args, stdout=f_out)
            sub_proc.wait()

        command_line = 'samtools view -b -o ' + output_bam_file + ' ' + intermediate_sam_file
        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        return output_bam_file
