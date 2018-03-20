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

from tool.common import cd
from utils import logger


class alignerUtils(object):
    """
    Functions for downloading and processing N-seq FastQ files. Functions
    provided allow for the downloading and indexing of the genome assemblies.
    """

    def __init__(self):
        """
        Initialise the module
        """
        logger.info("Alignment Utils")

    @staticmethod
    def replaceENAHeader(file_path, file_out):
        """
        The ENA header has pipes in the header as part of the stable_id. This
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

    @staticmethod
    def gem_index_genome(genome_file, gem_file):
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

    @staticmethod
    def bowtie_index_genome(genome_file, index_name=None):
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

    @staticmethod
    def bwa_index_genome(genome_file):
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

           from tool.aligner_utils import alignerUtils
           au_handle = alignerUtils()

           indexes = au_handle.bwa_index_genome('/<data_dir>/human_GRCh38.fa.gz')
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

    @staticmethod
    def bowtie2_align_reads(
            genome_file, bam_loc, params, reads_file_1, reads_file_2=None):
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

        reads = ["-U", reads_file_1]
        if reads_file_2 is not None:
            reads = [
                "-1", reads_file_1,
                "-2", reads_file_2
            ]

        logger.info(genome_file)
        logger.info(' '.join(params))

        g_idx = genome_file.split("/")
        g_idx[-1] = g_idx[-1].replace(".fasta", "")

        cmd_aln = ' '.join([
            'bowtie2',
            '-p 4',
            '-x', '/'.join(g_idx),
            ' '.join(params),
        ] + reads)

        cmd_view = ' '.join([
            'samtools view',
            '-b',
            '-o', bam_loc,
            reads_file_1 + '.sam'
        ])

        try:
            with open(reads_file_1 + '.sam', "w") as f_out:
                logger.info("BOWTIE2 COMMAND: " + cmd_aln)
                process = subprocess.Popen(cmd_aln, shell=True, stdout=f_out)
                process.wait()
        except (IOError, OSError) as msg:
            logger.info("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, cmd_aln))
            return False

        try:
            logger.info("BOWTIE2 COMMAND: " + cmd_view)
            process = subprocess.Popen(cmd_view, shell=True)
            process.wait()
        except (IOError, OSError) as msg:
            logger.info("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, cmd_view))
            return False

        return True

    @staticmethod
    def bwa_aln_align_reads_single(genome_file, reads_file, bam_loc, params):
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

        cmd_aln = ' '.join([
            'bwa aln',
            '-q', '5',
            ' '.join(params),
            '-f', reads_file + '.sai',
            genome_file, reads_file
        ])

        cmd_samse = ' '.join([
            'bwa samse',
            '-f', reads_file + '.sam',
            genome_file, reads_file + '.sai', reads_file
        ])

        cmd_view = ' '.join([
            'samtools view',
            '-b',
            '-o', bam_loc,
            reads_file + '.sam'
        ])

        command_lines = [cmd_aln, cmd_samse, cmd_view]

        # print("BWA COMMAND LINES:", command_lines)
        try:
            for command_line in command_lines:
                logger.info("BWA ALN COMMAND: " + command_line)
                process = subprocess.Popen(
                    command_line, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                process.wait()
                proc_out, proc_err = process.communicate()
                logger.info("BWA ALN stdout" + proc_out)
                logger.info("BWA ALN stderr" + proc_err)
        except (IOError, OSError) as msg:
            logger.info("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, command_line))
            proc_out, proc_err = process.communicate()
            logger.fatal("BWA ALN stdout" + proc_out)
            logger.fatal("BWA ALN stderr" + proc_err)
            return False

        return True

    @staticmethod
    def bwa_aln_align_reads_paired(genome_file, reads_file_1, reads_file_2, bam_loc, params):
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

        cmd_aln_1 = ' '.join([
            'bwa aln',
            '-q', '5',
            ' '.join(params),
            '-f', reads_file_1 + '.sai',
            genome_file, reads_file_1
        ])

        cmd_aln_2 = ' '.join([
            'bwa aln',
            '-q', '5',
            ' '.join(params),
            '-f', reads_file_2 + '.sai',
            genome_file, reads_file_2
        ])


        cmd_samse = ' '.join([
            'bwa sampe',
            '-f', reads_file_1 + '.sam',
            genome_file,
            reads_file_1 + '.sai', reads_file_2 + '.sai',
            reads_file_1, reads_file_2
        ])

        cmd_view = ' '.join([
            'samtools view',
            '-b',
            '-o', bam_loc,
            reads_file_1 + '.sam'
        ])

        command_lines = [cmd_aln_1, cmd_aln_2, cmd_samse, cmd_view]

        try:
            for command_line in command_lines:
                logger.info("BWA ALN COMMAND: " + command_line)
                process = subprocess.Popen(command_line, shell=True)
                process.wait()
        except (IOError, OSError) as msg:
            logger.info("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, command_line))
            return False

        return True

    @staticmethod
    def bwa_mem_align_reads(
            genome_file, bam_loc, params, reads_file_1, reads_file_2=None):
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

        reads = [reads_file_1]
        if reads_file_2 is not None:
            reads.append(reads_file_2)

        cmd_aln = ' '.join([
            'bwa mem -t 4',
            ' '.join(params),
            genome_file
        ] + reads)

        cmd_view = ' '.join([
            'samtools view',
            '-b',
            '-o', bam_loc,
            reads_file_1 + '.sam'
        ])

        try:
            with open(reads_file_1 + '.sam', "w") as f_out:
                logger.info("BWA MEM COMMAND: " + cmd_aln)
                process = subprocess.Popen(cmd_aln, shell=True, stdout=f_out)
                process.wait()
        except (IOError, OSError) as msg:
            logger.info("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, cmd_aln))
            return False

        try:
            logger.info("BWA MEM COMMAND: " + cmd_view)
            process = subprocess.Popen(cmd_view, shell=True)
            process.wait()
        except (IOError, OSError) as msg:
            logger.info("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, cmd_view))
            return False

        return True
