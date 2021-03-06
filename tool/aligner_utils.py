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
import os
import os.path
import shutil
import tarfile

from utils import logger
from tool.common import cd
from tool.bam_utils import bamUtils


class alignerUtils(object):  # pylint: disable=invalid-name
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
    def gem_index_genome(genome_file, index_name=None):
        """
        Create an index of the genome FASTA file with GEM. These are saved
        alongside the assembly file.

        Parameters
        ----------
        genome_file : str
            Location of the assembly file in the file system

        """

        if not index_name:
            index_name = genome_file
        command_line = 'gem-indexer -i ' + genome_file + ' -o ' + index_name

        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        return True

    @staticmethod
    def bowtie_index_genome(genome_file):
        """
        Create an index of the genome FASTA file with Bowtie2. These are saved
        alongside the assembly file.

        Parameters
        ----------
        genome_file : str
            Location of the assembly file in the file system
        """
        file_name = os.path.split(genome_file)

        bt2_1_name = genome_file + ".1.bt2"
        bt2_2_name = genome_file + ".2.bt2"
        bt2_3_name = genome_file + ".3.bt2"
        bt2_4_name = genome_file + ".4.bt2"
        rev1_bt2_name = genome_file + ".rev.1.bt2"
        rev2_bt2_name = genome_file + ".rev.2.bt2"

        with cd(file_name[0]):
            command_line = 'bowtie2-build ' + genome_file + ' ' + genome_file
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()

        return (bt2_1_name, bt2_2_name, bt2_3_name, bt2_4_name, rev1_bt2_name, rev2_bt2_name)

    def bowtie2_untar_index(self, genome_name, tar_file,  # pylint: disable=too-many-arguments
                            bt2_1_file, bt2_2_file, bt2_3_file, bt2_4_file,
                            bt2_rev1_file, bt2_rev2_file):
        """
        Extracts the BWA index files from the genome index tar file.

        Parameters
        ----------
        genome_file_name : str
            Location string of the genome fasta file
        tar_file : str
            Location of the Bowtie2 index file
        bt2_1_file : str
            Location of the amb index file
        bt2_2_file : str
            Location of the ann index file
        bt2_3_file : str
            Location of the bwt index file
        bt2_4_file : str
            Location of the pac index file
        bt2_rev1_file : str
            Location of the sa index file
        bt2_rev2_file : str
            Location of the sa index file

        Returns
        -------
        bool
            Boolean indicating if the task was successful
        """
        index_files = {
            "1.bt2": bt2_1_file,
            "2.bt2": bt2_2_file,
            "3.bt2": bt2_3_file,
            "4.bt2": bt2_4_file,
            "rev.1.bt2": bt2_rev1_file,
            "rev.2.bt2": bt2_rev2_file,
        }

        return self._untar_index(genome_name, tar_file, index_files)

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

    def bwa_untar_index(self, genome_name, tar_file,  # pylint: disable=too-many-arguments
                        amb_file, ann_file, bwt_file, pac_file, sa_file):
        """
        Extracts the BWA index files from the genome index tar file.

        Parameters
        ----------
        genome_file_name : str
            Location string of the genome fasta file
        genome_idx : str
            Location of the BWA index file
        amb_file : str
            Location of the amb index file
        ann_file : str
            Location of the ann index file
        bwt_file : str
            Location of the bwt index file
        pac_file : str
            Location of the pac index file
        sa_file : str
            Location of the sa index file

        Returns
        -------
        bool
            Boolean indicating if the task was successful
        """
        index_files = {
            "amb": amb_file,
            "ann": ann_file,
            "bwt": bwt_file,
            "pac": pac_file,
            "sa": sa_file
        }

        return self._untar_index(genome_name, tar_file, index_files)

    @staticmethod
    def _untar_index(genome_name, tar_file, index_files):
        """
        Untar the specified files for a genomic index into the specified
        location.

        Parameters
        ----------
        genome_name : str
            Name of the genome for the folder within the tar file
        tar_file : str
            Location of the tarred index files
        index_files : dict
            Dictionary object of the suffix and final index file location
        """
        try:
            g_dir = os.path.split(tar_file)[0]

            tar = tarfile.open(tar_file)
            tar.extractall(path=g_dir)
            tar.close()

            gidx_folder = tar_file.replace('.tar.gz', '/') + genome_name

            piece_size = 5120000  # 500MB

            for suffix in list(index_files.keys()):
                with open(index_files[suffix], "wb") as f_out:
                    with open(gidx_folder + "." + suffix, "rb") as f_in:
                        while True:
                            piece = f_in.read(piece_size)

                            if not piece:
                                break  # end of file

                            print("PIECE:", piece[0:5])
                            f_out.write(piece)

            shutil.rmtree(tar_file.replace('.tar.gz', ''))
        except (OSError, IOError) as error:
            logger.fatal("UNTAR: I/O error({0}): {1}".format(error.errno, error.strerror))
            return False

        return True

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

        cmd_aln = ' '.join([
            'bowtie2',
            '-p 4',
            '-x', genome_file,
            '-S', reads_file_1 + '.sam',
            ' '.join(params),
        ] + reads)

        try:
            logger.info("BOWTIE2 COMMAND: " + cmd_aln)
            process = subprocess.Popen(cmd_aln, shell=True)
            process.wait()
        except (OSError, IOError) as msg:
            logger.info("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, cmd_aln))
            return False

        bu_handle = bamUtils()
        return_val = bu_handle.sam_to_bam(reads_file_1 + '.sam', bam_loc)

        if return_val:
            os.remove(reads_file_1 + '.sam')
        else:
            logger.warn("IO error with {} to {}".format(reads_file_1 + '.sam', bam_loc))

        return return_val

    def _bwa_aln_sai(self, genome_file, reads_file, params, single=True):  # pylint: disable=no-self-use
        """
        Generate the sai files required for creating the sam file.

        Parameters
        ----------
        genome_file : str
            Location of the assembly file in the file system
        reads_file : str
            Location of the reads file in the file system
        params : dict
            Dictionary of the parameters for bwa aln
        single : bool
            True for single ended, will use 4 threads for processing
            False for paired end, will use 2 threads for processing
        """
        threads = "2"
        if single:
            threads = "4"

        cmd_aln_sai = ' '.join([
            'bwa aln',
            '-t', threads,
            '-q', '5',
            ' '.join(params),
            '-f', reads_file + '.sai',
            genome_file, reads_file
        ])

        try:
            logger.info("BWA ALN COMMAND: " + cmd_aln_sai)
            process = subprocess.Popen(cmd_aln_sai, shell=True)
            process.wait()
        except (OSError, IOError) as msg:
            logger.info("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, cmd_aln_sai))

    def bwa_aln_align_reads_single(self, genome_file, reads_file, bam_loc, params):
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

        cmd_samse = ' '.join([
            'bwa samse',
            '-f', reads_file + '.sam',
            genome_file, reads_file + '.sai', reads_file
        ])

        self._bwa_aln_sai(genome_file, reads_file, params, True)

        try:
            logger.info("BWA ALN COMMAND: " + cmd_samse)
            process = subprocess.Popen(
                cmd_samse, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            process.wait()
            proc_out, proc_err = process.communicate()  # pylint: disable=unused-variable
        except (OSError, IOError) as msg:
            logger.info("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, cmd_samse))
            proc_out, proc_err = process.communicate()
            logger.fatal("BWA ALN stderr" + proc_err)
            return False

        bu_handle = bamUtils()
        return_val = bu_handle.sam_to_bam(reads_file + '.sam', bam_loc)

        if return_val:
            os.remove(reads_file + '.sam')
            os.remove(reads_file + '.sai')
        else:
            logger.warn("IO error with {} to {}".format(reads_file + '.sam', bam_loc))

        return return_val

    def bwa_aln_align_reads_paired(self, genome_file, reads_file_1, reads_file_2, bam_loc, params):  # pylint: disable=too-many-arguments
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

        cmd_samse = ' '.join([
            'bwa sampe',
            '-f', reads_file_1 + '.sam',
            genome_file,
            reads_file_1 + '.sai', reads_file_2 + '.sai',
            reads_file_1, reads_file_2
        ])

        try:
            import multiprocessing

            f1_proc = multiprocessing.Process(
                name='fastq_1', target=self._bwa_aln_sai,
                args=(genome_file, reads_file_1, params, False)
            )
            f2_proc = multiprocessing.Process(
                name='fastq_2', target=self._bwa_aln_sai,
                args=(genome_file, reads_file_2, params, False)
            )

            f1_proc.start()
            f2_proc.start()

            f1_proc.join()
            f2_proc.join()
        except (OSError, IOError) as msg:
            logger.info("SAI ERROR: I/O error({0}): {1}".format(
                msg.errno, msg.strerror))
            return False

        try:
            logger.info("BWA ALN COMMAND: " + cmd_samse)
            process = subprocess.Popen(cmd_samse, shell=True)
            process.wait()
        except (OSError, IOError) as msg:
            logger.info("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, cmd_samse))
            return False

        bu_handle = bamUtils()
        return_val = bu_handle.sam_to_bam(reads_file_1 + '.sam', bam_loc)

        if return_val:
            os.remove(reads_file_1 + '.sam')
            os.remove(reads_file_1 + '.sai')
            os.remove(reads_file_2 + '.sai')
        else:
            logger.warn("IO error with {} to {}".format(reads_file_1 + '.sam', bam_loc))

        return return_val

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

        try:
            with open(reads_file_1 + '.sam', "w") as f_out:
                logger.info("BWA MEM COMMAND: " + cmd_aln)
                process = subprocess.Popen(cmd_aln, shell=True, stdout=f_out)
                process.wait()
        except (OSError, IOError) as msg:
            logger.info("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, cmd_aln))
            return False

        bu_handle = bamUtils()
        return_val = bu_handle.sam_to_bam(reads_file_1 + '.sam', bam_loc)

        if return_val:
            os.remove(reads_file_1 + '.sam')
        else:
            logger.warn("IO error with {} to {}".format(reads_file_1 + '.sam', bam_loc))

        return return_val
