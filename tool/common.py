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

import gzip
import shutil
import shlex
import subprocess
import os.path

try:
    import pysam
except ImportError:
    print("[Error] Cannot import \"pysam\" package. Have you installed it?")
    exit(-1)

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

    def run_indexers(self, file_name):
        """
        Assembly Index Manager

        Manges the creation of indexes for a given genome assembly file. If the
        downloaded file has not been unzipped then it will get unzipped here.
        There are then 3 indexers that are available including BWA, Bowtie2 and
        GEM. If the indexes already exist for the given file then the indexing
        is not rerun.

        Parameters
        ----------
        file_name : str
            Location of the assembly FASTA file

        Returns
        -------
        dict
            bowtie : str
                Location of the Bowtie index file
            bwa : str
                Location of the BWA index file
            gem : str
                Location of the gem index file

        Example
        -------
        .. code-block:: python
           :linenos:

           from tool.common import common
           cf = common()

           indexes = cf.run_indexers('/<data_dir>/human_GRCh38.fa.gz')
           print(indexes)


        """

        file_name_unzipped = file_name.replace('.fa.gz', '.fa')
        file_name_nofa = file_name_unzipped.replace('.fa', '')

        if os.path.isfile(file_name_unzipped) is False:
            print("Unzipping")
            with gzip.open(file_name, 'rb') as f_in, open(file_name_unzipped, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        print("Indexing - BWA")
        amn_file, ann_file, bwt_file, pac_file, sa_file = self.bwa_index_genome(file_name)
        bwa_indexes = {
            'amn' : amn_file,
            'ann' : ann_file,
            'bwt' : bwt_file,
            'pac' : pac_file,
            'sa'  : sa_file
        }

        if os.path.isfile(file_name_nofa + '.1.bt2') is False:
            print("Indexing - Bowtie")
            self.bowtie_index_genome(file_name_unzipped)

        if os.path.isfile(file_name_nofa + '.gem') is False:
            print("Indexing - GEM")
            self.gem_index_genome(file_name_unzipped, file_name_nofa + '.gem')

        return {
            'bowtie' : file_name_unzipped + '.1.bt2',
            'bwa' : bwa_indexes + '.bwt',
            'gem' : file_name_unzipped + '.gem'
        }

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

    def bwa_align_reads(self, genome_file, reads_file, bam_loc):
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

    def merge_bam(self, data_dir, project_id, final_id, run_ids=None):
        """
        Merge together all the bams in a directory and sort to create the final
        bam ready to be filtered

        If run_ids is blank then the function looks for all bam files in the
        data_dir
        """
        out_bam_file = data_dir + project_id + '/' + final_id + '.bam'

        if len(run_ids) is None:
            bam_files = [f for f in os.path.listdir(data_dir + project_id) if f.endswith(("sai"))]
        else:
            bam_files = [f + ".bam" for f in run_ids]

        bam_sort_files = []
        bam_merge_files = []
        for bam in bam_files:
            bam_loc = data_dir + project_id + '/' + bam
            bam_sort_files.append(bam_loc)
            bam_merge_files.append(bam_loc)

        for bam_sort_file in bam_sort_files:
            print(bam_sort_file)
            pysam.sort("-o", str(bam_sort_file), str(bam_sort_file))

        if len(bam_sort_files) == 1:
            pysam.sort("-o", str(out_bam_file), str(bam_sort_files[0]))
        else:
            pysam.merge(out_bam_file, *bam_merge_files)
            pysam.sort(
                "-o", str(out_bam_file),
                "-T", str(out_bam_file) + ".bam_sort",
                str(out_bam_file)
            )

        pysam.index(str(out_bam_file))
