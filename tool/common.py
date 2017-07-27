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

import gzip
import shutil
import shlex
import subprocess
import os.path

try:
    from urllib2 import urlopen
    from urllib2 import Request
except ImportError:
    from urllib.request import urlopen
    from urllib.request import Request

from socket import error as SocketError
import errno

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
    Functions for downloading and processing \*-seq FastQ files. Functions
    provided allow for the downloading andindexing of the genome assemblies.
    """

    def __init__ (self):
        """
        Initialise the module
        """
        print("Common functions")


    def getGenomeFile(self, file_loc, species, assembly, index=True):
        """
        Function for downloading and extracting the DNA files from the ensembl
        FTP

        Parameters
        ----------
        file_loc : str
            Location of the zipped genome file
        species : str
            Species name
        assembly : str
            Genome assembly for the given species
        index : bool
            Flag to indicate if the indexes are to be built for the downloaded
            assembly (default : True)

        Returns
        -------
        dict
            zipped : str
                This is the location of the
        """

        file_loc_unzipped = file_loc.replace('.fa.gz', '.fa')

        if os.path.isfile(file_loc) is False and os.path.isfile(file_loc_unzipped) is False:
            ftp_url = 'ftp://ftp.ensembl.org/pub/current_fasta/'
            ftp_url += species.lower() + '/dna/' + species[0].upper() + species[1:] + '.'
            ftp_url += assembly + '.dna.toplevel.fa.gz'

            self.download_file(file_loc, ftp_url)

        if os.path.isfile(file_loc_unzipped) is False:
            with gzip.open(file_loc, 'rb') as f_in, open(file_loc_unzipped, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

#        indexes = {}
        if index is True:
            indexes = self.run_indexers(file_loc_unzipped)
            #bowtie_index = indexes['bowtie']
            #bwa_index = indexes['bwa']
            #gem_index = indexes['gem']

        #return {'zipped': file_loc, 'unzipped': file_loc_unzipped, 'index' : indexes}
        return True


    def getGenomeFromENA(self, data_dir, species, assembly, index = True):
        """
        Obtain the Assembly from the ENA. This means that it is possible to get
        the exact assembly matching that defined by the user.
        """

        file_name = data_dir + species + '_' + assembly + '/' + species + '.' + assembly + '.fa'

        if os.path.isfile(file_name) is False:
            ftp_list_url = 'ftp://ftp.ebi.ac.uk/pub/databases/ena/assembly/' + assembly[0:7] + '/' + assembly[0:10] + '/' + assembly + '_sequence_report.txt'
            res_list = urlopen(ftp_list_url)
            table = res_list.read()
            table = table.split("\n")

            chr_list = []
            for row in table:
                col = row.split("\t")
                if len(col) > 5 and col[3] == 'assembled-molecule':
                    chr_list.append(col[0])

            ftp_url = 'http://www.ebi.ac.uk/ena/data/view/' + ','.join(chr_list) + '&display=fasta'
            self.download_file(file_name, ftp_url)

        indexes = {}
        if index == True:
            self.replaceENAHeader(file_name)
            indexes = self.run_indexers(file_name)

        return {'unzipped': file_name, 'index' : indexes}


    def replaceENAHeader(self, file_path, file_out):
        """
        The ENA header has pipes in the header as part of teh stable_id. This
        function removes the ENA stable_id and replaces it with the final
        section after splitting the stable ID on the pipe.
        """
        #from tempfile import mkstemp
        #from shutil import move
        #from os import remove, close

        #Create temp file
        #file_handle, abs_path = mkstemp()
        with open(file_out, 'w') as new_file:
            with open(file_path) as old_file:
                for line in old_file:
                    if line[0] == '>':
                        space_line = line.split(" ")
                        new_file.write(">" + space_line[0].split("|")[-1].replace(">", "") + "\n")
                    else:
                        new_file.write(line)
        #close(file_handle)
        #Remove original file
        #remove(file_path)
        #Move new file
        #move(abs_path, file_path)

        return True


    def getcDNAFiles(self, data_dir, species, assembly, e_release):
        """
        Function for downloading and extracting the CDNA files from the ensembl FTP
        """

        file_name = data_dir + species + '_' + assembly + '/' + species + '.' + assembly + '.cdna.all.fa.gz'

        if os.path.isfile(file_name) is False:
            ftp_url = 'ftp://ftp.ensembl.org/pub/release-' + e_release + '/' + species.lower() + '/cdna/' + species[0].upper() + species[1:] + '.' + assembly + '.cdna.all.fa.gz'

            self.download_file(file_name, ftp_url)

        return file_name


    def getFastqFiles(self, ena_err_id, data_dir, ena_srr_id = None):
        """
        Function for downloading and extracting the FastQ files from the ENA

        Parameters
        ----------
        ena_err_id : str
        data_dir : str
        ena_srr_id : str

        Returns
        -------
        files_list : list
            Locations of the downloaded FASTQ files
        """

        f_index = urlopen(
            'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=' + str(ena_err_id) + '&result=read_run&fields=study_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp&download=txt'
        )
        data = f_index.read()
        rows = data.split("\n")
        row_count = 0
        files = []
        gzfiles = []
        for row in rows:
            if row_count == 0:
                row_count += 1
                continue

            row = row.rstrip()
            row = row.split("\t")

            if len(row) < 6:
                continue

            project = row[0]
            srr_id = row[1]
            if ena_srr_id != None and srr_id != ena_srr_id:
                continue
            fastq_files = row[6].split(';')
            row_count += 1

            for fastq_file in fastq_files:
                file_name = fastq_file.split("/")

                if (
                        os.path.isfile(data_dir + '/' + project + "/" + file_name[-1]) is False and
                        os.path.isfile(
                            data_dir + '/' + project + "/" + file_name[-1].replace(
                                '.fastq.gz', '.fastq'
                            )
                        ) is False
                    ):
                    file_location = data_dir + '/' + project + "/" + file_name[-1]
                    ftp_url = "ftp://" + fastq_file
                    self.download_file(file_location, ftp_url)

                if os.path.isfile(data_dir + '/' + project + "/" + file_name[-1]) is True:
                    gzfiles.append(data_dir + '/' + project + "/" + file_name[-1])
                files.append(
                    data_dir + '/' + project + "/" + file_name[-1].replace('.fastq.gz', '.fastq')
                )

        for gzf in gzfiles:
            with gzip.open(gzf, 'rb') as f_in, open(gzf.replace('.fastq.gz', '.fastq'), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(gzf)

        return files


    def download_file(self, file_location, url):
        """
        Function to download a file to a given location and file name. Will
        attempt to restart the download up to 5 times if the connection is
        closed by the remote server.
        """

        restart_counter = 0

        while True:
            meta = urlopen(url)
            info = meta.info()
            if restart_counter > 0:
                byte_start = int(os.stat(file_location).st_size) + 1
                request = Request(url, headers={"Range" : "bytes=" + str(byte_start) + "-" + str(info["content-length"])})
            else:
                request = Request(url)

            req = urlopen(request)
            chunk_size = 16 * 1024

            try:
                if restart_counter == 0:
                    with open(file_location, 'wb') as file_handle:
                        while True:
                            chunk = req.read(chunk_size)
                            if not chunk:
                                break
                            file_handle.write(chunk)
                else:
                    with open(file_location, 'ab') as file_handle:
                        while True:
                            chunk = req.read(chunk_size)
                            if not chunk:
                                break
                            file_handle.write(chunk)

                req.close()
            except SocketError as err:
                if err.errno != errno.ECONNRESET:
                    raise
                print(err)
                print("Attempting to restart download ...")
                restart_counter += 1

            if restart_counter >= 5:
                return False
            break

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
            self.gem_index_genome(file_name_unzipped)

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
        file_name = genome_file.split("/")

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

        for command_line in command_lines:
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()

        return output_bam_file

    def merge_bam(self, data_dir, project_id, final_id, run_ids=[]):
        """
        Merge together all the bams in a directory and sort to create the final
        bam ready to be filtered

        If run_ids is blank then the function looks for all bam files in the
        data_dir
        """
        out_bam_file = data_dir + project_id + '/' + final_id + '.bam'

        if len(run_ids) == 0:
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
