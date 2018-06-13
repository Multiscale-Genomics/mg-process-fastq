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
import subprocess
import pysam

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, FILE_INOUT
    from pycompss.api.task import task
    from pycompss.api.api import barrier, compss_delete_file
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, FILE_INOUT  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import barrier, compss_delete_file  # pylint: disable=ungrouped-imports


# ------------------------------------------------------------------------------

class bamUtils(object):
    """
    Tool for handling bam files
    """

    def __init__(self):
        """
        Initialise the tool with its configuration.
        """
        logger.info("BAM Utils")

    @staticmethod
    def bam_sort(bam_file):
        """
        Wrapper for the pysam SAMtools sort function

        Parameters
        ----------
        bam_file : str
            Location of the bam file to sort
        """
        try:
            pysam.sort(  # pylint: disable=no-member
                "-o", bam_file, "-T", bam_file + "_sort", bam_file)
        except IOError:
            return False
        return True

    @staticmethod
    def bam_merge(*args):
        """
        Wrapper for the pysam SAMtools merge function

        Parameters
        ----------
        bam_file_1 : str
            Location of the bam file to merge into
        bam_file_2 : str
            Location of the bam file that is to get merged into bam_file_1
        """
        logger.info("Merging:")

        if isinstance(args[0], list):
            final_bam = args[0][0]
            tmp_bam = final_bam + "_merge.bam"
            pysam.merge("-f", tmp_bam, *args[0])  # pylint: disable=no-member
        else:
            final_bam = args[0]
            tmp_bam = final_bam + "_merge.bam"
            pysam.merge("-f", tmp_bam, *args)  # pylint: disable=no-member

        try:
            with open(tmp_bam, "rb") as f_in:
                with open(final_bam, "wb") as f_out:
                    f_out.write(f_in.read())
        except IOError:
            return False

        os.remove(tmp_bam)

        return True

    @staticmethod
    def bam_copy(bam_in, bam_out):
        """
        Wrapper function to copy from one bam file to another

        Parameters
        ----------
        bam_in : str
            Location of the input bam file
        bam_out : str
            Location of the output bam file
        """
        logger.info("Copying: " + bam_in + " - " + bam_out)

        try:
            with open(bam_in, "rb") as f_in:
                with open(bam_out, "wb") as f_out:
                    f_out.write(f_in.read())
        except IOError:
            return False

        return True

    @staticmethod
    def bam_index(bam_file, bam_idx_file):
        """
        Wrapper for the pysam SAMtools index function

        Parameters
        ----------
        bam_file : str
            Location of the bam file that is to be indexed
        bam_idx_file : str
            Location of the bam index file (.bai)
        """

        cmd_view = ' '.join([
            'samtools index',
            '-b',
            bam_file,
            bam_file + "_tmp.bai"
        ])

        try:
            logger.info("INDEX BAM COMMAND: " + cmd_view)
            process = subprocess.Popen(cmd_view, shell=True)
            process.wait()
        except (IOError, OSError) as msg:
            logger.info("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, cmd_view))
            return False

        try:
            with open(bam_file + "_tmp.bai", "rb") as f_in:
                with open(bam_idx_file, "wb") as f_out:
                    f_out.write(f_in.read())
        except IOError:
            return False

        return True

    @staticmethod
    def bam_list_chromosomes(bam_file):
        """
        Wrapper to list the chromosome names that are present within the bam file

        Parameters
        ----------
        bam_file : str
            Location of the bam file

        Returns
        -------
        list
            List of the names of the chromosomes that are present in the bam file
        """
        bam_file_handle = pysam.AlignmentFile(bam_file, "rb")  # pylint: disable=no-member
        if "SQ" not in bam_file_handle.header:
            return []

        chromosome_list = [chromosome["SN"] for chromosome in bam_file_handle.header["SQ"]]
        bam_file_handle.close()

        return chromosome_list

    @staticmethod
    def bam_split(bam_file_in, bai_file, chromosome, bam_file_out):  # pylint: disable=unused-argument
        """
        Wrapper to extract a single chromosomes worth of reading into a new bam
        file

        Parameters
        ----------
        bam_file_in : str
            Location of the input bam file
        bai_file : str
            Location of the bam index file. This needs to be in the same directory
            as the bam_file_in
        chromosome : str
            Name of the chromosome whose alignments are to be extracted
        bam_file_out : str
            Location of the output bam file
        """

        # Extract the subsection from the bam file using samtools
        cmd_view_1 = ' '.join([
            'samtools view',
            '-h',
            '-o', bam_file_in + ".sam",
            bam_file_in,
            chromosome
        ])

        #
        cmd_view_2 = ' '.join([
            'samtools view',
            '-b',
            '-o', bam_file_out,
            bam_file_in + ".sam"
        ])

        try:
            logger.info("EXTRACT SAM COMMAND: " + cmd_view_1)
            process = subprocess.Popen(cmd_view_1, shell=True)
            process.wait()
        except (IOError, OSError) as msg:
            logger.info("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, cmd_view_1))
            return False

        try:
            logger.info("CREATE BAM COMMAND: " + cmd_view_2)
            process = subprocess.Popen(cmd_view_2, shell=True)
            process.wait()
        except (IOError, OSError) as msg:
            logger.info("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, cmd_view_2))
            return False

        return True

    @staticmethod
    def bam_stats(bam_file):
        """
        Wrapper for the pysam SAMtools flagstat function

        Parameters
        ----------
        bam_file : str
            Location of the bam file

        Returns
        -------
        list : dict
            qc_passed : int
            qc_failed : int
            description : str
        """
        results = pysam.flagstat(bam_file)  # pylint: disable=no-member
        separate_results = results.strip().split("\n")
        return [
            {
                "qc_passed": int(element[0]),
                "qc_failed": int(element[2]),
                "description": " ".join(element[3:])
            } for element in [row.split(" ") for row in separate_results]
        ]

    @staticmethod
    def check_header(bam_file):
        """
        Wrapper for the pysam SAMtools for checking if a bam file is sorted

        Parameters
        ----------
        bool
            True if the file has been sorted
        """
        output = True

        bam_file_handle = pysam.AlignmentFile(bam_file, "rb")  # pylint: disable=no-member
        if ("SO" not in bam_file_handle.header["HD"] or
                bam_file_handle.header["HD"]["SO"] == "unsorted"):
            output = False
        bam_file_handle.close()

        return output


class bamUtilsTask(object):
    """
    Wrappers so that the function above can be used as part of a @task within
    COMPSs avoiding the files being copied around the infrastructure too many
    times
    """

    def __init__(self):
        """
        Init function
        """
        logger.info("BAM @task Utils")

    @task(returns=list, bam_file=FILE_IN)
    def bam_list_chromosomes(self, bam_file):  # pylint: disable=no-self-use
        """
        Wrapper to get the list of chromosomes in a given bam file

        Parameters
        ----------
        bam_file : str
            Location of the bam file

        Returns
        -------
        chromosome_list : list
            List of the chromosomes in the bam file
        """
        bam_handle = bamUtils()
        return bam_handle.bam_list_chromosomes(bam_file)

    @task(bam_file=FILE_INOUT)
    def bam_sort(self, bam_file):  # pylint: disable=no-self-use
        """
        Wrapper for the pysam SAMtools sort function

        Parameters
        ----------
        bam_file : str
            Location of the bam file to sort
        """
        bam_handle = bamUtils()
        return bam_handle.bam_sort(bam_file)

    def bam_merge(self, in_bam_job_files):  # pylint: disable=too-many-branches
        """
        Wrapper task taking any number of bam files and merging them into a
        single bam file.

        Parameters
        ----------
        bam_job_files : list
            List of the locations of the separate bam files that are to be merged
            The first file in the list will be taken as the output file name
        """
        merge_round = -1

        if len(in_bam_job_files) == 1:
            return in_bam_job_files[0]

        bam_job_files = [i for i in in_bam_job_files]

        cleanup_files = []

        while True:
            merge_round += 1
            if len(bam_job_files) > 1:
                tmp_alignments = []

                if bam_job_files:
                    while len(bam_job_files) >= 10:
                        current_list_len = len(bam_job_files)
                        for i in range(0, current_list_len-9, 10):  # pylint: disable=unused-variable
                            bam_out = bam_job_files[0] + "_merge_" + str(merge_round) + ".bam"
                            tmp_alignments.append(bam_out)
                            cleanup_files.append(bam_out)

                            self.bam_merge_10(
                                bam_job_files.pop(0), bam_job_files.pop(0), bam_job_files.pop(0),
                                bam_job_files.pop(0), bam_job_files.pop(0), bam_job_files.pop(0),
                                bam_job_files.pop(0), bam_job_files.pop(0), bam_job_files.pop(0),
                                bam_job_files.pop(0), bam_out
                            )

                    bam_out = bam_job_files[0] + "_merge_" + str(merge_round) + ".bam"
                    if len(bam_job_files) >= 5:
                        tmp_alignments.append(bam_out)
                        cleanup_files.append(bam_out)
                        self.bam_merge_5(
                            bam_job_files.pop(0), bam_job_files.pop(0), bam_job_files.pop(0),
                            bam_job_files.pop(0), bam_job_files.pop(0), bam_out
                        )
                        bam_out = bam_job_files[0] + "_merge_" + str(merge_round) + ".bam"

                    if len(bam_job_files) == 4:
                        tmp_alignments.append(bam_out)
                        cleanup_files.append(bam_out)
                        self.bam_merge_4(
                            bam_job_files.pop(0), bam_job_files.pop(0), bam_job_files.pop(0),
                            bam_job_files.pop(0), bam_out
                        )
                    elif len(bam_job_files) == 3:
                        tmp_alignments.append(bam_out)
                        cleanup_files.append(bam_out)
                        self.bam_merge_3(
                            bam_job_files.pop(0), bam_job_files.pop(0), bam_job_files.pop(0),
                            bam_out
                        )
                    elif len(bam_job_files) == 2:
                        tmp_alignments.append(bam_out)
                        cleanup_files.append(bam_out)
                        self.bam_merge_2(
                            bam_job_files.pop(0), bam_job_files.pop(0), bam_out
                        )
                    else:
                        tmp_alignments.append(bam_job_files[0])

                barrier()

                bam_job_files = []
                bam_job_files = [new_bam for new_bam in tmp_alignments]

            else:
                break

        return_value = self.bam_copy(bam_job_files[0], in_bam_job_files[0])
        for tmp_bam_file in cleanup_files.append(bam_out):
            compss_delete_file(tmp_bam_file)

        return return_value

    @task(bam_file_1=FILE_IN, bam_file_2=FILE_IN, bam_file_out=FILE_OUT)
    def bam_merge_2(self, bam_file_1, bam_file_2, bam_file_out):  # pylint: disable=no-self-use
        """
        Wrapper for the pysam SAMtools merge function

        Parameters
        ----------
        bam_file_1 : str
            Location of the bam file to merge into
        bam_file_2 : str
            Location of the bam file that is to get merged into bam_file_1
        """
        bam_handle = bamUtils()
        bam_handle.bam_merge(bam_file_1, bam_file_2)
        bam_handle.bam_copy(bam_file_1, bam_file_out)

    @task(
        bam_file_1=FILE_IN, bam_file_2=FILE_IN, bam_file_3=FILE_IN,
        bam_file_out=FILE_OUT)
    def bam_merge_3(self, bam_file_1, bam_file_2, bam_file_3, bam_file_out):  # pylint: disable=no-self-use
        """
        Wrapper for the pysam SAMtools merge function

        Parameters
        ----------
        bam_file_1 : str
            Location of the bam file to merge into
        bam_file_2 : str
            Location of the bam file that is to get merged into bam_file_1
        bam_file_3 : str
            Location of the bam file that is to get merged into bam_file_1
        """
        bam_handle = bamUtils()
        bam_handle.bam_merge([
            bam_file_1, bam_file_2, bam_file_3
        ])
        return bam_handle.bam_copy(bam_file_1, bam_file_out)

    @task(
        bam_file_1=FILE_IN, bam_file_2=FILE_IN, bam_file_3=FILE_IN,
        bam_file_4=FILE_IN, bam_file_out=FILE_OUT)
    def bam_merge_4(self, bam_file_1, bam_file_2, bam_file_3, bam_file_4, bam_file_out):  # pylint: disable=no-self-use
        """
        Wrapper for the pysam SAMtools merge function

        Parameters
        ----------
        bam_file_1 : str
            Location of the bam file to merge into
        bam_file_2 : str
            Location of the bam file that is to get merged into bam_file_1
        bam_file_3 : str
            Location of the bam file that is to get merged into bam_file_1
        bam_file_4 : str
            Location of the bam file that is to get merged into bam_file_1
        """
        bam_handle = bamUtils()
        bam_handle.bam_merge([
            bam_file_1, bam_file_2, bam_file_3, bam_file_4
        ])
        return bam_handle.bam_copy(bam_file_1, bam_file_out)

    @task(
        bam_file_1=FILE_IN, bam_file_2=FILE_IN, bam_file_3=FILE_IN,
        bam_file_4=FILE_IN, bam_file_5=FILE_IN, bam_file_out=FILE_OUT
    )  # pylint: disable=no-self-use,too-many-arguments
    def bam_merge_5(self, bam_file_1, bam_file_2, bam_file_3, bam_file_4, bam_file_5, bam_file_out):
        """
        Wrapper for the pysam SAMtools merge function

        Parameters
        ----------
        bam_file_1 : str
            Location of the bam file to merge into
        bam_file_2 : str
            Location of the bam file that is to get merged into bam_file_1
        bam_file_3 : str
            Location of the bam file that is to get merged into bam_file_1
        bam_file_4 : str
            Location of the bam file that is to get merged into bam_file_1
        bam_file_5 : str
            Location of the bam file that is to get merged into bam_file_1
        """
        bam_handle = bamUtils()
        bam_handle.bam_merge([
            bam_file_1, bam_file_2, bam_file_3, bam_file_4, bam_file_5
        ])
        return bam_handle.bam_copy(bam_file_1, bam_file_out)

    @task(
        bam_file_1=FILE_IN, bam_file_2=FILE_IN, bam_file_3=FILE_IN,
        bam_file_4=FILE_IN, bam_file_5=FILE_IN, bam_file_6=FILE_IN,
        bam_file_7=FILE_IN, bam_file_8=FILE_IN, bam_file_9=FILE_IN,
        bam_file_10=FILE_IN, bam_file_out=FILE_OUT
    )  # pylint: disable=no-self-use,too-many-arguments
    def bam_merge_10(self,
                     bam_file_1, bam_file_2, bam_file_3, bam_file_4, bam_file_5,
                     bam_file_6, bam_file_7, bam_file_8, bam_file_9, bam_file_10,
                     bam_file_out):
        """
        Wrapper for the pysam SAMtools merge function

        Parameters
        ----------
        bam_file_1 : str
            Location of the bam file to merge into
        bam_file_2 : str
            Location of the bam file that is to get merged into bam_file_1
        bam_file_3 : str
            Location of the bam file that is to get merged into bam_file_1
        bam_file_4 : str
            Location of the bam file that is to get merged into bam_file_1
        bam_file_5 : str
            Location of the bam file that is to get merged into bam_file_1
        bam_file_6 : str
            Location of the bam file that is to get merged into bam_file_1
        bam_file_7 : str
            Location of the bam file that is to get merged into bam_file_1
        bam_file_8 : str
            Location of the bam file that is to get merged into bam_file_1
        bam_file_9 : str
            Location of the bam file that is to get merged into bam_file_1
        bam_file_10 : str
            Location of the bam file that is to get merged into bam_file_1
        """
        bam_handle = bamUtils()
        bam_handle.bam_merge([
            bam_file_1, bam_file_2, bam_file_3, bam_file_4, bam_file_5,
            bam_file_6, bam_file_7, bam_file_8, bam_file_9, bam_file_10
        ])
        return bam_handle.bam_copy(bam_file_1, bam_file_out)

    @task(bam_in=FILE_IN, bam_out=FILE_OUT)
    def bam_copy(self, bam_in, bam_out):  # pylint: disable=no-self-use
        """
        Wrapper function to copy from one bam file to another

        Parameters
        ----------
        bam_in : str
            Location of the input bam file
        bam_out : str
            Location of the output bam file
        """
        bam_handle = bamUtils()
        return bam_handle.bam_copy(bam_in, bam_out)

    @task(bam_file=FILE_IN, bam_idx_file=FILE_OUT)
    def bam_index(self, bam_file, bam_idx_file):  # pylint: disable=no-self-use
        """
        Wrapper for the pysam SAMtools merge function

        Parameters
        ----------
        bam_file : str
            Location of the bam file that is to be indexed
        bam_idx_file : str
            Location of the bam index file (.bai)
        """
        bam_handle = bamUtils()
        return bam_handle.bam_index(bam_file, bam_idx_file)

    @task(bam_file=FILE_IN)
    def bam_stats(self, bam_file):  # pylint: disable=no-self-use
        """
        Wrapper for the pysam SAMtools flagstat function

        Parameters
        ----------
        bam_file : str
            Location of the bam file that is to be indexed
        bam_idx_file : str
            Location of the bam index file (.bai)
        """
        bam_handle = bamUtils()
        return bam_handle.bam_stats(bam_file)

    @task(bam_file=FILE_IN)
    def check_header(self, bam_file):  # pylint: disable=no-self-use
        """
        Wrapper for the pysam SAMtools merge function

        Parameters
        ----------
        bam_file_1 : str
            Location of the bam file to merge into
        bam_file_2 : str
            Location of the bam file that is to get merged into bam_file_1
        """
        bam_handle = bamUtils()
        return bam_handle.check_header(bam_file)
