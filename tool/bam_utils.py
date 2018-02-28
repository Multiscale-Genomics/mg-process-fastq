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
import pysam

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, FILE_INOUT
    from pycompss.api.task import task
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, FILE_INOUT # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task # pylint: disable=ungrouped-imports

# ------------------------------------------------------------------------------

class bamUtils(object):
    """
    Tool for aligning sequence reads to a genome using BWA
    """

    def __init__(self):
        """
        Initialise the tool with its configuration.


        Parameters
        ----------
        configuration : dict
            a dictionary containing parameters that define how the operation
            should be carried out, which are specific to each Tool.
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
            pysam.sort("-o", bam_file, "-T", bam_file + "_sort", bam_file)
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
            pysam.merge("-f", tmp_bam, *args[0])
        else:
            final_bam = args[0]
            tmp_bam = final_bam + "_merge.bam"
            pysam.merge("-f", tmp_bam, *args)

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
        Wrapper for the pysam SAMtools merge function

        Parameters
        ----------
        bam_file : str
            Location of the bam file that is to be indexed
        bam_idx_file : str
            Location of the bam index file (.bai)
        """
        pysam.index(bam_file, bam_file + "_tmp.bai")

        try:
            with open(bam_file + "_tmp.bai", "rb") as f_in:
                with open(bam_idx_file, "wb") as f_out:
                    f_out.write(f_in.read())
        except IOError:
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
        results = pysam.flagstat(bam_file)
        separate_results = results.strip().split("\n")
        return [
            {
                "qc_passed" : int(element[0]),
                "qc_failed" : int(element[2]),
                "description" : " ".join(element[3:])
            } for element in [row.split(" ") for row in separate_results]
        ]


    @staticmethod
    def check_header(bam_file):
        """
        Wrapper for the pysam SAMtools merge function

        Parameters
        ----------
        bam_file_1 : str
            Location of the bam file to merge into
        bam_file_2 : str
            Location of the bam file that is to get merged into bam_file_1
        """
        output = True

        bam_file_handle = pysam.AlignmentFile(bam_file, "rb")
        if ("SO" not in bam_file_handle.header["HD"] or
                bam_file_handle.header["HD"]["SO"] == "unsorted"):
            output = False
        bam_file_handle.close()

        return output


class bamUtilsTask(object):
    """
    Wrappers so that the function above can be used as part of a @task within
    COMPSs avoiding the files being copied around the infrstructure too many
    times
    """

    def __init__(self):
        """
        Init function
        """
        logger.info("BAM @task Utils")

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

    @task(bam_file_1=FILE_INOUT, bam_file_2=FILE_IN)
    def bam_merge(self, *args):  # pylint: disable=no-self-use
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
        return bam_handle.bam_merge(*args)

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
