#!/usr/bin/env python

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

import os
import errno
import re


class fastqreader(object):  # pylint: disable=too-many-instance-attributes,invalid-name
    """
    Module for reading single end and paired end FASTQ files
    """

    def __init__(self):
        """
        Initialise the module
        """
        self.fastq1 = ''
        self.fastq2 = ''

        self.f1_file = None
        self.f2_file = None

        self.f1_eof = False
        self.f2_eof = False

        self.f1_output_file_loc = None
        self.f2_output_file_loc = None

        self.f1_output_file = None
        self.f2_output_file = None

        self.output_tag = ''
        self.output_file_count = 0

        self.paired = False

    def openFastQ(self, file1, file2=None):  # pylint: disable=invalid-name
        """
        Create file handles for reading the FastQ files

        Parameters
        ----------
        file1 : str
            Location of the first FASTQ file
        file2 : str
            Location of a paired end FASTQ file.
        """
        self.fastq1 = file1
        self.f1_file = open(self.fastq1, "r")
        self.f1_eof = False

        if file2 is not None:
            self.fastq2 = file2
            self.f2_file = open(self.fastq2, "r")
            self.f2_eof = False
            self.paired = True

    def closeFastQ(self):  # pylint: disable=invalid-name
        """
        Close file handles for the FastQ files.
        """
        self.f1_file.close()

        if self.paired is True:
            self.f2_file.close()

    def eof(self, side=1):
        """
        Indicate if the end of the file has been reached

        Parameters
        ----------
        side : int
            1 or 2
        """
        if side == 1:
            return self.f1_eof
        elif side == 2:
            return self.f2_eof

        return "ERROR"

    def next(self, side=1):
        """
        Get the next read element for the specific FastQ file pair

        Parameters
        ----------
        side : int
            1 or 2 to get the element from the relevant end (DEFAULT: 1)

        Returns
        -------
        dict
            id : str
                Sequence ID
            seq : str
                Called sequence
            add : str
                Plus sign
            score : str
                Base call score
        """
        read_id = ''
        read_seq = ''
        read_addition = ''
        read_score = ''

        if side == 1:
            try:
                start_posn = self.f1_file.tell()
                read_id = self.f1_file.readline()
                read_seq = self.f1_file.readline()
                read_addition = self.f1_file.readline()
                read_score = self.f1_file.readline()
                end_posn = self.f1_file.tell()
                if read_id == "":
                    raise EOFError
            except EOFError:
                self.f1_eof = True
                return False

        elif side == 2:
            try:
                start_posn = self.f2_file.tell()
                read_id = self.f2_file.readline()
                read_seq = self.f2_file.readline()
                read_addition = self.f2_file.readline()
                read_score = self.f2_file.readline()
                end_posn = self.f2_file.tell()
                if read_id == "":
                    raise EOFError
            except EOFError:
                self.f2_eof = True
                return False
        else:
            return 'ERROR'

        return {
            'id': read_id.rstrip(),
            'seq': read_seq,
            'add': read_addition,
            'score': read_score,
            'start_posn': start_posn,
            'end_posn': end_posn
        }

    def createOutputFiles(self, tag=''):  # pylint: disable=invalid-name
        """
        Create and open the file handles for the output files

        Parameters
        ----------
        tag : str
            Tag to identify the output files (DEFAULT: '')
        """
        if tag != '' and self.output_tag != tag:
            self.output_tag = tag

        fq1 = self.fastq1.split("/")
        new_suffix = "." + str(self.output_tag) + "_" + str(self.output_file_count) + ".fastq"
        fq1[-1] = re.sub('.fastq$', new_suffix, fq1[-1])
        fq1.insert(-1, "tmp")

        if os.path.isdir("/".join(fq1[0:-1])) is False:
            try:
                os.mkdir("/".join(fq1[0:-1]))
            except OSError as oserror:
                if oserror.errno != errno.EEXIST:
                    raise OSError

        self.f1_output_file = open("/".join(fq1), "w")
        self.f1_output_file_loc = "/".join(fq1)

        if self.paired is True:
            fq2 = self.fastq2.split("/")
            new_suffix = "." + str(self.output_tag) + "_" + str(self.output_file_count) + ".fastq"
            fq2[-1] = re.sub('.fastq$', new_suffix, fq2[-1])
            fq2.insert(-1, "tmp")
            self.f2_output_file = open("/".join(fq2), "w")
            self.f2_output_file_loc = "/".join(fq2)

    def writeOutput(self, read, side=1):  # pylint: disable=invalid-name
        """
        Writer to print the extracted lines

        Parameters
        ----------
        read : dict
            Read is the dictionary object returned from self.next()
        side : int
            The side that the read has coe from (DEFAULT: 1)

        Returns
        -------
        bool
            False if a value other than 1 or 2 is entered for the side.
        """
        line = read["id"] + "\n" + read["seq"] + read["add"] + read["score"]
        if side == 1:
            self.f1_output_file.write(line)
        elif side == 2:
            self.f2_output_file.write(line)
        else:
            return False
        return True

    def closeOutputFiles(self):  # pylint: disable=invalid-name
        """
        Close the output file handles
        """
        self.f1_output_file.close()

        if self.paired is True:
            self.f2_output_file.close()

    def incrementOutputFiles(self):  # pylint: disable=invalid-name
        """
        Increment the counter and create new files for splitting the original
        FastQ paired end files.
        """
        self.closeOutputFiles()

        self.output_file_count += 1

        self.createOutputFiles(self.output_tag)
