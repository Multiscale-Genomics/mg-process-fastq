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

from mg_process_fastq.tool.fastqreader import fastqreader


class fastqUtils(object):  # pylint: disable=invalid-name
    """
    Set of methods to help with the management of FastQ files.
    """

    @staticmethod
    def fastq_match_paired_ends(fastq_1, fastq_2):
        """
        Take 2 fastq files and remove ends that don't have a matching pair.
        Requires that the FastQ files are ordered correctly.

        Mismatches can occur if there is a filtering step that removes one of
        the paired ends

        Parameters
        ----------
        fastq_1 : str
            Location of FastQ file
        fastq_2 : str
            Location of FastQ file
        """
        fqr = fastqreader()
        fqr.openFastQ(fastq_1, fastq_2)
        fqr.createOutputFiles('match')

        record1 = fqr.next(1)
        record2 = fqr.next(2)

        count_r1 = 0
        count_r2 = 0
        count_r3 = 0

        while fqr.eof(1) is False and fqr.eof(2) is False:
            r1_id = record1["id"].split(" ")
            r2_id = record2["id"].split(" ")

            if r1_id[0] == r2_id[0]:
                fqr.writeOutput(record1, 1)
                fqr.writeOutput(record2, 2)

                record1 = fqr.next(1)
                record2 = fqr.next(2)

                count_r1 += 1
                count_r2 += 1
                count_r3 += 1
            elif r1_id[0] < r2_id[0]:
                record1 = fqr.next(1)
                count_r1 += 1
            else:
                record2 = fqr.next(2)
                count_r2 += 1

        fqr.closeFastQ()
        fqr.closeOutputFiles()

        with open(fqr.f1_output_file_loc, "rb") as f1_match:
            with open(fastq_1, "wb") as f_old:
                f_old.write(f1_match.read())

        with open(fqr.f2_output_file_loc, "rb") as f2_match:
            with open(fastq_2, "wb") as f_old:
                f_old.write(f2_match.read())

        os.remove(fqr.f1_output_file_loc)
        os.remove(fqr.f2_output_file_loc)

    @staticmethod
    def fastq_sort_file(fastq, output=None):
        """
        Sorting of a FastQ file

        Parameters
        ----------
        fastq : str
            Location of the FastQ file to sort
        output : str
            [OPTIONAL] Location of the output FastQ file. If left blank then the
            sorted output is saved to the same location as `fastq`
        """
        fqr = fastqreader()
        fqr.openFastQ(fastq)

        record = fqr.next()

        fq_tmp = {}
        while fqr.eof() is False:
            r_id = record["id"].split(" ")
            fq_tmp[r_id[0]] = [record["start_posn"], record["end_posn"]]
            record = fqr.next()

        fqr.closeFastQ()

        with open(fastq, "r") as fq_file:
            with open(fastq + ".sorted", "w") as tmp_fq_file:
                for r_id in sorted(fq_tmp.keys()):
                    fq_file.seek(fq_tmp[r_id][0])
                    tmp_fq_file.write(fq_file.read(fq_tmp[r_id][1]-fq_tmp[r_id][0]))

        output_fastq = fastq
        if output is not None:
            output_fastq = output
        with open(fastq + ".sorted", "rb") as fq_sorted:
            with open(output_fastq, "wb") as f_out:
                f_out.write(fq_sorted.read())

        if output is not None:
            os.remove(fastq + ".sorted")

    @staticmethod
    def fastq_randomise(fastq, output=None):
        """
        Randomising the order of reads withim a FastQ file

        Parameters
        ----------
        fastq : str
            Location of the FastQ file to randomise
        output : str
            [OPTIONAL] Location of the output FastQ file. If left blank then the
            randomised output is saved to the same location as `fastq`
        """
        fqr = fastqreader()
        fqr.openFastQ(fastq)

        record = fqr.next()

        fq_tmp = {}
        while fqr.eof() is False:
            r_id = record["id"].split(" ")
            fq_tmp[r_id[0]] = [record["start_posn"], record["end_posn"]]
            record = fqr.next()

        fqr.closeFastQ()

        fastq_suffix = "." + fastq[-1].split(".")[-1]

        with open(fastq, "r") as fq_file:
            with open(fastq.replace(fastq_suffix, ".random"), "w") as tmp_fq_file:
                for r_id in list(fq_tmp.keys()):
                    fq_file.seek(fq_tmp[r_id][0])
                    tmp_fq_file.write(fq_file.read(fq_tmp[r_id][1]-fq_tmp[r_id][0]))

        output_fastq = fastq
        if output is not None:
            output_fastq = output
        with open(fastq + ".random", "rb") as fq_random:
            with open(output_fastq, "wb") as f_out:
                f_out.write(fq_random.read())

        if output is not None:
            os.remove(fastq + ".random")
