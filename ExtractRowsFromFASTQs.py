#!/usr/bin/python

"""
Script to extract FASTQ entry pairs from paired end FASTQ files given a list of
FASTQ IDs.
"""

import argparse

from fastqreader import fastqreader

# Extract the rows that matched to chr21
# grep chr21 hic_test/map/SRR1658573_1_full_1-100.map | tr "\t" "~" | cut -d"~" -f1 -f5 | tr ":" "\t" | awk '(NR==1) || (($4>15000000) && ($4<20000000))' | tr "\t" "~" | cut -d "~" -f1 > SRR1658573_1_chr21_1-100.row
# grep chr21 hic_test/map/SRR1658573_2_full_1-100.map | tr "\t" "~" | cut -d"~" -f1 -f5 | tr ":" "\t" | awk '(NR==1) || (($4>15000000) && ($4<20000000))' | tr "\t" "~" | cut -d "~" -f1 > SRR1658573_2_chr21_1-100.row

# grep -Fx -f SRR1658573_1_chr21_1-100.row SRR1658573_2_chr21_1-100.row > SRR1658573_chr21.row

def paired_selector(in_file1, in_file2, rows, tag='tmp'):
    """
    Function to divide the FastQ files into separate sub files of 1000000
    sequences so that the aligner can run in parallel.

    Parameters
    ----------
    in_file1 : str
        Location of first paired end FASTQ file
    in_file2 : str
        Location of second paired end FASTQ file
    tag : str
        DEFAULT = tmp
        Tag used to identify the files. Useful if this is getting run
        manually on a single machine multiple times to prevent collisions of
        file names


    Returns
    -------
    Returns: Returns a list of lists of the files that have been generated.
             Each sub list containing the two paired end files for that
             subset.
    paired_files : list
        List of lists of pair end files. Each sub list containing the two
        paired end files for that subset.
    """

    fqr = fastqreader()
    fqr.openFastQ(in_file1, in_file2)
    fqr.createOutputFiles(tag)

    record1 = fqr.next(1)
    record2 = fqr.next(2)

    count_r1 = 0
    count_r2 = 0
    count_r3 = 0
    count_r4 = 0

    file_loc_1 = fqr.fastq1.split("/")
    file_loc_1[-1] = file_loc_1[-1].replace(
        ".fastq",
        "." + str(fqr.output_tag) + "_" + str(fqr.output_file_count) + ".fastq")
    file_loc_1.insert(-1, tag)

    file_loc_2 = fqr.fastq2.split("/")
    file_loc_2[-1] = file_loc_2[-1].replace(
        ".fastq",
        "." + str(fqr.output_tag) + "_" + str(fqr.output_file_count) + ".fastq")
    file_loc_2.insert(-1, tag)
    files_out = [["/".join(file_loc_1), "/".join(file_loc_2)]]

    while fqr.eof(1) is False and fqr.eof(2) is False:
        r1_id = record1["id"].split(" ")
        r2_id = record2["id"].split(" ")

        if r1_id[0] == r2_id[0]:
            if r1_id[0][1:] in rows:
                fqr.writeOutput(record1, 1)
                fqr.writeOutput(record2, 2)
                count_r4 += 1

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

        if count_r3 % 1000 == 0:
            print(count_r3, count_r4)

    fqr.closeFastQ()
    fqr.closeOutputFiles()

    return files_out


def get_if_list(row_file):
    """
    Get the IDs from the specified file
    """
    fid = open(row_file, 'r')

    id_list = fid.readlines()
    id_list = [i.rstrip() for i in id_list]

    from sets import Set
    id_set = Set(id_list)

    fid.close()

    return list(id_set)


if __name__ == "__main__":

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="Load adjacency list into HDF5 file")
    PARSER.add_argument("--input_1", help="File 1")
    PARSER.add_argument("--input_2", help="File 2")
    PARSER.add_argument("--rows", help="Row File")
    PARSER.add_argument("--output_tag", help="Inserted before the file descriptor and after the file name: e.g. 'matching' would convert file_id-1.fastq to file_id-1.matching.fastq")

    ARGS = PARSER.parse_args()
    FILE_01 = ARGS.input_1
    FILE_02 = ARGS.input_2
    ROW_FILE = ARGS.rows
    TAG = ARGS.output_tag

    VALID_IDS = get_if_list(ROW_FILE)

    paired_selector(FILE_01, FILE_02, VALID_IDS, TAG)
