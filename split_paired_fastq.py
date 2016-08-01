#!/usr/bin/env python

import argparse, os.path, sys

class fastqreader:
    
    def __init__(self):
        """
        Initialise the module and 
        """
        self.fastq1 = ''
        self.fastq2 = ''
        
        self.f1_file = None
        self.f2_file = None
        
        self.f1_eof = False
        self.f2_eof = False
        
        self.f1_output_file = None
        self.f2_output_file = None
        
        self.output_tag = ''
        self.output_file_count = 0
    
    def openPairedFastQ(self, file1, file2):
        """
        Create file handles for reading the FastQ files
        """
        self.fastq1 = file1
        self.fastq2 = file2
        
        self.f1_file = open(self.fastq1, "r")
        self.f2_file = open(self.fastq2, "r")
        
        self.f1_eof = False
        self.f2_eof = False
    
    def closePairedFastQ(self):
        """
        Close file handles for the FastQ files.
        """
        self.f1_file.close()
        self.f2_file.close()
    
    def eof(self, side = 1):
        """
        Indicate if the end of the file has been reached
        """
        if side == 1:
            return self.f1_eof
        elif side == 2:
            return self.f2_eof
        else:
            return "ERROR"
    
    def next(self, side = 1):
        """
        Get the next read element for the specific FastQ file pair
        """
        read_id = ''
        read_seq = ''
        read_addition = ''
        read_score = ''
        
        if side == 1:
            read_id = self.f1_file.readline()
            read_seq = self.f1_file.readline()
            read_addition = self.f1_file.readline()
            read_score = self.f1_file.readline()
        elif side == 2:
            read_id = self.f2_file.readline()
            read_seq = self.f2_file.readline()
            read_addition = self.f2_file.readline()
            read_score = self.f2_file.readline()
        else:
            return 'ERROR'
        
        if read_id == '':
            if side == 1:
                self.f1_eof = True
                return False
            if side == 2:
                self.f2_eof = True
                return False
        return {'id': read_id.rstrip(), 'seq': read_seq.rstrip(), 'add': read_addition.rstrip(), 'score': read_score.rstrip()}
    
    def createOutputFiles(self, tag = ''):
        if tag != '' and self.output_tag != tag:
            self.output_tag = tag
        
        """
        Create and open the file handles for the output files
        """
        f1 = self.fastq1.split("/")
        f1[-1] = f1[-1].replace(".fastq", "." + str(self.output_tag) + "_" + str(self.output_file_count) + ".fastq")
        f1.insert(-1, "tmp")
        
        f2 = self.fastq2.split("/")
        f2[-1] = f2[-1].replace(".fastq", "." + str(self.output_tag) + "_" + str(self.output_file_count) + ".fastq")
        f2.insert(-1, "tmp")
        
        self.f1_output_file = open("/".join(f1), "w")
        self.f2_output_file = open("/".join(f2), "w")
    
    def writeOutput(self, read, side = 1):
        """
        Writer to print the extracted lines
        """
        line = read["id"] + "\n" + read["seq"] + "\n" + read["add"] + "\n" + read["score"] + "\n"
        if side == 1:
            self.f1_output_file.write(line)
        elif side == 2:
            self.f2_output_file.write(line)
        else:
            return False
        return True
    
    def closeOutputFiles(self):
        """
        Close the output file handles
        """
        self.f1_output_file.close()
        self.f2_output_file.close()
    
    def incrementOutputFiles(self):
        """
        Increment the counter and create new files for splitting the original
        FastQ paired end files.
        """
        self.closeOutputFiles()
        
        self.output_file_count+=1
        
        self.createOutputFiles(self.output_tag)
    
        
# Set up the command line parameters
parser = argparse.ArgumentParser(description="Load adjacency list into HDF5 file")
parser.add_argument("--input_1", help="File 1")
parser.add_argument("--input_2", help="File 2")
parser.add_argument("--output_tag", help="Inserted before the file descriptor and after the file name: e.g. 'matching' would convert file_id-1.fastq to file_id-1.matching.fastq")

args = parser.parse_args()
file1 = args.input_1
file2 = args.input_2
tag = args.output_tag

fqr = fastqreader()
fqr.openPairedFastQ(file1, file2)
fqr.createOutputFiles(tag)

r1 = fqr.next(1)
r2 = fqr.next(2)

count_r1 = 1
count_r2 = 2
count_r3 = 0

while fqr.eof(1) == False and fqr.eof(2) == False:
    r1_id = r1["id"].split(" ")
    r2_id = r2["id"].split(" ")
    
    if r1_id[0] == r2_id[0]:
        fqr.writeOutput(r1, 1)
        fqr.writeOutput(r2, 2)
        
        r1 = fqr.next(1)
        r2 = fqr.next(2)
        
        count_r1 += 1
        count_r2 += 1
        count_r3 += 1
    elif r1_id[0] < r2_id[0]:
        r1 = fqr.next(1)
        count_r1 += 1
    else:
        r2 = fqr.next(2)
        count_r2 += 1
    
    if count_r3 % 1000000 == 0:
        print count_r1, r1_id[0], count_r2, r2_id[0], count_r3
        fqr.incrementOutputFiles()

fqr.closePairedFastQ()
fqr.closeOutputFiles()

print count_r1, r1_id[0], count_r2, r2_id[0], count_r3
