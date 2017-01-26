#!/usr/bin/env python

"""
Copyright 2017 EMBL-European Bioinformatics Institute

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
