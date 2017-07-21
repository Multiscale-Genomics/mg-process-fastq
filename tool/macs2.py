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

import os
import shlex
import subprocess
import sys

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    #from pycompss.api.api import compss_wait_on
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT, IN
    from dummy_pycompss import task

#from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

#from tool.common import common

# ------------------------------------------------------------------------------

class macs2(Tool):
    """
    Tool for peak calling for ChIP-seq data
    """

    def __init__(self, configuration=None):
        """
        Init function
        """
        print("MACS2 Peak Caller")
        Tool.__init__(self)

    @task(
        returns=bool,
        name=IN,
        bam_file=FILE_IN,
        bam_file_bgd=FILE_IN,
        narrowpeak=FILE_OUT,
        summit_bed=FILE_OUT,
        broadpeak=FILE_OUT,
        gappedpeak=FILE_OUT,
        isModifier=False)
    def macs2_peak_calling(
            self, name, bam_file, bam_file_bgd,
            narrowpeak, summits_bed, broadpeak, gappedpeak): # pylint: disable=unused-argument
        """
        Function to run MACS2 for peak calling on aligned sequence files and
        normalised against a provided background set of alignments.

        Parameters
        ----------
        name : str
            Name to be used to identify the files
        bam_file : str
            Location of the aligned FASTQ files as a bam file
        bam_file_bgd : str
            Location of the aligned FASTQ files as a bam file representing
            background values for the cell

        Returns
        -------
        narrowPeak : file
            BED6+4 file - ideal for transcription factor binding site
            identification
        summitPeak : file
            BED4+1 file - Contains the peak summit locations for everypeak
        broadPeak : file
            BED6+3 file - ideal for histone binding site identification
        gappedPeak : file
            BED12+3 file - Contains a merged set of the broad and narrow peak
            files

        Definitions defined for each of these files have come from the MACS2
        documentation described in the docs at https://github.com/taoliu/MACS
        """

        od_list = bam_file.split("/")
        output_dir = "/".join(od_list[0:-1])

        bgd_command = ''
        if bam_file_bgd is not None:
            bgd_command = '-c ' + bam_file_bgd

        nomodel = ''
        if name == 'macs2.Human.DRR000150.22.filtered':
            # This is for when running the test data
            nomodel = '--nomodel'

        command_param = [
            'macs2 callpeak', '-t', bam_file, 'n', name, bgd_command,
            ' --outdir ', output_dir, nomodel
        ]
        command_line = ' '.join(command_param)

        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        print('Process Results 1:', process)
        print('LIST DIR 1:', os.listdir(output_dir))

        # Might not be an issue with PyCOMPSs v2.1
        #out_peaks_narrow = bam_file + '_peaks.narrowPeak'
        #out_peaks_broad = bam_file + '_summits.bed'
        #out_peaks_gapped = bam_file + '_summits.bed'
        #out_summits = bam_file + '_summits.bed'

        #os.rename(out_peaks_narrow, narrowpeak)
        #os.rename(out_summits, summits_bed)
        #os.rename(out_peaks_broad, broadpeak)
        #os.rename(out_peaks_gapped, gappedpeak)

        return 0

    @task(
        returns=bool,
        name=IN,
        bam_file=FILE_IN,
        narrowPeak=FILE_OUT,
        summit_bed=FILE_OUT,
        broadPeak=FILE_OUT,
        gappedPeak=FILE_OUT,
        isModifier=False)
    def macs2_peak_calling_nobgd( # pylint: disable=too-many-arguments
            self, name, bam_file,
            narrowpeak, summits_bed, broadpeak, gappedpeak): # pylint: disable=unused-argument
        """
        Function to run MACS2 for peak calling on aligned sequence files without
        a background dataset for normalisation.

        Parameters
        ----------
        name : str
            Name to be used to identify the files
        bam_file : str
            Location of the aligned FASTQ files as a bam file

        Returns
        -------
        narrowPeak : file
            BED6+4 file - ideal for transcription factor binding site
            identification
        summitPeak : file
            BED4+1 file - Contains the peak summit locations for everypeak
        broadPeak : file
            BED6+3 file - ideal for histone binding site identification
        gappedPeak : file
            BED12+3 file - Contains a merged set of the broad and narrow peak
            files

        Definitions defined for each of these files have come from the MACS2
        documentation described in the docs at https://github.com/taoliu/MACS
        """
        od_list = bam_file.split("/")
        #od_list = narrowpeak.split("/")
        output_dir = "/".join(od_list[0:-1])

        command_line = 'macs2 callpeak -t ' + bam_file + ' -n ' + name + ' --outdir ' + output_dir

        if name == 'macs2.Human.DRR000150.22.filtered':
            # This is for when running the test data
            command_line = command_line + ' --nomodel'

        print('Output Files:', narrowpeak, summits_bed, broadpeak, gappedpeak)
        print(command_line)

        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        print('Process Results 2:', process.returncode)
        print('LIST DIR 2:', os.listdir(output_dir))
        

        # Might not be an issue with PyCOMPSs v2.1
        # out_peaks_narrow = bam_file + '_peaks.narrowPeak'
        # out_peaks_broad = bam_file + '_peaks.broadPeak'
        # out_peaks_gapped = bam_file + '_peaks.gappedPeak'
        # out_summits = bam_file + '_summits.bed'

        # os.rename(out_peaks_narrow, narrowpeak)
        # os.rename(out_summits, summits_bed)
        # os.rename(out_peaks_broad, broadpeak)
        # os.rename(out_peaks_gapped, gappedpeak)

        if process.returncode is 1:
            return False
        return True


    def run(self, input_files, output_files, metadata=None):
        """
        The main function to run MACS 2 for peak calling over a given BAM file
        and matching background BAM file.

        Parameters
        ----------
        input_files : list
            List of input bam file locations where 0 is the bam data file and 1
            is the matching background bam file
        metadata : dict


        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects

        """

        bam_file = input_files[0]

        print("Check bam file exists:", os.path.isfile(bam_file), os.path.getsize(bam_file))

        bam_file_bgd = None
        if len(input_files) == 2 and input_files[1] is not None:
            bam_file_bgd = input_files[1]

        root_name = bam_file.split("/")
        root_name[-1] = root_name[-1].replace('.bam', '')

        name = root_name[-1]
        #name = '/'.join(root_name)

        out_peaks_narrow = '/'.join(root_name) + '_peaks.narrowPeak'
        out_peaks_broad = '/'.join(root_name) + '_peaks.broadPeak'
        out_peaks_gapped = '/'.join(root_name) + '_peaks.gappedPeak'
        out_summits = '/'.join(root_name) + '_summits.bed'

        print("MACS2 output files:", name, out_peaks_narrow, out_peaks_broad, out_peaks_gapped, out_summits)

        # input and output share most metadata
        output_metadata = {}
        output_metadata["bed_types"] = ["bed4+1", "bed6+4", "bed6+3", "bed12+3"]

        # handle error
        if bam_file_bgd is None:
            results = self.macs2_peak_calling_nobgd(
                name, bam_file,
                out_peaks_narrow, out_summits, out_peaks_broad, out_peaks_gapped)
        else:
            results = self.macs2_peak_calling(
                name, bam_file, bam_file_bgd,
                out_peaks_narrow, out_summits, out_peaks_broad, out_peaks_gapped)
        #results = compss_wait_on(results)

        if results is False:
            return (
                [], []
            )

        print('Results:', results)

        return (
            [out_peaks_narrow, out_summits, out_peaks_broad, out_peaks_gapped],
            output_metadata
        )

# ------------------------------------------------------------------------------
