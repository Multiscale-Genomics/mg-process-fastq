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
import shlex
import subprocess
import sys

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task
    from utils.dummy_pycompss import compss_wait_on

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

# ------------------------------------------------------------------------------

class macs2(Tool):
    """
    Tool for peak calling for ChIP-seq data
    """

    def __init__(self, configuration=None):
        """
        Init function
        """
        logger.info("MACS2 Peak Caller")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @staticmethod
    @task(
        returns=int,
        name=IN,
        macs_params=IN,
        bam_file=FILE_IN,
        bam_file_bgd=FILE_IN,
        narrowpeak=FILE_OUT,
        summits_bed=FILE_OUT,
        broadpeak=FILE_OUT,
        gappedpeak=FILE_OUT,
        isModifier=False)
    def macs2_peak_calling(
            name, bam_file, bam_file_bgd, macs_params,
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

        bgd_command = '-c ' + bam_file_bgd

        command_param = [
            'macs2 callpeak', " ".join(macs_params), '-t', bam_file, '-n', name, bgd_command,
            ' --outdir ', output_dir
        ]
        command_line = ' '.join(command_param)

        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        if process.returncode is not 0:
            logger.fatal("MACS2 ERROR", process.returncode)
            return process.returncode

        logger.info('Process Results 1:', process)
        logger.info('LIST DIR 1:', os.listdir(output_dir))

        out_suffix = ['peaks.narrowPeak', 'peaks.broadPeak', 'peaks.gappedPeak', 'summits.bed']
        for f_suf in out_suffix:
            output_tmp = output_dir + '/' + name + '_out_' + f_suf
            logger.info(output_tmp, os.path.isfile(output_tmp))
            if os.path.isfile(output_tmp) is True and os.path.getsize(output_tmp) > 0:
                if f_suf == 'peaks.narrowPeak':
                    with open(narrowpeak, "wb") as f_out:
                        with open(output_tmp, "rb") as f_in:
                            f_out.write(f_in.read())
                elif f_suf == 'summits.bed':
                    with open(summits_bed, "wb") as f_out:
                        with open(output_tmp, "rb") as f_in:
                            f_out.write(f_in.read())
                elif f_suf == 'peaks.broadPeak':
                    with open(broadpeak, "wb") as f_out:
                        with open(output_tmp, "rb") as f_in:
                            f_out.write(f_in.read())
                elif f_suf == 'peaks.gappedPeak':
                    with open(gappedpeak, "wb") as f_out:
                        with open(output_tmp, "rb") as f_in:
                            f_out.write(f_in.read())
            else:
                if f_suf == 'peaks.narrowPeak':
                    with open(narrowpeak, "w") as f_out:
                        f_out.write("")
                elif f_suf == 'summits.bed':
                    with open(summits_bed, "w") as f_out:
                        f_out.write("")
                elif f_suf == 'peaks.broadPeak':
                    with open(broadpeak, "w") as f_out:
                        f_out.write("")
                elif f_suf == 'peaks.gappedPeak':
                    with open(gappedpeak, "w") as f_out:
                        f_out.write("")

        return 0

    @staticmethod
    @task(returns=int, name=IN, macs_params=IN, bam_file=FILE_IN,
          narrowpeak=FILE_OUT, summits_bed=FILE_OUT, broadpeak=FILE_OUT,
          gappedpeak=FILE_OUT, isModifier=False)
    def macs2_peak_calling_nobgd( # pylint: disable=too-many-arguments
            name, bam_file, macs_params,
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
        output_dir = "/".join(od_list[0:-1])

        command_line = "macs2 callpeak " + " ".join(macs_params) + " -t " + bam_file
        command_line = command_line + ' -n ' + name + '_out --outdir ' + output_dir

        logger.info("MACS2 - NAME: " + name)
        logger.info("Output Files: " + narrowpeak, summits_bed, broadpeak, gappedpeak)
        logger.info("MACS2 COMMAND LINE: " + command_line)

        if os.path.isfile(bam_file) is False or os.path.getsize(bam_file) == 0:
            logger.fatal("MISSING FILE: " + bam_file)

        try:
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()
        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, command_line))
            return msg.errno

        if process.returncode is not 0:
            logger.fatal("MACS2 ERROR: " + str(process.returncode))
            return process.returncode

        out_suffix = ['peaks.narrowPeak', 'peaks.broadPeak', 'peaks.gappedPeak', 'summits.bed']
        for f_suf in out_suffix:
            output_tmp = output_dir + '/' + name + '_out_' + f_suf
            logger.info(output_tmp, os.path.isfile(output_tmp))
            if os.path.isfile(output_tmp) is True and os.path.getsize(output_tmp) > 0:
                if f_suf == 'peaks.narrowPeak':
                    with open(narrowpeak, "wb") as f_out:
                        with open(output_tmp, "rb") as f_in:
                            f_out.write(f_in.read())
                elif f_suf == 'summits.bed':
                    with open(summits_bed, "wb") as f_out:
                        with open(output_tmp, "rb") as f_in:
                            f_out.write(f_in.read())
                elif f_suf == 'peaks.broadPeak':
                    with open(broadpeak, "wb") as f_out:
                        with open(output_tmp, "rb") as f_in:
                            f_out.write(f_in.read())
                elif f_suf == 'peaks.gappedPeak':
                    with open(gappedpeak, "wb") as f_out:
                        with open(output_tmp, "rb") as f_in:
                            f_out.write(f_in.read())
            else:
                if f_suf == 'peaks.narrowPeak':
                    with open(narrowpeak, "w") as f_out:
                        f_out.write("")
                elif f_suf == 'summits.bed':
                    with open(summits_bed, "w") as f_out:
                        f_out.write("")
                elif f_suf == 'peaks.broadPeak':
                    with open(broadpeak, "w") as f_out:
                        f_out.write("")
                elif f_suf == 'peaks.gappedPeak':
                    with open(gappedpeak, "w") as f_out:
                        f_out.write("")

        return 0

    @staticmethod
    def get_macs2_params(params):
        """
        Function to handle to extraction of commandline parameters and formatting
        them for use in the aligner for BWA ALN

        Parameters
        ----------
        params : dict

        Returns
        -------
        list
        """
        command_params = []

        command_parameters = {
            "macs_gsize_param" : ["--gsize", True],
            "macs_tsize_param" : ["--tsize", True],
            "macs_bw_param" : ["--bw", True],
            "macs_qvalue_param" : ["--qvalue", True],
            "macs_pvalue_param" : ["--pvalue", True],
            "macs_mfold_param" : ["--mfold", True],
            "macs_nolambda_param" : ["--nolambda", False],
            "macs_slocal_param" : ["--slocal", True],
            "macs_llocal_param" : ["--llocal", True],
            "macs_fix-bimodal_param" : ["--fix-bimodal", False],
            "macs_nomodel_param" : ["--nomodel", False],
            "macs_extsize_param" : ["--extsize", True],
            "macs_shift_param" : ["--shift", True],
            "macs_keep-dup_param" : ["--keep-dup", True],
            "macs_broad_param" : ["--broad", False],
            "macs_broad-cutoff_param" : ["--broad-cutoff", True],
            "macs_to-large_param" : ["--to-large", False],
            "macs_down-sample_param" : ["--down-sample", False],
            "macs_bdg_param" : ["--bdg", True],
            "macs_call-summits_param" : ["--call-summits", True],
        }

        for param in params:
            if param in command_parameters:
                if command_parameters[param][1]:
                    command_params = command_params + [command_parameters[param][0], params[param]]
                else:
                    if command_parameters[param][0]:
                        command_params.append(command_parameters[param][0])

        return command_params

    def run(self, input_files, input_metadata, output_files):
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
        root_name = input_files['input'].split("/")
        root_name[-1] = root_name[-1].replace('.bam', '')
        name = root_name[-1]

       # input and output share most metadata
        output_bed_types = {
            'narrow_peak': "bed4+1",
            'summits': "bed6+4",
            'broad_peak': "bed6+3",
            'gapped_peak': "bed12+3"
        }

        command_params = self.get_macs2_params(self.configuration)

        logger.info("MACS2 COMMAND PARAMS:" + ", ".join(command_params))

        # handle error
        if 'background' in input_files:
            results = self.macs2_peak_calling(
                name, str(input_files['input']), str(input_files['background']),
                command_params,
                str(output_files['narrow_peak']), str(output_files['summits']),
                str(output_files['broad_peak']), str(output_files['gapped_peak']))
        else:
            results = self.macs2_peak_calling_nobgd(
                name, str(input_files['input']), command_params,
                str(output_files['narrow_peak']), str(output_files['summits']),
                str(output_files['broad_peak']), str(output_files['gapped_peak']))
        results = compss_wait_on(results)

        if results > 0:
            logger.fatal("MACS2 failed with error: {}", results)
            return (
                {}, {}
            )

        output_files_created = {}
        output_metadata = {}
        for result_file in output_files:
            if (
                    os.path.isfile(output_files[result_file]) is True
                    and os.path.getsize(output_files[result_file]) > 0
                ):
                output_files_created[result_file] = output_files[result_file]

                output_metadata[result_file] = Metadata(
                    data_type="data_chip_seq",
                    file_type="BED",
                    file_path=output_files[result_file],
                    sources=[input_metadata['input'].file_path],
                    taxon_id=input_metadata["input"].taxon_id,
                    meta_data={
                        "assembly": input_metadata["input"].meta_data["assembly"],
                        "tool": "macs2",
                        "bed_type": output_bed_types[result_file]
                    }
                )

        logger.info('MACS2: GENERATED FILES:', output_files)

        return (output_files_created, output_metadata)

# ------------------------------------------------------------------------------
