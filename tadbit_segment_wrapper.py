#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

import os.path
import shutil
import argparse
import sys
import json
import multiprocessing
import tarfile

from random import random
from string import ascii_letters as letters

from basic_modules.workflow import Workflow
from basic_modules.metadata import Metadata
from utils import logger
from utils import remap

from mg_process_fastq.tool.common import CommandLineParser
from mg_process_fastq.tool.common import format_utils
from mg_process_fastq.tool.tb_segment import tbSegmentTool


# ------------------------------------------------------------------------------

class tadbit_segment(Workflow):  # pylint: disable=invalid-name, too-few-public-methods
    """
    Wrapper for the VRE form TADbit segment.
    It detects TADs and compartments from a BAM file.
    """
    configuration = {}

    def __init__(self, configuration=None):
        """
        Initialise the tool with its configuration.


        Parameters
        ----------
        configuration : dict
            a dictionary containing parameters that define how the operation
            should be carried out, which are specific to each Tool.
        """

        tool_extra_config = json.load(open(os.path.dirname(
            os.path.abspath(__file__))+'/tadbit_wrappers_config.json'))
        os.environ["PATH"] += os.pathsep + format_utils.convert_from_unicode(
            tool_extra_config["bin_path"])

        if configuration is None:
            configuration = {}

        self.configuration.update(format_utils.convert_from_unicode(configuration))

        # Number of cores available
        num_cores = multiprocessing.cpu_count()
        self.configuration["ncpus"] = num_cores

        tmp_name = ''.join([letters[int(random()*52)]for _ in range(5)])
        if 'execution' in self.configuration:
            self.configuration['project'] = self.configuration['execution']
        self.configuration['workdir'] = self.configuration['project']+'/_tmp_tadbit_'+tmp_name
        if not os.path.exists(self.configuration['workdir']):
            os.makedirs(self.configuration['workdir'])
            open(self.configuration['workdir']+'/trace.db', 'a').close()

        self.configuration.update(
            {(key.split(':'))[-1]: val for key, val in self.configuration.items()}
        )

    def run(self, input_files, metadata, output_files):
        """
        Parameters
        ----------
        files_ids : list
            List of file locations
        metadata : list
            Required meta data
        output_files : list
            List of output file locations

        Returns
        -------
        outputfiles : list
            List of locations for the output files
        """

        logger.info(
            "PROCESS SEGMENT - FILES PASSED TO TOOLS: {0}".format(str(input_files["bamin"]))
        )
        m_results_files = {}
        m_results_meta = {}

        input_metadata = remap(self.configuration,
                               "resolution", "callers", "workdir", "ncpus")
        assembly = format_utils.convert_from_unicode(
            metadata['bamin'].meta_data['assembly'])
        if 'refGenomes_folder' in input_files and os.path.isfile(format_utils.convert_from_unicode(
                input_files['refGenomes_folder'])+assembly+'/'+assembly+'.fa'):
            input_metadata["fasta"] = format_utils.convert_from_unicode(
                input_files['refGenomes_folder'])+assembly+'/'+assembly+'.fa'
        if "chromosome_names" in self.configuration:
            input_metadata["chromosomes"] = self.configuration["chromosome_names"]

        in_files = [format_utils.convert_from_unicode(input_files['bamin'])]
        if 'biases' in input_files:
            in_files.append(format_utils.convert_from_unicode(input_files['biases']))

        ts_handler = tbSegmentTool()
        ts_files, _ = ts_handler.run(in_files, input_metadata, [])

        m_results_files["tads_compartments"] = self.configuration['project'] + \
            "/tads_compartments.tar.gz"

        tar = tarfile.open(m_results_files["tads_compartments"], "w:gz")
        if '1' in self.configuration['callers'] and '2' in self.configuration['callers']:
            tar.add(ts_files[0], arcname='tads')
            tar.add(ts_files[1], arcname='compartments')
        elif '1' in self.configuration['callers']:
            tar.add(ts_files[0], arcname='tads')
        elif '2' in self.configuration['callers']:
            tar.add(ts_files[0], arcname='compartments')

        tar.close()

        m_results_meta["tads_compartments"] = Metadata(
            data_type="tool_statistics",
            file_type="TAR",
            file_path=m_results_files["tads_compartments"],
            sources=in_files,
            meta_data={
                "description": "TADbit HiC tads and compartments statistics",
                "visible": True
            })
        # List of files to get saved
        logger.info("TADBIT RESULTS:" + m_results_files["tads_compartments"])
        # clean_temps(os.path.dirname(ts_files[0]))
        # clean_temps(os.path.join(self.configuration['workdir'],"06_segmentation"))
        # cleaning
        clean_temps(self.configuration['workdir'])

        return m_results_files, m_results_meta


# ------------------------------------------------------------------------------

def main(args):
    """
    Main function
    """
    from apps.jsonapp import JSONApp
    app = JSONApp()
    result = app.launch(tadbit_segment,
                        args.config,
                        args.in_metadata,
                        args.out_metadata)

    return result


# ------------------------------------------------------------------------------

def clean_temps(working_path):
    """Cleans the workspace from temporal folder and scratch files"""
    for the_file in os.listdir(working_path):
        file_path = os.path.join(working_path, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except OSError:
            pass
    try:
        os.rmdir(working_path)
    except OSError:
        pass
    logger.info('[CLEANING] Finished')


# ------------------------------------------------------------------------------

if __name__ == "__main__":
    sys._run_from_cmdl = True  # pylint: disable=protected-access

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="TADbit map")
    # Config file
    PARSER.add_argument("--config", help="Configuration JSON file",
                        type=CommandLineParser.valid_file, metavar="config", required=True)

    # Metadata
    PARSER.add_argument("--in_metadata", help="Project metadata",
                        metavar="in_metadata", required=True)
    # Output metadata
    PARSER.add_argument("--out_metadata", help="Output metadata",
                        metavar="output_metadata", required=True)
    # Log file
    PARSER.add_argument("--log_file", help="Log file",
                        metavar="log_file", required=True)

    IN_ARGS = PARSER.parse_args()

    RESULTS = main(IN_ARGS)
