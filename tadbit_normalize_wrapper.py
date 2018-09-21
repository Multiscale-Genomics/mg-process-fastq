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
import argparse
import sys
import json
import multiprocessing
import tarfile
from random import random
from string import ascii_letters as letters

from pysam import AlignmentFile  # pylint: disable=no-name-in-module

from basic_modules.workflow import Workflow
from basic_modules.metadata import Metadata
from utils import logger
from utils import remap

from tool.common import CommandLineParser
from tool.common import format_utils
from tool.tb_normalize import tbNormalizeTool


# ------------------------------------------------------------------------------

class tadbit_normalize(Workflow):  # pylint: disable=invalid-name, too-few-public-methods
    """
    Wrapper for the VRE form TADbit normalize.
    It normalizes a BAM file at a given resolution.
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
        tool_extra_config = json.load(file(os.path.dirname(
            os.path.abspath(__file__))+'/tadbit_wrappers_config.json'))
        if os.path.isdir(format_utils.convert_from_unicode(
                tool_extra_config["bin_path"])):
            os.environ["PATH"] += os.pathsep + format_utils.convert_from_unicode(
                tool_extra_config["bin_path"])

        if configuration is None:
            configuration = {}

        self.configuration.update(format_utils.convert_from_unicode(configuration))

        # Number of cores available
        num_cores = multiprocessing.cpu_count()
        self.configuration["ncpus"] = num_cores

        if "normalization" not in self.configuration:
            self.configuration["normalization"] = "Vanilla"

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
            "PROCESS NORMALIZE - FILES PASSED TO TOOLS: {0}".format(str(input_files["bamin"]))
        )

        bamin = format_utils.convert_from_unicode(input_files['bamin'])
        input_metadata = remap(self.configuration,
                               "normalization", "resolution", "min_perc", "workdir",
                               "max_perc", "ncpus")

        assembly = format_utils.convert_from_unicode(metadata['bamin'].meta_data['assembly'])
        if self.configuration["normalization"] == 'oneD':
            if 'refGenomes_folder' in input_files and os.path.isfile(
                    format_utils.convert_from_unicode(
                        input_files['refGenomes_folder'])+assembly+'/'+assembly+'.fa'):
                input_metadata["fasta"] = format_utils.convert_from_unicode(
                    input_files['refGenomes_folder'])+assembly+'/'+assembly+'.fa'
                input_metadata["mappability"] = format_utils.convert_from_unicode(
                    input_files['refGenomes_folder'])+assembly+'/MAPPABILITY/'+assembly+'.bedGraph'
                if 'rest_enzyme' in metadata['bamin'].meta_data:
                    input_metadata["rest_enzyme"] = format_utils.convert_from_unicode(
                        metadata['bamin'].meta_data['rest_enzyme'])

            if 'rest_enzyme' not in input_metadata or 'fasta' not in input_metadata \
                    or 'mappability' not in input_metadata:
                logger.fatal('Error: missing parameters for oneD normalization.\
                    Please check that the BAM input file has been generated with the VRE tool.')
                return {}, {}

        bamfile = AlignmentFile(bamin, 'rb')
        if len(bamfile.references) == 1:
            input_metadata["min_count"] = "10"
        bamfile.close()

        m_results_meta = {}

        tn_handler = tbNormalizeTool()
        tn_files, _ = tn_handler.run([bamin], input_metadata, [])

        m_results_files = {}
        try:
            m_results_files["hic_biases"] = self.configuration['project'] + "/" + \
                os.path.basename(tn_files[0])
            os.rename(tn_files[0], m_results_files["hic_biases"])
        except OSError:
            pass

        m_results_files["normalize_stats"] = self.configuration['project'] + "/normalize_stats.tar.gz"  # pylint: disable=line-too-long

        with tarfile.open(m_results_files["normalize_stats"], "w:gz") as tar:
            if len(tn_files) > 1:
                tar.add(tn_files[1], arcname=os.path.basename(tn_files[1]))
            if len(tn_files) > 2:
                tar.add(tn_files[2], arcname=os.path.basename(tn_files[2]))

        # List of files to get saved
        logger.info("TADBIT RESULTS: " + m_results_files["normalize_stats"])

        m_results_meta["hic_biases"] = Metadata(
            data_type="hic_biases",
            file_type="PICKLE",
            file_path=m_results_files["hic_biases"],
            sources=[bamin],
            meta_data={
                "description": "HiC biases for normalization",
                "visible": True,
                "assembly": format_utils.convert_from_unicode(
                    metadata['bamin'].meta_data['assembly']),
                "norm": self.configuration["normalization"]
            },
            taxon_id=metadata['bamin'].taxon_id)
        m_results_meta["normalize_stats"] = Metadata(
            data_type="tool_statistics",
            file_type="TAR",
            file_path=m_results_files["normalize_stats"],
            sources=[bamin],
            meta_data={
                "description": "TADbit normalize statistics",
                "visible": False
            })

        # cleaning
        clean_temps(self.configuration['workdir']+"/04_normalization")
        clean_temps(self.configuration['workdir'])

        return m_results_files, m_results_meta


# ------------------------------------------------------------------------------

def main(args):
    """
    Main function
    """
    from apps.jsonapp import JSONApp
    app = JSONApp()
    result = app.launch(tadbit_normalize,
                        args.config,
                        args.in_metadata,
                        args.out_metadata)

    return result


def clean_temps(working_path):
    """Cleans the workspace from temporal folder and scratch files"""
    for the_file in os.listdir(working_path):
        file_path = os.path.join(working_path, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            # elif os.path.isdir(file_path): shutil.rmtree(file_path)
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
