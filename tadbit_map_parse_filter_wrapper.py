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
import tarfile
import multiprocessing
import json
from random import random
from string import ascii_letters as letters

from basic_modules.workflow import Workflow
from basic_modules.metadata import Metadata
from utils import logger
from utils import remap

from mg_process_fastq.tool.common import CommandLineParser
from mg_process_fastq.tool.common import format_utils
from mg_process_fastq.tool.tb_full_mapping import tbFullMappingTool
from mg_process_fastq.tool.tb_parse_mapping import tbParseMappingTool
from mg_process_fastq.tool.tb_filter import tbFilterTool

if '/opt/COMPSs/Bindings/python' in sys.path:
    sys.path.pop(sys.path.index('/opt/COMPSs/Bindings/python'))


# ------------------------------------------------------------------------------

class tadbit_map_parse_filter(Workflow):  # pylint: disable=invalid-name,too-few-public-methods
    """
    Wrapper for the VRE form TADbit map, parse and filter.
    It combines different tools to map, merge and filter
    two fastq files corresponding to read1 and read2.
    """
    configuration = {}

    def __init__(self, configuration=None):  # pylint: disable=too-many-branches
        """
        Initialise the tool with its configuration.


        Parameters
        ----------
        configuration : dict
            a dictionary containing parameters that define how the operation
            should be carried out, which are specific to each Tool.
        """
        tool_extra_config = json.load(open(os.path.dirname(os.path.abspath(__file__))
                                           + '/tadbit_wrappers_config.json'))
        if os.path.isdir(format_utils.convert_from_unicode(tool_extra_config["bin_path"])):
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

        if 'mapping:refGenome' in self.configuration:
            self.configuration['mapping_refGenome'] = self.configuration['mapping:refGenome']
        if 'parsing:refGenome' in self.configuration:
            self.configuration['parsing_refGenome'] = self.configuration['parsing:refGenome']
        self.configuration.update(
            {(key.split(':'))[-1]: val for key, val in self.configuration.items()}
        )
        if 'filters' in self.configuration:
            self.configuration['filters'] = [int(f) for f in self.configuration['filters']]
        if 'windows' in self.configuration:
            if self.configuration['windows']:
                windows1 = self.configuration['windows'].split(" ")
                self.configuration['windows'] = [tuple(map(int, x.split(':'))) for x in windows1]
            else:
                self.configuration['windows'] = ''
        else:
            self.configuration['windows'] = ''
        if 'chromosomes' not in self.configuration:
            self.configuration["chromosomes"] = ''
        if 'rest_enzyme' not in self.configuration:
            self.configuration["rest_enzyme"] = ''

    def run(self, input_files, metadata, output_files):  # pylint: disable=too-many-locals,too-many-statements
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
            List of locations for the output bam files
        """
        logger.info(
            "PROCESS MAP - FILES PASSED TO TOOLS:",
            remap(input_files, "read1", "read2")
        )

        m_results_meta = {}
        m_results_files = {}

        assembly = "UNK"
        if 'parsing:refGenome' in input_files:
            genome_fa = format_utils.convert_from_unicode(input_files['parsing:refGenome'])
        elif 'parsing_refGenome' in self.configuration:
            genome_fa = self.configuration['public_dir'] + \
                format_utils.convert_from_unicode(self.configuration['parsing_refGenome'])

        if 'mapping:refGenome' in input_files:
            genome_gem = format_utils.convert_from_unicode(input_files['mapping:refGenome'])
            assembly = format_utils.convert_from_unicode(
                metadata['mapping:refGenome'].meta_data['assembly'])
        elif 'mapping_refGenome' in self.configuration:
            genome_gem = self.configuration['public_dir'] + \
                format_utils.convert_from_unicode(self.configuration['mapping_refGenome'])
            assembly = os.path.basename(genome_gem).split('.')[0]

        fastq_file_1 = format_utils.convert_from_unicode(input_files['read1'])
        fastq_file_2 = format_utils.convert_from_unicode(input_files['read2'])
        input_metadata = remap(self.configuration, "ncpus", "iterative_mapping",
                               "workdir", "windows", enzyme_name="rest_enzyme")
        input_metadata['quality_plot'] = True
        summary_file = input_metadata["workdir"]+'/'+'summary.txt'

        logger.info("MAPPING")
        logger.info("Read 1 of 2")
        tfm1 = tbFullMappingTool()
        tfm1_files, tfm1_meta = tfm1.run([genome_gem, fastq_file_1], input_metadata, [])
        with open(summary_file, 'w') as outfile:
            outfile.write("Mapping read 1\n--------------\n")
            with open(tfm1_files[-2]) as infile:
                outfile.write(infile.read())

        logger.info("Read 2 of 2")
        tfm2 = tbFullMappingTool()
        tfm2_files, tfm2_meta = tfm2.run([genome_gem, fastq_file_2], input_metadata, [])
        with open(summary_file, 'a') as outfile:
            outfile.write("\n\nMapping read 2\n--------------\n")
            with open(tfm2_files[-2]) as infile:
                outfile.write(infile.read())

        logger.info("PARSING")
        tpm = tbParseMappingTool()
        files = [genome_fa] + tfm1_files[:-2] + tfm2_files[:-2]

        input_metadata = remap(self.configuration, "ncpus", "chromosomes",
                               "workdir", enzyme_name="rest_enzyme")
        input_metadata['mapping'] = [tfm1_meta['func'], tfm2_meta['func']]
        input_metadata['expt_name'] = 'vre'

        logger.info("TB MAPPED FILES:", files)
        logger.info("TB PARSE METADATA:", input_metadata)
        tpm_files, tpm_meta = tpm.run(files, input_metadata, [])

        if 'error' in tpm_meta:
            m_results_meta["paired_reads"] = Metadata(
                data_type="hic_sequences",
                file_type="BAM",
                file_path=None,
                sources=[metadata['read1'].file_path, metadata['read2'].file_path],
                meta_data={
                    "tool": "tadbit",
                    "description": "Paired end reads",
                    "visible": True
                })
            m_results_meta["paired_reads"].error = True
            m_results_meta["paired_reads"].exception = tpm_meta['error']

        for file_path in tfm1_files[:-2]:
            if os.path.isfile(file_path):
                os.unlink(file_path)

        for file_path in tfm2_files[:-2]:
            if os.path.isfile(file_path):
                os.unlink(file_path)

        logger.info("TB PARSED FILES:", tpm_files)

        input_metadata = remap(self.configuration, "ncpus", "chromosomes", "workdir", 'filters')
        if 'min_dist_RE' in self.configuration:
            input_metadata['min_dist_RE'] = self.configuration['min_dist_RE']
        if 'min_fragment_size' in self.configuration:
            input_metadata['min_fragment_size'] = self.configuration['min_fragment_size']
        if 'max_fragment_size' in self.configuration:
            input_metadata['max_fragment_size'] = self.configuration['max_fragment_size']
        input_metadata['expt_name'] = 'vre'
        input_metadata['outbam'] = 'paired_reads'
        input_metadata['root_dir'] = self.configuration['project']
        input_metadata['custom_filter'] = True
        input_metadata['histogram'] = True

        logger.info("FILTERING")
        tbf = tbFilterTool()
        tf_files, _ = tbf.run(tpm_files, input_metadata, [])
        with open(summary_file, 'a') as outfile:
            outfile.write("\n\nFiltering\n--------------\n")
            with open(tf_files[-2]) as infile:
                outfile.write(infile.read())

        logger.info("TB FILTER FILES:", tf_files[0])

        m_results_files["paired_reads"] = tf_files[0]+'.bam'
        m_results_files["paired_reads_bai"] = m_results_files["paired_reads"]+'.bai'
        m_results_files["map_parse_filter_stats"] = self.configuration['project'] + \
            "/map_parse_filter_stats.tar.gz"

        with tarfile.open(m_results_files["map_parse_filter_stats"], "w:gz") as tar:
            tar.add(summary_file, arcname=os.path.basename(summary_file))
            tar.add(tfm1_files[-1], arcname=os.path.basename(tfm1_files[-1]))
            tar.add(tfm2_files[-1], arcname=os.path.basename(tfm2_files[-1]))
            tar.add(tf_files[-1], arcname=os.path.basename(tf_files[-1]))

        # List of files to get saved
        logger.info("TADBIT RESULTS:", m_results_files)

        m_results_meta["paired_reads"] = Metadata(
            data_type="hic_sequences",
            file_type="BAM",
            file_path=m_results_files["paired_reads"],
            sources=[metadata['read1'].file_path, metadata['read2'].file_path],
            meta_data={
                "description": "Paired end reads",
                "visible": True,
                "assembly": assembly,
                "paired": "paired",
                "sorted": "sorted",
                "associated_files": [
                    m_results_files["paired_reads"]+'.bai'
                ],
                "func": tfm1_meta['func'],
                "rest_enzyme": self.configuration["rest_enzyme"]
            },
            taxon_id=metadata['read1'].taxon_id)

        m_results_meta["paired_reads_bai"] = Metadata(
            data_type="hic_sequences",
            file_type="BAI",
            file_path=m_results_files["paired_reads"]+'.bai',
            sources=[metadata['read1'].file_path, metadata['read2'].file_path],
            meta_data={
                "description": "Paired end reads index",
                "visible": False,
                "assembly": assembly,
                "associated_master": m_results_files["paired_reads"]
            },
            taxon_id=metadata['read1'].taxon_id)

        m_results_meta["map_parse_filter_stats"] = Metadata(
            data_type="tool_statistics",
            file_type="TAR",
            file_path=m_results_files["map_parse_filter_stats"],
            sources=[metadata['read1'].file_path, metadata['read2'].file_path],
            meta_data={
                "description": "TADbit mapping, parsing and filtering statistics",
                "visible": False
            })

        clean_temps(self.configuration['workdir'])

        return m_results_files, m_results_meta


# ------------------------------------------------------------------------------

def main(args):
    """
    Main function
    """
    from apps.jsonapp import JSONApp
    app = JSONApp()
    result = app.launch(tadbit_map_parse_filter,
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
    PARSER.add_argument("--log_file", help="Log file", metavar="log_file", required=True)

    IN_ARGS = PARSER.parse_args()

    RESULTS = main(IN_ARGS)
