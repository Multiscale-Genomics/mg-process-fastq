#!/usr/bin/env python
#!/home/pmes/miniconda2/bin/python

# -*- coding: utf-8 -*-

from __future__ import print_function

import os.path
import shutil
import argparse
import sys
import json
import multiprocessing
import collections
import tarfile

from random import random
from string import ascii_letters as letters

# Required for ReadTheDocs
from functools import wraps # pylint: disable=unused-import

from basic_modules.workflow import Workflow
from basic_modules.metadata import Metadata

from tool.tb_segment import tbSegmentTool

class CommandLineParser(object):
    """Parses command line"""
    @staticmethod
    def valid_file(file_name):
        if not os.path.exists(file_name):
            raise argparse.ArgumentTypeError("The file does not exist")
        return file_name

    @staticmethod
    def valid_integer_number(ivalue):
        try:
            ivalue = int(ivalue)
        except:
            raise argparse.ArgumentTypeError("%s is an invalid value" % ivalue)
        if ivalue <= 0:
            raise argparse.ArgumentTypeError("%s is an invalid value" % ivalue)
        return ivalue
# ------------------------------------------------------------------------------
class tadbit_segment(Workflow):
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

        tool_extra_config = json.load(file(os.path.dirname(os.path.abspath(__file__))+'/tadbit_wrappers_config.json'))
        os.environ["PATH"] += os.pathsep + convert_from_unicode(tool_extra_config["bin_path"])

        if configuration is None:
            configuration = {}

        self.configuration.update(convert_from_unicode(configuration))

        #self.configuration['public_dir'] = '/orozco/services/MuG/MuG_public/refGenomes/'
        #self.configuration['public_dir'] = '/scratch/genomes/'

        # Number of cores available
        num_cores = multiprocessing.cpu_count()
        self.configuration["ncpus"] = num_cores

        tmp_name = ''.join([letters[int(random()*52)]for _ in xrange(5)])
        self.configuration['workdir'] = self.configuration['project']+'/_tmp_tadbit_'+tmp_name
        if not os.path.exists(self.configuration['workdir']):
            os.makedirs(self.configuration['workdir'])

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

        print(
            "PROCESS SEGMENT - FILES PASSED TO TOOLS:",
            remap(input_files, "bamin")
        )
        m_results_files = {}
        m_results_meta = {}
        #hic_data = load_hic_data_from_reads('/home/dcastillo/workspace/vre/mg-process-fastq-tadbit/tests/data/raw_None:0-13381_10kb.abc', resolution=10000)
        #exp = Experiment("vre", resolution=10000, hic_data=hic_data)

        input_metadata = remap(self.configuration, "resolution", "callers", "workdir", "ncpus")
        assembly = convert_from_unicode(metadata['bamin'].meta_data['assembly'])
        if 'refGenomes_folder' in input_files and os.path.isfile(convert_from_unicode(input_files['refGenomes_folder'])+assembly+'/'+assembly+'.fa'):
            input_metadata["fasta"] = convert_from_unicode(input_files['refGenomes_folder'])+assembly+'/'+assembly+'.fa'
        if "chromosome_names" in self.configuration:
            input_metadata["chromosomes"] = self.configuration["chromosome_names"]

        in_files = [convert_from_unicode(input_files['bamin'])]
        if 'biases' in input_files:
            in_files.append(convert_from_unicode(input_files['biases']))

        #hic_data = HiC_data((), len(bins_dict), sections, bins_dict, resolution=int(input_metadata['resolution']))
        ts = tbSegmentTool()
        ts_files, ts_meta = ts.run(in_files, [], input_metadata)

        m_results_files["tads_compartments"] = self.configuration['project']+"/tads_compartments.tar.gz"

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
        print("TADBIT RESULTS:", m_results_files)
        #clean_temps(os.path.dirname(ts_files[0]))
        #clean_temps(os.path.join(self.configuration['workdir'],"06_segmentation"))
        #cleaning
        clean_temps(self.configuration['workdir'])

        return m_results_files, m_results_meta

# ------------------------------------------------------------------------------

def remap(indict, *args, **kwargs):
    """
    Re-map keys of indict using information from arguments.
    Non-keyword arguments are keys of input dictionary that are passed
    unchanged to the output. Keyword arguments must be in the form
    old="new"
    and act as a translation table for new key names.
    """
    outdict = {role: indict[role] for role in args}
    outdict.update(
        {new: indict[old] for old, new in kwargs.items()}
    )
    return outdict

# ------------------------------------------------------------------------------

def convert_from_unicode(data):
    if isinstance(data, basestring):
        return str(data)
    if isinstance(data, collections.Mapping):
        return dict(map(convert_from_unicode, data.iteritems()))
    if isinstance(data, collections.Iterable):
        return type(data)(map(convert_from_unicode, data))
    return data

# ------------------------------------------------------------------------------

def main(args):

    from apps.jsonapp import JSONApp
    app = JSONApp()
    result = app.launch(tadbit_segment,
                        args.config,
                        args.in_metadata,
                        args.out_metadata)

    return result

#===============================================================================

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
    print('[CLEANING] Finished')

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    sys._run_from_cmdl = True # pylint: disable=protected-access

    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="TADbit map")
    # Config file
    parser.add_argument("--config", help="Configuration JSON file",
                        type=CommandLineParser.valid_file, metavar="config", required=True)

    # Metadata
    parser.add_argument("--in_metadata", help="Project metadata", metavar="in_metadata", required=True)
    # Output metadata
    parser.add_argument("--out_metadata", help="Output metadata", metavar="output_metadata", required=True)
    # Log file
    parser.add_argument("--log_file", help="Log file", metavar="log_file", required=True)

    in_args = parser.parse_args()

    RESULTS = main(in_args)

    print(RESULTS)
    