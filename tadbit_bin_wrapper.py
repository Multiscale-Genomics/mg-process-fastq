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
import urllib2
import collections
import tarfile

from random import random
from string import ascii_letters as letters

# Required for ReadTheDocs
from functools import wraps # pylint: disable=unused-import

from basic_modules.workflow import Workflow
from basic_modules.metadata import Metadata

from tool.tb_bin import tbBinTool

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
class tadbit_bin(Workflow):
    """
    Wrapper for the VRE form TADbit bin.
    It extracts a section of a matrix from a BAM file.
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
            "PROCESS BIN - FILES PASSED TO TOOLS:",
            remap(input_files, "bamin")
        )
        m_results_files = {}
        m_results_meta = {}
        #hic_data = load_hic_data_from_reads('/home/dcastillo/workspace/vre/mg-process-fastq-tadbit/tests/data/raw_None:0-13381_10kb.abc', resolution=10000)
        #exp = Experiment("vre", resolution=10000, hic_data=hic_data)

        input_metadata = remap(self.configuration, "resolution", "workdir", "ncpus")
        if "coord1" in self.configuration:
            input_metadata["coord1"] = self.configuration["coord1"]
        if "coord2" in self.configuration:
            input_metadata["coord2"] = self.configuration["coord2"]
        input_metadata["norm"] = ['raw']
        in_files = [convert_from_unicode(input_files['bamin'])]
        if 'hic_biases' in input_files:
            in_files.append(convert_from_unicode(input_files['hic_biases']))
            input_metadata["norm"] = ['raw', 'norm']
        input_metadata["species"] = "Unknown"
        input_metadata["assembly"] = "Unknown"
        if "assembly" in metadata['bamin'].meta_data:
            input_metadata["assembly"] = metadata['bamin'].meta_data["assembly"]
        if metadata['bamin'].taxon_id:
            dt = json.load(urllib2.urlopen("http://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/tax-id/"+str(metadata['bamin'].taxon_id)))
            input_metadata["species"] = dt['scientificName']
        #hic_data = HiC_data((), len(bins_dict), sections, bins_dict, resolution=int(input_metadata['resolution']))
        tb = tbBinTool()
        tb_files, tb_meta = tb.run(in_files, [], input_metadata)

        m_results_files["bin_stats"] = self.configuration['project']+"/bin_stats.tar.gz"
        m_results_files["hic_contacts_matrix_raw"] = self.configuration['project']+"/"+os.path.basename(tb_files[0])
        os.rename(tb_files[0], m_results_files["hic_contacts_matrix_raw"])
        if len(input_metadata["norm"]) > 1:
            m_results_files["hic_contacts_matrix_norm"] = self.configuration['project']+"/"+os.path.basename(tb_files[2])
            os.rename(tb_files[2], m_results_files["hic_contacts_matrix_norm"])

        with tarfile.open(m_results_files["bin_stats"], "w:gz") as tar:
            tar.add(tb_files[1], arcname=os.path.basename(tb_files[1]))
            if len(input_metadata["norm"]) > 1:
                tar.add(tb_files[3], arcname=os.path.basename(tb_files[3]))
            tar.add(tb_files[-1], arcname=os.path.basename(tb_files[-1]))
        if os.path.isdir(tb_files[1]):
            clean_temps(tb_files[1])
        if len(input_metadata["norm"]) > 1:
            if os.path.isdir(tb_files[3]):
                clean_temps(tb_files[3])

        # List of files to get saved
        print("TADBIT RESULTS:", m_results_files)

        m_results_meta["hic_contacts_matrix_raw"] = Metadata(
            data_type="hic_contacts_matrix",
            file_type="TXT",
            file_path=m_results_files["hic_contacts_matrix_raw"],
            sources=in_files,
            meta_data={
                "description": "HiC contact matrix raw",
                "visible": True,
                "assembly": convert_from_unicode(metadata['bamin'].meta_data['assembly']),
                "norm" : 'raw'
            },
            taxon_id=metadata['bamin'].taxon_id)
        m_results_meta["bin_stats"] = Metadata(
            data_type="tool_statistics",
            file_type="TAR",
            file_path=m_results_files["bin_stats"],
            sources=in_files,
            meta_data={
                "description": "TADbit HiC matrices in png format",
                "visible": False
            })
        if len(input_metadata["norm"]) > 1:
            m_results_meta["hic_contacts_matrix_norm"] = Metadata(
                data_type="hic_contacts_matrix",
                file_type="TXT",
                file_path=m_results_files["hic_contacts_matrix_norm"],
                sources=in_files,
                meta_data={
                    "description": "HiC contact matrix normalized",
                    "visible": True,
                    "assembly": convert_from_unicode(metadata['bamin'].meta_data['assembly']),
                    "norm" : 'norm'
                },
                taxon_id=metadata['bamin'].taxon_id)
        
        m_results_files["tadkit_matrix"] = self.configuration['project']+"/"+os.path.basename(tb_files[-1])
        os.rename(tb_files[-1], m_results_files["tadkit_matrix"])
        m_results_meta["tadkit_matrix"] = Metadata(
            data_type="chromatin_3dmodel_ensemble",
            file_type="JSON",
            file_path=m_results_files["tadkit_matrix"],
            sources=in_files,
            meta_data={
                "description": "Ensemble of chromatin 3D structures",
                "visible": True,
                "assembly": input_metadata["assembly"]
            },
            taxon_id=metadata['bamin'].taxon_id)
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
    result = app.launch(tadbit_bin,
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
    