#!/usr/bin/env python
#!/home/pmes/.pyenv/shims/python

# -*- coding: utf-8 -*-

from __future__ import print_function

import os.path
import argparse
import sys
import tarfile
import multiprocessing
import json
import urllib2
import shutil
import collections

from random import random
from string import ascii_letters as letters

# Required for ReadTheDocs
from functools import wraps # pylint: disable=unused-import

from basic_modules.workflow import Workflow
from basic_modules.metadata import Metadata

from tool.tb_model import tbModelTool

if '/opt/COMPSs/Bindings/python' in sys.path:
    sys.path.pop(sys.path.index('/opt/COMPSs/Bindings/python'))

class ResultObj(dict):

    error = False
    metadata = {}

    def __init__(self, error, metadata):
        self.error = error
        self.metadata = metadata

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
class tadbit_model(Workflow):
    """
    Wrapper for the VRE form TADbit model.
    It has two main sections:
        - looks for optimal parameters for modeling a region
        - models a region for a given optimal parameters
    .
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

        self.configuration["optimize_only"] = not "generation:num_mod_comp" in self.configuration
        if "optimization:max_dist" in self.configuration and not self.configuration["optimize_only"]:
            del self.configuration["optimization:max_dist"]
            del self.configuration["optimization:upper_bound"]
            del self.configuration["optimization:lower_bound"]
            del self.configuration["optimization:cutoff"]
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
            List of locations for the output bam files
        """

        print(
            "PROCESS MODEL - FILES PASSED TO TOOLS:",
            remap(input_files, "hic_contacts_matrix_norm")
        )

        m_results_meta = {}
        m_results_files = {}

        input_metadata = remap(self.configuration, "optimize_only", "gen_pos_chrom_name", "resolution", "gen_pos_begin",
                               "gen_pos_end", "max_dist", "upper_bound", "lower_bound", "cutoff", "workdir", "ncpus")
        in_files = [convert_from_unicode(input_files['hic_contacts_matrix_norm'])]
        input_metadata["species"] = "Unknown"
        input_metadata["assembly"] = "Unknown"
        if "assembly" in metadata['hic_contacts_matrix_norm'].meta_data:
            input_metadata["assembly"] = metadata['hic_contacts_matrix_norm'].meta_data["assembly"]
        if metadata['hic_contacts_matrix_norm'].taxon_id:
            dt = json.load(urllib2.urlopen("http://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/tax-id/"+str(metadata['hic_contacts_matrix_norm'].taxon_id)))
            input_metadata["species"] = dt['scientificName']

        input_metadata["num_mod_comp"] = self.configuration["num_mod_comp"]
        input_metadata["num_mod_keep"] = self.configuration["num_mod_keep"]

        tm = tbModelTool()
        tm_files, tm_meta = tm.run(in_files, [], input_metadata)

        m_results_files["modeling_stats"] = self.configuration['project']+"/model_stats.tar.gz"

        tar = tarfile.open(m_results_files["modeling_stats"], "w:gz")
        tar.add(tm_files[0], arcname='modeling_files_and_stats')
        tar.close()

        if not self.configuration["optimize_only"]:
            m_results_files["tadkit_models"] = self.configuration['project']+"/"+os.path.basename(tm_files[1])
            os.rename(tm_files[1], m_results_files["tadkit_models"])
            m_results_meta["tadkit_models"] = Metadata(
                data_type="chromatin_3dmodel_ensemble",
                file_type="JSON",
                file_path=m_results_files["tadkit_models"],
                sources=in_files,
                meta_data={
                    "description": "Ensemble of chromatin 3D structures",
                    "visible": True,
                    "assembly": input_metadata["assembly"]
                },
                taxon_id=metadata['hic_contacts_matrix_norm'].taxon_id)

        # List of files to get saved
        print("TADBIT RESULTS:", m_results_files)

        m_results_meta["modeling_stats"] = Metadata(
            data_type="tool_statistics",
            file_type="TAR",
            file_path=m_results_files["modeling_stats"],
            sources=in_files,
            meta_data={
                "description": "TADbit modeling statistics and result files",
                "visible": True
            })

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
    result = app.launch(tadbit_model,
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

def make_absolute_path(files, root):
    """Make paths absolute."""
    for role, path in files.items():
        files[role] = os.path.join(root, path)
    return files

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
