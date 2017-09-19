#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import print_function

import os.path
import argparse
import sys

import collections
# Required for ReadTheDocs
from functools import wraps # pylint: disable=unused-import

from basic_modules.workflow import Workflow

from tool.tb_normalize import tbNormalizeTool

# ------------------------------------------------------------------------------
class tadbit_normalize(Workflow):
    
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
        if configuration is None:
            configuration = {}

        self.configuration.update(convert_from_unicode(configuration))
        
        self.configuration['workdir'] = os.path.abspath('tests/data/tmp/') #Not clear where this will be passed


    def run(self, input_files,metadata, output_files):
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
            "PROCESS MAP - FILES PASSED TO TOOLS:",
            remap(input_files, "reads1", "reads2", "ref_genome")
        )
        
        bamin = convert_from_unicode(input_files['bamin'])
        input_metadata = remap(self.configuration, "min_perc","workdir", "max_perc")
        m_results_meta = {}
        
        tn = tbNormalizeTool()
        tn_files, tn_meta = tn.run([bamin], [], input_metadata)
        
        m_results_meta['normalize'] = tn_meta
        m_results_meta['normalize']['error'] = ''
        
        m_results_files = {}
        m_results_files["hic_biases"] = tn_files[0]+'.cpickle'
        
         # List of files to get saved
        print("TADBIT RESULTS:", m_results_files)
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
    elif isinstance(data, collections.Mapping):
        return dict(map(convert_from_unicode, data.iteritems()))
    elif isinstance(data, collections.Iterable):
        return type(data)(map(convert_from_unicode, data))
    else:
        return data
# ------------------------------------------------------------------------------

def main(input_files, output_files, input_metadata):
    """
    Main function
    -------------

    This function launches the app.
    """

    # import pprint  # Pretty print - module for dictionary fancy printing

    # 1. Instantiate and launch the App
    print("1. Instantiate and launch the App")
    from apps.workflowapp import WorkflowApp
    app = WorkflowApp()
    result = app.launch(tadbit_normalize, input_files, input_metadata, output_files,
                        {})

    # 2. The App has finished
    print("2. Execution finished")
    print(result)
    return result

def main_json():
    """
    Alternative main function
    -------------
    This function launches the app using configuration written in
    two json files: config.json and input_metadata.json.
    """
    # 1. Instantiate and launch the App
    print("1. Instantiate and launch the App")
    from apps.jsonapp import JSONApp
    app = JSONApp()
    root_path = os.path.dirname(os.path.abspath(__file__))
    result = app.launch(tadbit_normalize,
                        root_path,
                        "tests/json/config_tadbit_normalize.json",
                        "tests/json/input_tadbit_normalize.json")

    # 2. The App has finished
    print("2. Execution finished; see " + root_path + "/results.json")
    print(result)

    return result

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    sys._run_from_cmdl = True

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="TADbit normalize")
    PARSER.add_argument("--bamin", help="path to a TADbit-generated BAM file with all reads")
    PARSER.add_argument("--min_perc", help="lower percentile from which consider bins as good.")
    PARSER.add_argument("--max_perc", help="upper percentile until which consider bins as good.")
    PARSER.add_argument("--json",
                        help="Use defined JSON config files",
                        action='store_const', const=True, default=False)
    
    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    
    BAMIN = ARGS.bamin
    MIN_PERC = ARGS.min_perc
    MAX_PERC = ARGS.max_perc
    JSON_CONFIG = ARGS.json
    
    METADATA = {
        'min_perc' : MIN_PERC,
        'max_perc' : MAX_PERC
    }
    FILES = [
        BAMIN
    ]

    

    if JSON_CONFIG is True:
        RESULTS = main_json()
    else:
        # 3. Instantiate and launch the App
        RESULTS = main(FILES, [], METADATA)

    print(RESULTS)
