#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import print_function

import os.path
import argparse
import sys
import tarfile

import collections
# Required for ReadTheDocs
from functools import wraps # pylint: disable=unused-import

from basic_modules.workflow import Workflow

from tool.tb_full_mapping import tbFullMappingTool
from tool.tb_parse_mapping import tbParseMappingTool
from tool.tb_filter import tbFilterTool

# ------------------------------------------------------------------------------
class tadbit_map_parse_filter(Workflow):
    
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
        if 'windows' in self.configuration:
            if self.configuration['windows']:
                w1 = self.configuration['windows'].split(" ")
                self.configuration['windows'] = [tuple(map(int, x.split(':'))) for x in w1]
            else:
                self.configuration['windows'] = ''

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
            List of locations for the output bam files
        """
        
        print(
            "PROCESS MAP - FILES PASSED TO TOOLS:",
            remap(input_files, "reads1", "reads2", "ref_genome")
        )
        
        genome_fa = convert_from_unicode(input_files['ref_genome'])
        genome_gem = convert_from_unicode(input_files['ref_genome_gem'])
        fastq_file_1 = convert_from_unicode(input_files['reads1'])
        fastq_file_2 = convert_from_unicode(input_files['reads2'])
        input_metadata = remap(self.configuration, "iterative_mapping","workdir", "windows",rest_enzyme="enzyme_name")
        input_metadata['quality_plot'] = True        
        m_results_meta = {}
        
        tfm1 = tbFullMappingTool()
        tfm1_files, tfm1_meta = tfm1.run([genome_gem, fastq_file_1], [], input_metadata)
        
        m_results_meta['map1'] = tfm1_meta
        m_results_meta['map1']['error'] = ''
        
        tfm2 = tbFullMappingTool()
        tfm2_files, tfm2_meta = tfm2.run([genome_gem, fastq_file_2], [], input_metadata)
        
        m_results_meta['map2'] = tfm2_meta
        m_results_meta['map2']['error'] = ''
        
        tpm = tbParseMappingTool()
        files = [genome_fa] + tfm1_files[:-2] + tfm2_files[:-2]

        input_metadata = remap(self.configuration, "chromosomes","workdir",rest_enzyme="enzyme_name")                        
        input_metadata['mapping'] = [tfm1_meta['func'], tfm2_meta['func']]
        input_metadata['expt_name'] = 'vre' 
        
         
        print("TB MAPPED FILES:", files)
        print("TB PARSE METADATA:", input_metadata)
        tpm_files, tpm_meta = tpm.run(files, [], input_metadata)
  
        m_results_meta['parse'] = tpm_meta
        m_results_meta['parse']['error'] = ''
          
        print("TB PARSED FILES:", tpm_files)
         
        input_metadata = remap(self.configuration, "chromosomes","workdir",'self-circle','dangling-end', 'error', 
                               'extra dangling-end','too close from REs', 'too short','too large', 
                               'over-represented' ,'duplicated', 'random breaks','min_dist_RE','min_fragment_size','max_fragment_size')                
        input_metadata['expt_name'] = 'vre'
        input_metadata['outbam'] = 'vre_filtered_reads'
        input_metadata['custom_filter'] = True
        input_metadata['histogram'] = True
         
        tbf = tbFilterTool()
        tf_files, tf_meta = tbf.run(tpm_files, [], input_metadata)
   
        m_results_meta['filter'] = tf_meta
        m_results_meta['filter']['error'] = ''
        #adjlist_loc = f2a.save_hic_data()
         
        print("TB FILTER FILES:", tf_files[0])
          
        m_results_files = {}
        m_results_files["paired_reads"] = tf_files[0]+'.bam'
        m_results_files["map_parse_filter_stats"] = input_metadata['workdir']+"/map_parse_filter_stats.tar.gz"
         
        with tarfile.open(m_results_files["map_parse_filter_stats"], "w:gz") as tar:
            tar.add(tfm1_files[-1],arcname=os.path.basename(tfm1_files[-1]))
            tar.add(tfm1_files[-2],arcname=os.path.basename(tfm1_files[-2]))
            tar.add(tfm2_files[-1],arcname=os.path.basename(tfm2_files[-1]))
            tar.add(tfm2_files[-2],arcname=os.path.basename(tfm2_files[-2]))
            tar.add(tf_files[-1],arcname=os.path.basename(tf_files[-1]))
            tar.add(tf_files[-2],arcname=os.path.basename(tf_files[-2]))
        
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
    result = app.launch(tadbit_map_parse_filter, input_files, input_metadata, output_files,
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
    result = app.launch(tadbit_map_parse_filter,
                        root_path,
                        "tests/json/config_tadbit_map_parse_filter.json",
                        "tests/json/input_tadbit_map_parse_filter.json")

    # 2. The App has finished
    print("2. Execution finished; see " + root_path + "/results.json")
    print(result)

    return result

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    sys._run_from_cmdl = True

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="TADbit map")
    PARSER.add_argument("--genome_gem", help="Genome assembly GEM file")
    PARSER.add_argument("--file1", help="Location of FASTQ file 1")
    PARSER.add_argument("--file2", help="Location of FASTQ file 2")
    PARSER.add_argument("--rest_enzyme", help="Enzyme used to digest the DNA")
    PARSER.add_argument("--iterative_mapping", help="Iterative or fragment based mapping", default="false")
    PARSER.add_argument(
        "--windows",
        help="FASTQ windowing - start locations",
        default=None)
    PARSER.add_argument("--workdir",help="Working directory",default='')
    PARSER.add_argument("--json",
                        help="Use defined JSON config files",
                        action='store_const', const=True, default=False)
    

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    
    GENOME_GEM = ARGS.genome_gem
    FASTQ_01_FILE = ARGS.file1
    FASTQ_02_FILE = ARGS.file2
    REST_ENZYME = ARGS.rest_enzyme
    WINDOWS = ARGS.windows
    ITERATIVE_MAPPING = ARGS.iterative_mapping
    WORKDIR = ARGS.workdir
    JSON_CONFIG = ARGS.json
    
    print("ENZYME_NAME:", REST_ENZYME)

    
    METADATA = {
        'enzyme_name' : REST_ENZYME,
        'windows' : None,
        'iterative_mapping' : ITERATIVE_MAPPING,
        'workdir': WORKDIR
    }
    
    if WINDOWS is not None:
        W1 = WINDOWS.split(" ")
        WINDOWSARG = [tuple(map(int, x.split(':'))) for x in W1]
        print("WINDOWS:", WINDOWSARG, WINDOWS)
        METADATA['windows'] = WINDOWSARG
    
    FILES = [
        GENOME_GEM,
        FASTQ_01_FILE,
        FASTQ_02_FILE
    ]

    

    if JSON_CONFIG is True:
        RESULTS = main_json()
    else:
        # 3. Instantiate and launch the App
        RESULTS = main(FILES, [], METADATA)

    print(RESULTS)
