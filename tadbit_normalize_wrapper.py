#!/usr/bin/env python
#!/home/pmes/.pyenv/shims/python

# -*- coding: utf-8 -*-

from __future__ import print_function

import os.path
import argparse
import sys
import json
import multiprocessing
from random import random
from pysam                                import AlignmentFile
from string import ascii_letters as letters

import collections
import tarfile
# Required for ReadTheDocs
from functools import wraps # pylint: disable=unused-import

from basic_modules.workflow import Workflow
from basic_modules.metadata import Metadata

from tool.tb_normalize import tbNormalizeTool

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
            "PROCESS NORMALIZE - FILES PASSED TO TOOLS:",
            remap(input_files, "bamin")
        )
        
        try:
            bamin = convert_from_unicode(input_files['bamin'])
            input_metadata = remap(self.configuration, "resolution","min_perc","workdir", "max_perc", "ncpus")
            
            bamfile = AlignmentFile(bamin, 'rb')
            if len(bamfile.references) == 1:
                input_metadata["min_count"] = "10"
            bamfile.close()
                
            m_results_meta = {}
            
            tn = tbNormalizeTool()
            tn_files, tn_meta = tn.run([bamin], [], input_metadata)
            
            m_results_files = {}
            try:
                m_results_files["hic_biases"] = self.configuration['project']+"/"+os.path.basename(tn_files[0])
                os.rename(tn_files[0], m_results_files["hic_biases"])
            except:
                pass
            
            m_results_files["normalize_stats"] = self.configuration['project']+"/normalize_stats.tar.gz"
              
            with tarfile.open(m_results_files["normalize_stats"], "w:gz") as tar:
                tar.add(tn_files[1],arcname=os.path.basename(tn_files[1]))
                if len(tn_files) > 2:
                    tar.add(tn_files[2],arcname=os.path.basename(tn_files[2]))
                
            # List of files to get saved
            print("TADBIT RESULTS:", m_results_files)

            m_results_meta["hic_biases"] = Metadata(
                    data_type="hic_biases",
                    file_type="PICKEL",
                    file_path=m_results_files["hic_biases"],
                    sources=[""],
                    meta_data={
                        "description": "HiC biases for normalization",
                        "visible": True,
                        "assembly": ""
                    },
                    taxon_id=metadata['bamin'].taxon_id)
            m_results_meta["normalize_stats"] = Metadata(
                    data_type="tool_statistics",
                    file_type="TAR",
                    file_path=m_results_files["normalize_stats"],
                    sources=[""],
                    meta_data={
                        "description": "TADbit normalize statistics",
                        "visible": True
                    })
            
         
        except Exception as e:
            m_results_meta["hic_biases"] = Metadata(
                    data_type="hic_biases",
                    file_type="PICKEL",
                    file_path=None,
                    sources=[""],
                    meta_data={
                        "description": "HiC biases for normalization",
                        "visible": True
                    })
            m_results_meta["hic_biases"].error = True
            m_results_meta["hic_biases"].exception = str(e)
        
        #cleaning
        clean_temps(self.configuration['workdir']+"/04_normalization")
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
       
def _read_config(json_path):
    """
    Read config.json to obtain:
    input_IDs: dict containing IDs of tool input files
    arguments: dict containing tool arguments
    output_files: dict containing absolute paths of tool outputs

    For more information see the schema for config.json.
    """
    configuration = json.load(file(json_path))
    input_IDs = {}
    for input_ID in configuration["input_files"]:
        input_IDs[input_ID["name"]] = input_ID["value"]

    output_files = {}
    if "output_files" in configuration:
        for output_file in configuration["output_files"]:
            output_files[output_file["name"]] = output_file["file"]

    arguments = {}
    for argument in configuration["arguments"]:
        arguments[argument["name"]] = argument["value"]

    return input_IDs, arguments, output_files

def _read_metadata(json_path):
    """
    Read input_metadata.json to obtain input_metadata_IDs, a dict
    containing metadata on each of the tool input files,
    arranged by their ID.

    For more information see the schema for input_metadata.json.
    """
    metadata = json.load(file(json_path))
    input_metadata = {}
    for input_file in metadata:
        input_metadata[input_file["_id"]] = Metadata(
            data_type=input_file["data_type"],
            file_type=input_file["file_type"],
            file_path=input_file["file_path"],
            source_id=input_file["source_id"],
            meta_data=input_file["meta_data"],
            data_id=input_file["_id"])
    taxon_id =  metadata[0]["taxon_id"]
    return input_metadata, taxon_id

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
def make_absolute_path(files, root):
    """Make paths absolute."""
    for role, path in files.items():
        files[role] = os.path.join(root, path)
    return files
# ------------------------------------------------------------------------------

def _write_json(
                input_files, input_metadata,
                output_files, output_metadata, json_path):
    """
    Write results.json using information from input_files and output_files:
    input_files: dict containing absolute paths of input files
    input_metadata: dict containing metadata on input files
    output_files: dict containing absolute paths of output files
    output_metadata: dict containing metadata on output files

    For more information see the schema for results.json.
    """
    results = []
    for role, path in output_files.items():
        results.append({
            "name": role,
            "file_path": path,
            "data_type": output_metadata[role].data_type,
            "file_type": output_metadata[role].file_type,
            "source_id": output_metadata[role].source_id,
            "taxon_id": output_metadata[role].taxon_id,
            "meta_data": output_metadata[role].meta_data
        })
    json.dump({"output_files": results}, file(json_path, 'w'))
    return True

def main(args):
    
    from apps.jsonapp import JSONApp
    app = JSONApp()
    result = app.launch(tadbit_normalize,
                        args.config,
                        args.in_metadata,
                        args.out_metadata)
    
    
    return result
#===============================================================================
#     # 1. Instantiate and launch the App
#     print("1. Instantiate and launch the App")
#     from apps.workflowapp import WorkflowApp
#     app = WorkflowApp()
#     root_dir = args.root_dir
#     
#     print ("0) Unpack information from JSON")
#     input_IDs, arguments, output_files = _read_config(
#         args.config)
# 
#     input_metadata_IDs, taxon_id = _read_metadata(
#         args.metadata)
# 
#     # arrange by role
#     input_metadata = {}
#     for role, ID in input_IDs.items():
#         input_metadata[role] = input_metadata_IDs[ID]
# 
#     # get paths from IDs
#     input_files = {}
#     for role, metadata in input_metadata.items():
#         input_files[role] = metadata.file_path
# 
#     input_files = make_absolute_path(input_files, root_dir)
#     
#     tmp_name = ''.join([letters[int(random()*52)]for _ in xrange(5)])
#     workdir = os.path.dirname(os.path.abspath(args.out_metadata))+'/_tmp_tadbit_'+tmp_name
#     if not os.path.exists(workdir):
#         os.makedirs(workdir)
#     arguments.update({"ncpus":num_cores, "root_dir": args.root_dir, "public_dir": args.public_dir, "workdir": workdir, "taxon_id":taxon_id})
#     output_files, output_metadata = app.launch(tadbit_normalize, input_files, input_metadata, output_files, arguments, )
# 
#     print("4) Pack information to JSON")
#     #cleaning
#     clean_temps(workdir+"/04_normalization")
#     clean_temps(workdir)
#     
#     return _write_json(
#         input_files, input_metadata,
#         output_files, output_metadata,
#         args.out_metadata)
#===============================================================================
    
    
    
def clean_temps(working_path):
    """Cleans the workspace from temporal folder and scratch files"""
    for the_file in os.listdir(working_path):
        file_path = os.path.join(working_path, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            #elif os.path.isdir(file_path): shutil.rmtree(file_path)
        except:
            pass
    try:
        os.rmdir(working_path)
    except:
        pass
    print('[CLEANING] Finished')


# ------------------------------------------------------------------------------

if __name__ == "__main__":
    sys._run_from_cmdl = True

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

    args = parser.parse_args()

    RESULTS = main(args)
    
    print(RESULTS)
    
