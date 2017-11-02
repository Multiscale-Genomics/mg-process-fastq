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

import sys
from sys                             import stdout
from subprocess import CalledProcessError, PIPE, Popen
import glob, os
from cPickle                         import load

from pytadbit.parsers.hic_bam_parser import write_matrix
from pytadbit import HiC_data, Experiment, Chromosome
from pytadbit.parsers.hic_parser import load_hic_data_from_bam
from pytadbit.mapping.analyze import hic_map

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, FILE_INOUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
    # from pycompss.api.constraint import constraint
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_IN, FILE_OUT, FILE_INOUT, IN
    from dummy_pycompss import task
    from dummy_pycompss import compss_wait_on
    #from dummy_pycompss import constraint

from basic_modules.tool import Tool

# ------------------------------------------------------------------------------

class tbBinTool(Tool):
    """
    Tool for binning an adjacency matrix
    """

    def __init__(self):
        """
        Init function
        """
        print("TADbit - Bin")
        Tool.__init__(self)

    @task(bamin=FILE_IN, biases=FILE_IN, resolution=IN, c=IN, c2=IN, norm=IN, workdir=IN)
    # @constraint(ProcessorCoreCount=16)
    def tb_bin(self, bamin, biases, resolution, coord1, coord2, norm, workdir, ncpus="1"):
        """
        Function to bin to a given resolution the Hi-C
        matrix

        Parameters
        ----------
        bamin : str
            Location of the tadbit bam paired reads
        biases : str
            Location of the pickle hic biases
        resolution : int
            Resolution of the Hi-C
        coord1 : str
            Coordinate of the region to retrieve. By default all genome,
            arguments can be either one chromosome name, or the coordinate in
            the form: "-c chr3:110000000-120000000"                                             
        coord2 : str
            Coordinate of a second region to retrieve the matrix in the
            intersection with the first region. 
        norm : list
            [['raw']] normalization(s) to apply. Order matters. Choices: [norm,
            decay, raw]
        workdir : str
            Location of working directory
        ncpus : int
            Number of cpus to use
                
        Returns
        -------
        hic_contacts_matrix_raw : str
            Location of HiC raw matrix in text format
        hic_contacts_matrix_nrm : str
            Location of HiC normalized matrix in text format
        hic_contacts_matrix_raw_fig : str
            Location of HiC raw matrix in png format
        hic_contacts_matrix_norm_fig : str
            Location of HiC normalized matrix in png format

        """
        print("TB BIN:",bamin, resolution, workdir)
            
        output_metadata = {}
        
        if coord2 and not coord1:
            coord1, coord2 = coord2, coord1
    
        if not coord1:
            region1 = None
            start1  = None
            end1    = None
            region2 = None
            start2  = None
            end2    = None
        else:
            try:
                crm1, pos1   = coord1.split(':')
                start1, end1 = pos1.split('-')
                region1 = crm1
                start1  = int(start1)
                end1    = int(end1)
            except ValueError:
                region1 = coord1
                start1  = None
                end1    = None
            if coord2:
                try:
                    crm2, pos2   = coord2.split(':')
                    start2, end2 = pos2.split('-')
                    region2 = crm2
                    start2  = int(start2)
                    end2    = int(end2)
                except ValueError:
                    region2 = coord2
                    start2  = None
                    end2    = None
            else:
                region2 = None
                start2  = None
                end2    = None
    
        if region1:
            if region1:
                stdout.write('\nExtraction of %s' % (region1))
                if start1:
                    stdout.write(':%s-%s' % (start1, end1))
                else:
                    stdout.write(' (full chromosome)')
                if region2:
                    stdout.write(' intersection with %s' % (region2))
                    if start2:
                        stdout.write(':%s-%s\n' % (start2, end2))
                    else:
                        stdout.write(' (full chromosome)\n')
                else:
                    stdout.write('\n')
        else:
            stdout.write('\nExtraction of full genome\n')
            
        out_files = write_matrix(bamin, resolution,
                         load(open(biases)) if biases else None,
                         workdir, filter_exclude=[],
                         normalizations=norm,half_matrix=True,
                         region1=region1, start1=start1, end1=end1,
                         region2=region2, start2=start2, end2=end2,
                         tmpdir=workdir, append_to_tar=None, ncpus=ncpus,
                         verbose=True)
        output_files = [out_files['RAW']]
        imx = load_hic_data_from_bam(bamin, resolution, biases=biases if biases else None, ncpus=ncpus,tmpdir=workdir)
        hic_contacts_matrix_raw_fig = workdir+"/genomic_maps_raw.png"
        focus=None
        by_chrom=None
        if start1:
            focus=((start1/resolution),(end1/resolution))
        else:
            by_chrom='all'
        hic_map(imx, resolution, savefig=hic_contacts_matrix_raw_fig, normalized=False,by_chrom=by_chrom,focus=focus)
        output_files.append(hic_contacts_matrix_raw_fig)
        
        if len(norm) > 1:
            output_files.append(out_files['NRM'])
            hic_contacts_matrix_norm_fig = workdir+"/genomic_maps_nrm.png"    
            hic_map(imx, resolution, savefig=hic_contacts_matrix_norm_fig, normalized=True,by_chrom=by_chrom,focus=focus)
            output_files.append(hic_contacts_matrix_norm_fig)
        
        return (output_files, output_metadata)

    def run(self, input_files, output_files, metadata=None):
        """
        The main function to the predict TAD sites for a given resolution from
        the Hi-C matrix

        Parameters
        ----------
        input_files : list
            bamin : str
                Location of the tadbit bam paired reads
            biases : str
                Location of the pickle hic biases
        metadata : dict
            resolution : int
                Resolution of the Hi-C
            coord1 : str
                Coordinate of the region to retrieve. By default all genome,
                arguments can be either one chromosome name, or the coordinate in
                the form: "-c chr3:110000000-120000000"                                             
            coord2 : str
                Coordinate of a second region to retrieve the matrix in the
                intersection with the first region. 
            norm : str
                [['raw']] normalization(s) to apply. Order matters. Choices: [norm,
                decay, raw]
            workdir : str
                Location of working directory
            ncpus : int
                Number of cpus to use



        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects

        """

        bamin = input_files[0]
        
        if not os.path.isfile(bamin.replace('.bam','.bam.bai')):
            print('Creating bam index')
            _cmd = ['samtools', 'index', bamin]
            out, err = Popen(_cmd, stdout=PIPE, stderr=PIPE).communicate()
            print(out)
            print(err)
            
        biases = None
        if len(input_files) > 1:
            biases = input_files[1]
        resolution = 100000
        if 'resolution' in metadata:
            resolution = int(metadata['resolution'])

        coord1 = coord2 = norm = None 
        ncpus=1
        if 'ncpus' in metadata:
            ncpus = metadata['ncpus']
            
        if 'coord1' in metadata:
            coord1 = metadata['coord1']
        if 'coord2' in metadata:
            coord2 = metadata['coord2']
        if 'norm' in metadata:
            norm = metadata['norm']
        
        root_name = bamin.split("/")
        if 'workdir' in metadata:
            root_name = metadata['workdir']
        
        # input and output share most metadata
        
        
        output_files, output_metadata = self.tb_bin(bamin, biases, resolution, coord1, coord2, norm, root_name, ncpus)

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
