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
import glob, os
from subprocess import CalledProcessError, PIPE, Popen
from pytadbit                     import load_hic_data_from_bam
from pytadbit                            import Chromosome, HiC_data
import numpy as np

from cPickle import load, dump
from _pytest.tmpdir import tmpdir

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

class tbModelTool(Tool):
    """
    Tool for normalizing an adjacency matrix
    """

    def __init__(self):
        """
        Init function
        """
        print("TADbit - Modeling")
        Tool.__init__(self)

    @task(hic_contacts_matrix_norm=FILE_IN, resolution=IN, gen_pos_chrom_name=IN , gen_pos_begin=IN,
                                gen_pos_end=IN, num_mod_comp=IN, num_mod_keep=IN,
                                max_dist=IN, upper_bound=IN, lower_bound=IN, cutoff=IN, workdir=IN)
    # @constraint(ProcessorCoreCount=16)
    def tb_model(self, optimize_only, hic_contacts_matrix_norm, resolution, gen_pos_chrom_name , gen_pos_begin,
                                gen_pos_end, num_mod_comp, num_mod_keep,
                                max_dist, upper_bound, lower_bound, cutoff, workdir,metadata,
                                ncpus=1):
        """
        Function to normalize to a given resolution the Hi-C
        matrix

        Parameters
        ----------
        optimize_only: bool
            True if only optimize, False for computing the models and stats
        hic_contacts_matrix_norm : str
            Location of the tab-separated normalized matrix
        resolution : str
            Resolution of the Hi-C
        gen_pos_chrom_name : str
            Coordinates of the genomic region to model.                                                                                                
        gen_pos_begin : int
            Genomic coordinate from which to start modeling. 
        gen_pos_end : int
            Genomic coordinate where to end modeling.
        num_mod_comp : int
            Number of models to compute for each optimization step.
        num_mod_comp : int
            Number of models to keep.
        max_dist : str
            Range of numbers for optimal maxdist parameter, i.e. 400:1000:100; or just a single number e.g. 800; or a list of numbers e.g. 400 600 800 1000.
        upper_bound : int
            Range of numbers for optimal upfreq parameter, i.e. 0:1.2:0.3; or just a single number e.g. 0.8; or a list of numbers e.g. 0.1 0.3 0.5 0.9.
        lower_bound : int
            Range of numbers for optimal low parameter, i.e. -1.2:0:0.3; or just a single number e.g. -0.8; or a list of numbers e.g. -0.1 -0.3 -0.5 -0.9.
        cutoff : str
            Range of numbers for optimal cutoff distance. Cutoff is computed based on the resolution. This cutoff distance is calculated taking as reference the diameter 
            of a modeled particle in the 3D model. i.e. 1.5:2.5:0.5; or just a single number e.g. 2; or a list of numbers e.g. 2 2.5.
        workdir : str
            Location of working directory
        ncpus : str
            Number of cpus to use
                
        Returns
        -------
        tadkit_models : str
            Location of TADkit json file
        modeling_stats : str
            Location of the folder with the modeling files and stats

        """
        #chr_hic_data = read_matrix(matrix_file, resolution=int(resolution))

        print("TB MODELING:",hic_contacts_matrix_norm, resolution, gen_pos_chrom_name , gen_pos_begin,
                                gen_pos_end, num_mod_comp, num_mod_keep,
                                max_dist, upper_bound, lower_bound, cutoff, workdir)

        try:
            beg = int(float(gen_pos_begin) / int(resolution))
            end = int(float(gen_pos_end) / int(resolution))
            if end - beg <= 2:
                raise Exception('"beg" and "end" parameter should be given in ' +
                                'genomic coordinates, not bin')
        except TypeError:
            pass
        
        name = '{0}_{1}_{2}'.format(gen_pos_chrom_name, beg, end)
        if not os.path.exists(os.path.join(workdir, name)):
            os.makedirs(os.path.join(workdir, name))
        
        sqr_mat = convert_matrix(hic_contacts_matrix_norm)
        hic_contacts_matrix_square = os.path.join(workdir, os.path.basename(hic_contacts_matrix_norm)+'.sqr')
        np.savetxt(hic_contacts_matrix_square, sqr_mat)
        
        _cmd = [
                'model_and_analyze.py',
            '--nmodels_opt',str(num_mod_comp),
            '--nkeep_opt',str(num_mod_keep),
            '--norm', hic_contacts_matrix_square,
            '--res', resolution,
            '--crm', gen_pos_chrom_name,
            '--beg',str(gen_pos_begin),
            '--end',str(gen_pos_end),
            '--maxdist',max_dist,
            '--upfreq', upper_bound,
            '--lowfreq='+lower_bound,
            '--dcutoff', cutoff,
            '--ncpus', str(ncpus),
            '--assembly',metadata["assembly"],
            '--species',metadata["species"]
            ]
        if optimize_only:
            _cmd.append('--optimize_only')
            _cmd.append('--outdir')
            _cmd.append(workdir)
            
        else:

            #crm = load_hic_data(bamin, int(resolution), gen_pos_chrom_name, metadata["species"], metadata["assembly"], workdir, xbias=biases, ncpus=ncpus)
            chrom_name = gen_pos_chrom_name.split(':')[0]   
            # Start reading the data
            crm = Chromosome(chrom_name, species=(
                metadata["species"].split('_')[0].capitalize() + metadata["species"].split('_')[1]
                                  if '_' in metadata["species"] else metadata["species"]), assembly = metadata["assembly"]) # Create chromosome object
            crm.add_experiment(
                os.path.split(hic_contacts_matrix_norm)[-1], exp_type='Hi-C', 
                resolution=int(resolution),
                norm_data=hic_contacts_matrix_square)
            exp = crm.experiments[-1]
            if beg > crm.experiments[-1].size*int(resolution):
                raise Exception('ERROR: beg parameter is larger than chromosome size.')
            if end > crm.experiments[-1].size*int(resolution):
                raise Exception('WARNING: end parameter is larger than chromosome ' +
                             'size. Setting end to %s.\n' % (crm.experiments[-1].size *
                                                             int(resolution)))
                end = crm.experiments[-1].size
            STD_CONFIG = {
              'kforce'    : 5,
              'maxdist'   : float(max_dist), 
              'upfreq'    : float(upper_bound), 
              'lowfreq'   : float(lower_bound),
              'dcutoff'   : float(cutoff),
              'scale'     : 0.01
              }
             
            outfile = os.path.join(workdir, name, name + '.models')
            models = exp.model_region(start=beg, end=end, n_models=int(num_mod_comp), n_keep=int(num_mod_keep),
                     n_cpus=ncpus, verbose=1, keep_all=False, close_bins=1,
                     config=STD_CONFIG)
            #models._config['dcutoff']
            models.save_models(outfile)
            _cmd.append('--analyze_only')
            _cmd.append('--outdir')
            _cmd.append(workdir)
            
        output_metadata = {}
        
        out, err = Popen(_cmd, stdout=PIPE, stderr=PIPE).communicate()
        print(out)
        print(err)

        output_files = [os.path.join(workdir, name)]
        
        if not optimize_only:
            os.chdir(os.path.join(workdir, name))
            for fl in glob.glob("*.json"):
                output_files.append(fl)
                break
            
        
        
        return (output_files, output_metadata)

    def run(self, input_files, output_files, metadata=None):
        """
        The main function for the normalization of the Hi-C matrix to a given resolution

        Parameters
        ----------
        input_files : list
            hic_contacts_matrix_norm : str
                Location of the tab-separated normalized matrix
        metadata : dict
            optimize_only: bool
                True if only optimize, False for computing the models and stats
            gen_pos_chrom_name : str
                Coordinates of the genomic region to model.   
            resolution : str
                Resolution of the Hi-C                                                                                             
            gen_pos_begin : int
                Genomic coordinate from which to start modeling. 
            gen_pos_end : int
                Genomic coordinate where to end modeling.
            num_mod_comp : int
                Number of models to compute for each optimization step.
            num_mod_comp : int
                Number of models to keep.
            max_dist : str
                Range of numbers for optimal maxdist parameter, i.e. 400:1000:100; or just a single number e.g. 800; or a list of numbers e.g. 400 600 800 1000.
            upper_bound : int
                Range of numbers for optimal upfreq parameter, i.e. 0:1.2:0.3; or just a single number e.g. 0.8; or a list of numbers e.g. 0.1 0.3 0.5 0.9.
            lower_bound : int
                Range of numbers for optimal low parameter, i.e. -1.2:0:0.3; or just a single number e.g. -0.8; or a list of numbers e.g. -0.1 -0.3 -0.5 -0.9.
            cutoff : str
                Range of numbers for optimal cutoff distance. Cutoff is computed based on the resolution. This cutoff distance is calculated taking as reference the diameter 
                of a modeled particle in the 3D model. i.e. 1.5:2.5:0.5; or just a single number e.g. 2; or a list of numbers e.g. 2 2.5.
            workdir : str
                Location of working directory
            ncpus : str
                Number of cpus to use

        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects

        """

        hic_contacts_matrix_norm = input_files[0]
        
        
        ncpus=1
        if 'ncpus' in metadata:
            ncpus = metadata['ncpus']
        
        optimize_only = metadata["optimize_only"]
        gen_pos_chrom_name = metadata['gen_pos_chrom_name']
        resolution = metadata['resolution']
        gen_pos_begin = metadata['gen_pos_begin']
        gen_pos_end = metadata['gen_pos_end']
        num_mod_comp = metadata['num_mod_comp']
        num_mod_keep = metadata['num_mod_keep']
        max_dist = metadata['max_dist']
        upper_bound = metadata['upper_bound']
        lower_bound = metadata['lower_bound']
        cutoff = metadata['cutoff']
            
        root_name = hic_contacts_matrix_norm.split("/")
        if 'workdir' in metadata:
            root_name = metadata['workdir']

        project_metadata = {}
        project_metadata["species"] = metadata["species"]
        project_metadata["assembly"] = metadata["assembly"]
        
        # input and output share most metadata
        
        
        output_files, output_metadata = self.tb_model(optimize_only, hic_contacts_matrix_norm, resolution, gen_pos_chrom_name , gen_pos_begin,
                                gen_pos_end, num_mod_comp, num_mod_keep,
                                max_dist, upper_bound, lower_bound, cutoff, root_name,
                                project_metadata, ncpus)

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
def convert_matrix(hic_matrix_path):
    
    fh = open(hic_matrix_path)
    size = 0
    for line in fh:
        if not line.startswith('# CRM'):
            break
    offset_pos = line.split(' ')[1].split(':')[1].split('-')
    start = int(offset_pos[0])
    end = int(offset_pos[1])
    size = end
    
    fh.seek(0)
    
    for line in fh:
        if line.startswith('# BADS'):
            if len(line.split()) > 2:
                bads = set(map(int, line.split()[-1].split(',')))
            else:
                bads = {}
            break
        
    matrix = np.zeros((size, size))
    
    for line in fh:
        i, j, v = line.split()
        if i in bads or j in bads:
            continue
        matrix[int(i)+start][int(j)+start] = float(v)
    
    return matrix
    

def load_hic_data(xpath, resolution, gen_pos_chrom_name, species, assembly, workdir, xbias=None, ncpus=1):
    """
    Load Hi-C data
    """
    xnam = os.path.split(xpath)[-1]
    hic_raw = load_hic_data_from_bam(xpath, resolution, biases=xbias if xbias else None, tmpdir=workdir, ncpus=int(ncpus))    
    chrom_name = gen_pos_chrom_name.split(':')[0]   
    # Start reading the data
    crm = Chromosome(chrom_name, species=(
        species.split('_')[0].capitalize() + species.split('_')[1]
                          if '_' in species else species), assembly = assembly) # Create chromosome object

    crm.add_experiment(
        xnam, exp_type='Hi-C', 
        resolution=resolution,
        hic_data=hic_raw)
    if xbias:
        bias_ = load(open(xbias))
        bias = bias_['biases']
        bads = bias_['badcol']
        if bias_['resolution'] != resolution:
            raise Exception('ERROR: resolution of biases do not match to the '
                            'one wanted (%d vs %d)' % (
                                bias_['resolution'], resolution))
        
        def transform_value_norm(a, b, c):
            return c / bias[a] / bias[b]
        size = crm.experiments[-1].size
        xnorm = [HiC_data([(i + j * size, float(hic_raw[i, j]) /
                                bias[i] /
                                bias[j] * size)
                               for i in bias for j in bias if i not in bads and j not in bads], size)]
        crm.experiments[xnam]._normalization = 'visibility_factor:1'
        factor = sum(xnorm[0].values()) / (size * size)
        for n in xnorm[0]:
            xnorm[0][n] = xnorm[0][n] / factor
        crm.experiments[xnam].norm = xnorm
    
    return crm
        

    
def my_round(num, val=4):
    num = round(float(num), val)
    return str(int(num) if num == int(num) else num)