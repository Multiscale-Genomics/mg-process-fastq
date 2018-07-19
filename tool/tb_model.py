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
import glob
import os
from subprocess import PIPE, Popen

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, FILE_INOUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
    # from pycompss.api.constraint import constraint
except ImportError:
    logger.info("[Warning] Cannot import \"pycompss\" API packages.")
    logger.info("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import compss_wait_on # pylint: disable=ungrouped-imports
    #from utils.dummy_pycompss import constraint

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
        logger.info("TADbit - Modeling")
        Tool.__init__(self)

    @task(optimize_only=IN, hic_contacts_matrix_norm=FILE_IN, resolution=IN, gen_pos_chrom_name=IN,
          gen_pos_begin=IN, gen_pos_end=IN, num_mod_comp=IN, num_mod_keep=IN,
          max_dist=IN, upper_bound=IN, lower_bound=IN, cutoff=IN, workdir=IN,
          model_dir=FILE_OUT)
    # @constraint(ProcessorCoreCount=16)
    def tb_model(self, optimize_only, hic_contacts_matrix_norm, resolution, gen_pos_chrom_name, gen_pos_begin,
                 gen_pos_end, num_mod_comp, num_mod_keep,
                 max_dist, upper_bound, lower_bound, cutoff, workdir, metadata,
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
        
        logger.info("TB MODELING: {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}".format(hic_contacts_matrix_norm,
            resolution, gen_pos_chrom_name, gen_pos_begin,
            gen_pos_end, num_mod_comp, num_mod_keep,
            max_dist, upper_bound, lower_bound, cutoff, workdir))

        try:
            beg = int(float(gen_pos_begin) / int(resolution))
            end = int(float(gen_pos_end) / int(resolution))
            if end - beg <= 2:
                logger.fatal('"beg" and "end" parameter should be given in ' +
                                'genomic coordinates, not bin')
                raise Exception('"beg" and "end" parameter should be given in ' +
                                'genomic coordinates, not bin')
        except ValueError:
            pass
        except TypeError:
            pass

        name = '{0}_{1}_{2}'.format(gen_pos_chrom_name, gen_pos_begin, gen_pos_end)

        if not os.path.exists(os.path.join(workdir, name)):
            os.makedirs(os.path.join(workdir, name))

        #=======================================================================
        # _cmd = [
        #     'model_and_analyze.py',
        #     '--norm', hic_contacts_matrix_norm,
        #     '--res', resolution,
        #     '--maxdist', max_dist,
        #     '--upfreq', upper_bound,
        #     '--lowfreq='+lower_bound,
        #     '--dcutoff', cutoff,
        #     '--ncpus', str(ncpus),
        #     '--assembly', metadata["assembly"],
        #     '--species', metadata["species"],
        #     '--project', metadata["project"],
        #     '--fig_format', 'png'
        #     ]
        # if len(gen_pos_chrom_name) > 0:
        #     _cmd.append('--crm')
        #     _cmd.append(gen_pos_chrom_name)
        # if len(gen_pos_begin) > 0 and len(gen_pos_end) > 0:
        #     _cmd.append('--beg')
        #     _cmd.append(str(gen_pos_begin))
        #     _cmd.append('--end')
        #     _cmd.append(str(gen_pos_end))
        # if optimize_only:
        #     _cmd.append('--nmodels_opt')
        #     _cmd.append(str(num_mod_comp))
        #     _cmd.append('--nkeep_opt')
        #     _cmd.append(str(num_mod_keep))
        #     _cmd.append('--optimize_only')
        #     _cmd.append('--outdir')
        #     _cmd.append(workdir)
        # else:
        #     _cmd.append('--nmodels_opt')
        #     _cmd.append('0')
        #     _cmd.append('--nkeep_opt')
        #     _cmd.append('0')
        #     _cmd.append('--nmodels_mod')
        #     _cmd.append(str(num_mod_comp))
        #     _cmd.append('--nkeep_mod')
        #     _cmd.append(str(num_mod_keep))
        #     _cmd.append('--outdir')
        #     _cmd.append(workdir)
        #=======================================================================

        _cmd = [
            'tadbit', 'model',
            '--workdir', os.path.join(workdir, name),
            '--input_matrix', hic_contacts_matrix_norm,
            '--reso', resolution,
            '--maxdist', max_dist,
            '--upfreq', upper_bound,
            '--lowfreq='+lower_bound,
            '--dcutoff', cutoff,
            '--cpu', str(ncpus),
            '--assembly', metadata["assembly"],
            '--species', metadata["species"],
            '--project', metadata["project"],
            '--fig_format', 'png'
            ]
        if len(gen_pos_chrom_name) > 0:
            _cmd.append('--crm')
            _cmd.append(gen_pos_chrom_name)
        if gen_pos_begin and gen_pos_end:
            _cmd.append('--beg')
            _cmd.append(str(gen_pos_begin))
            _cmd.append('--end')
            _cmd.append(str(gen_pos_end))
        _cmd.append('--nmodels')
        _cmd.append(str(num_mod_comp))
        _cmd.append('--nkeep')
        _cmd.append(str(num_mod_keep))
        
        if optimize_only:
            _cmd.append('--optimize')
        else:
            _cmd.append('--model')
            _cmd.append('--force')
            _cmd.append('--analyze')
            
        output_metadata = {}
        logger.info(' '.join(_cmd))
        out, err = Popen(_cmd, stdout=PIPE, stderr=PIPE).communicate()
        logger.info(out)
        logger.info(err)

        output_folder = os.listdir(os.path.join(workdir, name, '06_model'))[0]
        output_files = [os.path.join(workdir, name,'06_model',output_folder)]

        if not optimize_only:
            os.chdir(os.path.join(workdir, name,'06_model',output_folder))
            for fl in glob.glob("*.json"):
                output_files.append(fl)
                break
            for fl in glob.glob("*optimal_params*"):
                os.unlink(fl)

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

        ncpus = 1
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

        root_name = os.path.dirname(os.path.abspath(hic_contacts_matrix_norm))
        if 'workdir' in metadata:
            root_name = metadata['workdir']

        project_metadata = {}
        project_metadata["species"] = metadata["species"]
        project_metadata["assembly"] = metadata["assembly"]
        project_metadata["project"] = os.path.basename(os.path.normpath(metadata["project"]))

        # input and output share most metadata

        output_files, output_metadata = self.tb_model(optimize_only, hic_contacts_matrix_norm, resolution, gen_pos_chrom_name, gen_pos_begin,
                                                      gen_pos_end, num_mod_comp, num_mod_keep,
                                                      max_dist, upper_bound, lower_bound, cutoff, root_name,
                                                      project_metadata, ncpus)

        return (output_files, output_metadata)
