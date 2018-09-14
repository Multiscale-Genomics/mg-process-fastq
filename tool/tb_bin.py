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
from sys import stdout
from subprocess import PIPE, Popen
import os
from cPickle import load

from pytadbit.parsers.hic_bam_parser import write_matrix
from pytadbit import Chromosome
from pytadbit.parsers.hic_parser import load_hic_data_from_bam
from pytadbit.mapping.analyze import hic_map

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
except ImportError:
    logger.info("[Warning] Cannot import \"pycompss\" API packages.")
    logger.info("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task  # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool


# ------------------------------------------------------------------------------

class tbBinTool(Tool):  # pylint: disable=invalid-name,too-few-public-methods
    """
    Tool for binning an adjacency matrix
    """

    def __init__(self):
        """
        Init function
        """
        logger.info("TADbit - Bin")
        Tool.__init__(self)

    @task(bamin=FILE_IN, biases=FILE_IN, resolution=IN, coord1=IN, coord2=IN,
          norm=IN, workdir=IN, raw_matrix=FILE_OUT, raw_fig=FILE_OUT,
          nrm_matrix=FILE_OUT, nrm_fig=FILE_OUT, json_matrix=FILE_OUT)
    def tb_bin(self, bamin, biases, resolution, coord1, coord2,  # pylint: disable=too-many-arguments,too-many-locals,too-many-statements,too-many-branches,no-self-use
               norm, workdir, ncpus="1", metadata=None):
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
        logger.info("TB BIN: {0}, {1}, {2}".format(bamin, resolution, workdir))

        output_metadata = {}

        if coord2 and not coord1:
            coord1, coord2 = coord2, coord1

        if not coord1:
            region1 = None
            start1 = None
            end1 = None
            region2 = None
            start2 = None
            end2 = None
        else:
            try:
                crm1, pos1 = coord1.split(':')
                start1, end1 = pos1.split('-')
                region1 = crm1
                start1 = int(start1)
                end1 = int(end1)
            except ValueError:
                region1 = coord1
                start1 = None
                end1 = None
            if coord2:
                try:
                    crm2, pos2 = coord2.split(':')
                    start2, end2 = pos2.split('-')
                    region2 = crm2
                    start2 = int(start2)
                    end2 = int(end2)
                except ValueError:
                    region2 = coord2
                    start2 = None
                    end2 = None
            else:
                region2 = None
                start2 = None
                end2 = None

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

        out_files = write_matrix(
            bamin, resolution,
            load(open(biases)) if biases else None,
            workdir, filter_exclude=[],
            normalizations=norm, half_matrix=True,
            region1=region1, start1=start1, end1=end1,
            region2=region2, start2=start2, end2=end2,
            tmpdir=workdir, append_to_tar=None, ncpus=ncpus,
            verbose=True)
        output_files = [out_files['RAW']]
        imx = load_hic_data_from_bam(bamin, resolution,
                                     biases=biases if biases else None,
                                     region=region1,
                                     ncpus=ncpus, tmpdir=workdir)
        hic_contacts_matrix_raw_fig = workdir+"/genomic_maps_raw.png"
        focus = None
        by_chrom = None
        if start1 is not None:
            focus = ((start1/resolution) if start1 != 0 else 1, (end1/resolution))
        else:
            if region1 is not None:
                focus = region1
            else:
                by_chrom = 'all'
                hic_contacts_matrix_raw_fig = workdir+"/genomic_maps_raw"
        hic_map(imx, resolution, savefig=hic_contacts_matrix_raw_fig,
                normalized=False, by_chrom=by_chrom, focus=focus)
        if by_chrom == 'all':
            hic_map(imx, resolution, savefig=hic_contacts_matrix_raw_fig+"/full_map.png",
                    normalized=False, by_chrom=None, focus=focus)
        output_files.append(hic_contacts_matrix_raw_fig)

        if len(norm) > 1:
            output_files.append(out_files['NRM'])
            hic_contacts_matrix_norm_fig = workdir+"/genomic_maps_nrm.png"
            if by_chrom == 'all':
                hic_contacts_matrix_norm_fig = workdir+"/genomic_maps_nrm"
            hic_map(imx, resolution, savefig=hic_contacts_matrix_norm_fig,
                    normalized=True, by_chrom=by_chrom, focus=focus)
            if by_chrom == 'all':
                hic_map(imx, resolution, savefig=hic_contacts_matrix_norm_fig+"/full_map.png",
                        normalized=True, by_chrom=None, focus=focus)
            output_files.append(hic_contacts_matrix_norm_fig)

        json_chr = Chromosome(name="VRE Chromosome", species=metadata["species"],
                              assembly=metadata["assembly"], max_tad_size=260000)
        json_chr.add_experiment(
            "exp1", resolution,
            hic_data=imx,
            silent=True)

        exp = json_chr.experiments[0]
        json_file_name = out_files['RAW'].replace(".abc", ".json")
        if start1 is not None:
            if focus[1]-focus[0] > 1200:
                logger.info("TADkit json too big, limiting to 1200 bins")
                focus = (focus[0], focus[0]+1200)
        else:
            if exp.size > 1200:
                logger.info("TADkit json too big, limiting to 1200 bins")
                focus = (1, 1200)
        exp.write_json(json_file_name, focus=focus)
        output_files.append(json_file_name)

        return (output_files, output_metadata)

    def run(self, input_files, input_metadata, output_files):  # pylint: disable=too-many-locals
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
        input_metadata : dict
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

        if not os.path.isfile(bamin.replace('.bam', '.bam.bai')):
            logger.info('Creating bam index')
            _cmd = ['samtools', 'index', bamin]
            out, err = Popen(_cmd, stdout=PIPE, stderr=PIPE).communicate()
            logger.info(out)
            logger.info(err)

        biases = None
        if len(input_files) > 1:
            biases = input_files[1]
        resolution = 100000
        if 'resolution' in input_metadata:
            resolution = int(input_metadata['resolution'])

        coord1 = coord2 = norm = None
        ncpus = 1
        if 'ncpus' in input_metadata:
            ncpus = input_metadata['ncpus']

        if 'coord1' in input_metadata:
            coord1 = input_metadata['coord1']
        if 'coord2' in input_metadata:
            coord2 = input_metadata['coord2']
        if 'norm' in input_metadata:
            norm = input_metadata['norm']

        root_name = os.path.dirname(os.path.abspath(bamin))
        if 'workdir' in input_metadata:
            root_name = input_metadata['workdir']

        # input and output share most metadata
        project_metadata = {}
        project_metadata["species"] = input_metadata["species"]
        project_metadata["assembly"] = input_metadata["assembly"]

        output_files, output_metadata = self.tb_bin(bamin, biases, resolution,
                                                    coord1, coord2, norm, root_name,
                                                    ncpus, project_metadata)

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
