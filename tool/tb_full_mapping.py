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

from os import path, unlink

from pytadbit.mapping.mapper import full_mapping
from pytadbit.utils.fastq_utils import quality_plot

from basic_modules.tool import Tool

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    # from pycompss.api.constraint import constraint
    from pycompss.api.api import compss_wait_on
except ImportError:
    logger.info("[Warning] Cannot import \"pycompss\" API packages.")
    logger.info("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task  # pylint: disable=ungrouped-imports
    # from utils.dummy_pycompss import constraint # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import compss_wait_on  # pylint: disable=ungrouped-imports


# ------------------------------------------------------------------------------

class tbFullMappingTool(Tool):  # pylint: disable=invalid-name
    """
    Tool for mapping fastq paired end files to the GEM index files
    """

    def __init__(self):
        """
        Init function
        """
        logger.info("TADbit full_mapping")
        Tool.__init__(self)

    @task(
        gem_file=FILE_IN, fastq_file=FILE_IN, windows=IN, window1=FILE_OUT,
        window2=FILE_OUT, window3=FILE_OUT, window4=FILE_OUT)
    def tb_full_mapping_iter(  # pylint: disable=too-many-locals,too-many-statements,unused-argument,no-self-use,too-many-arguments
            self, gem_file, fastq_file, windows,
            window1, window2, window3, window4, ncpus=1, workdir='/tmp/'):
        """
        Function to map the FASTQ files to the GEM file over different window
        sizes ready for alignment

        Parameters
        ----------
        gem_file : str
            Location of the genome GEM index file
        fastq_file_bgd : str
            Location of the FASTQ file
        windows : list
            List of lists with the window sizes to be computed
        window1 : str
            Location of the first window index file
        window2 : str
            Location of the second window index file
        window3 : str
            Location of the third window index file
        window4 : str
            Location of the fourth window index file


        Returns
        -------
        window1 : str
            Location of the first window index file
        window2 : str
            Location of the second window index file
        window3 : str
            Location of the third window index file
        window4 : str
            Location of the fourth window index file

        """
        logger.info("tb_full_mapping_iter")
        output_dir = workdir

        _ = full_mapping(
            gem_file, fastq_file, output_dir,
            windows=windows, frag_map=False, nthreads=ncpus, clean=True,
            temp_dir=workdir
        )

        return True

    @task(
        gem_file=FILE_IN, fastq_file=FILE_IN, enzyme_name=IN, windows=IN,
        full_file=FILE_OUT, frag_file=FILE_OUT)
    def tb_full_mapping_frag(  # pylint: disable=too-many-locals,too-many-statements,unused-argument,no-self-use,too-many-arguments
            self, gem_file, fastq_file, enzyme_name, windows,
            full_file, frag_file, ncpus=1, workdir='/tmp/'):
        """
        Function to map the FASTQ files to the GEM file based on fragments
        derived from the restriction enzyme that was used.

        Parameters
        ----------
        gem_file : str
            Location of the genome GEM index file
        fastq_file_bgd : str
            Location of the FASTQ file
        enzyme_name : str
            Restriction enzyme name (MboI)
        windows : list
            List of lists with the window sizes to be computed
        window_file : str
            Location of the first window index file


        Returns
        -------
        window_file : str
            Location of the window index file

        """
        logger.info("tb_full_mapping_frag")
        # od_loc = fastq_file.split("/")
        # output_dir = "/".join(od_loc[0:-1])
        output_dir = workdir

        logger.info("TB MAPPING - output_dir:", output_dir)
        logger.info("TB MAPPING - full_file dir:", full_file)
        logger.info("TB MAPPING - frag_file dir:", frag_file)
        gzipped = ''
        dsrc = ''
        if fastq_file.endswith('.fastq.gz') or fastq_file.endswith('.fq.gz'):
            gzipped = '.gz'
        if fastq_file.endswith('.fastq.dsrc'):
            dsrc = '.dsrc'

        file_name = path.basename(fastq_file)
        file_name = file_name.replace('.fastq'+dsrc+gzipped, '')
        file_name = file_name.replace('.fq'+dsrc+gzipped, '')
        fastq_file_tmp = workdir+'/'+file_name

        _ = full_mapping(
            gem_file, fastq_file, output_dir,
            r_enz=enzyme_name, windows=windows, frag_map=True, nthreads=ncpus,
            clean=True, temp_dir=workdir
        )

        with open(full_file, "wb") as f_out:
            with open(fastq_file_tmp + dsrc + "_full_1-end.map", "rb") as f_in:
                f_out.write(f_in.read())

        if path.isfile(fastq_file_tmp + dsrc + "_full_1-end.map"):
            unlink(fastq_file_tmp + dsrc + "_full_1-end.map")

        with open(frag_file, "wb") as f_out:
            with open(fastq_file_tmp + dsrc + "_frag_1-end.map", "rb") as f_in:
                f_out.write(f_in.read())

        if path.isfile(fastq_file_tmp + dsrc + "_frag_1-end.map"):
            unlink(fastq_file_tmp + dsrc + "_frag_1-end.map")

        return True

    def run(self, input_files, input_metadata, output_files):  # pylint: disable=too-many-locals,too-many-branches,too-many-statements,no-self-use,unused-argument
        """
        The main function to map the FASTQ files to the GEM file over different
        window sizes ready for alignment

        Parameters
        ----------
        input_files : list
            gem_file : str
                Location of the genome GEM index file
            fastq_file_bgd : str
                Location of the FASTQ file
        metadata : dict
            windows : list
                List of lists with the window sizes to be computed
            enzyme_name : str
                Restriction enzyme used [OPTIONAL]


        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects

        """

        gem_file = input_files[0]
        fastq_file = input_files[1]
        windows = input_metadata['windows']
        if not windows:
            windows = None
        if 'ncpus' not in input_metadata:
            input_metadata['ncpus'] = 8
        if 'iterative_mapping' in input_metadata:
            if isinstance(input_metadata['iterative_mapping'], basestring):
                frag_base = not input_metadata['iterative_mapping'].lower() in ("yes", "true", "t", "1")
            else:
                frag_base = not input_metadata['iterative_mapping']
        else:
            frag_base = (windows is None)
        if 'workdir' in input_metadata:
            root_path = input_metadata['workdir']
        else:
            root_path = path.dirname(path.abspath(fastq_file))

        gzipped = ''
        if fastq_file.lower().endswith('.fastq.gz') or fastq_file.lower().endswith('.fq.gz'):
            gzipped = '.gz'
        file_name = path.basename(fastq_file)
        file_name = (file_name.replace('.fastq'+gzipped, '')).replace('.FASTQ'+gzipped, '')
        file_name = (file_name.replace('.fq'+gzipped, '')).replace('.FQ'+gzipped, '')
        quality_plot_file = ''
        log_path = ''
        # name = root_name[-1]

        if 'quality_plot' in input_metadata:
            quality_plot_file = 'QC-plot_'+file_name + '.png'
            log_path = root_path+'/'+'mapping_log_'+file_name + '.txt'

            dangling_ends, ligated = quality_plot(
                fastq_file, r_enz=input_metadata['enzyme_name'],
                nreads=100000, paired=False,
                savefig=path.join(
                    root_path,
                    quality_plot_file))

            orig_stdout = sys.stdout
            f_handler = open(log_path, "w")
            sys.stdout = f_handler
            logger.info('Hi-C QC plot')
            for renz in dangling_ends:
                logger.info('  - Dangling-ends (sensu-stricto): ', dangling_ends[renz])
            for renz in ligated:
                logger.info('  - Ligation sites: ', ligated[renz])

            sys.stdout = orig_stdout
            f_handler.close()

        file_name = root_path+'/'+file_name
        output_metadata = {}
        # input and output share most metadata

        if frag_base:
            full_file = file_name + "_full.map"
            frag_file = file_name + "_frag.map"

            results = self.tb_full_mapping_frag(
                gem_file, fastq_file, input_metadata['enzyme_name'], None,
                full_file, frag_file, ncpus=input_metadata['ncpus'], workdir=root_path
            )
            results = compss_wait_on(results)

            output_metadata['func'] = 'frag'
            return_files = [full_file, frag_file]
            if 'quality_plot' in input_metadata:
                return_files += [log_path, root_path+'/'+quality_plot_file]
            return (return_files, output_metadata)

        window2 = window3 = window4 = None
        window1 = file_name + "_full_" + str(windows[0][0]) + "-" + str(windows[0][1]) + ".map"
        if len(windows) > 1:
            window2 = file_name + "_full_" + str(windows[1][0]) + "-" + str(windows[1][1]) + ".map"
        if len(windows) > 2:
            window3 = file_name + "_full_" + str(windows[2][0]) + "-" + str(windows[2][1]) + ".map"
        if len(windows) > 3:
            window4 = file_name + "_full_" + str(windows[3][0]) + "-" + str(windows[3][1]) + ".map"

        results = self.tb_full_mapping_iter(
            gem_file, fastq_file, windows,
            window1, window2, window3, window4, ncpus=input_metadata['ncpus'], workdir=root_path
        )
        results = compss_wait_on(results)

        output_metadata['func'] = 'iter'
        return_files = [window1, window2, window3, window4]
        if 'quality_plot' in input_metadata:
            return_files += [log_path, root_path+'/'+quality_plot_file]
        return (return_files, output_metadata)

# ------------------------------------------------------------------------------
