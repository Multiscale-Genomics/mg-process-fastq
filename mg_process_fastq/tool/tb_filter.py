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
import os.path

from utils import logger

from pytadbit.parsers.hic_bam_parser import bed2D_to_BAMhic

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    logger.info("[Warning] Cannot import \"pycompss\" API packages.")
    logger.info("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN
    from utils.dummy_pycompss import task
    from utils.dummy_pycompss import compss_wait_on

from basic_modules.tool import Tool

from pytadbit.mapping.filter import apply_filter, filter_reads
from pytadbit.mapping.analyze import insert_sizes

# ------------------------------------------------------------------------------

class tbFilterTool(Tool):
    """
    Tool for filtering out experimetnal artifacts from the aligned data
    """

    def __init__(self):
        """
        Init function
        """
        logger.info("TADbit filter aligned reads")
        Tool.__init__(self)

    @task(
        reads=FILE_IN, filter_reads_file=FILE_OUT, custom_filter=IN, conservative=IN,
        output_de=FILE_OUT, output_d=FILE_OUT, output_e=FILE_OUT,
        output_ed=FILE_OUT, output_or=FILE_OUT, output_rb=FILE_OUT,
        output_sc=FILE_OUT, output_tc=FILE_OUT, output_tl=FILE_OUT,
        output_ts=FILE_OUT, returns=int)
    def tb_filter(
            self, reads, filter_reads_file, custom_filter, min_dist_RE, min_fragment_size,
            max_fragment_size, conservative, output_de, output_d,
            output_e, output_ed, output_or, output_rb, output_sc, output_tc,
            output_tl, output_ts):
        """
        Function to filter out expoerimental artifacts

        Parameters
        ----------
        reads : str
            Location of the reads thats that has a matching location at both
            ends of the paired reads
        filtered_reads_file : str
            Location of the filtered reads
        conservative : bool
            Level of filtering [DEFAULT : True]


        Returns
        -------
        filtered_reads : str
            Location of the filtered reads

        """

        reads_tmp = reads.replace(".tsv", '')

        with open(reads_tmp + "_tmp.tsv", "wb") as f_out:
            with open(reads, "rb") as f_in:
                f_out.write(f_in.read())
        masked = filter_reads(
            reads_tmp + "_tmp.tsv",
            max_molecule_length=610,
            min_dist_to_re=min_dist_RE,
            over_represented=0.005,
            max_frag_size=max_fragment_size,
            min_frag_size=min_fragment_size,
            re_proximity=4)

        filter_reads_file_tmp = filter_reads_file.replace(".tsv", '')
        filters_suffixes = ['self-circle', 'dangling-end', 'error', 'extra dangling-end', 'too close from REs', 'too short', 'too large', 'over-represented', 'duplicated', 'random breaks']

        if custom_filter:
            applied_filters = custom_filter
            filters_suffixes = [filters_suffixes[i-1] for i in applied_filters]
        else:
            if conservative is True:
                applied_filters = [1, 2, 3, 4, 6, 7, 8, 9, 10]
                # Ignore filter 5 (based on docs) as not very helpful
            else:
                # Less conservative option
                applied_filters = [1, 2, 3, 9, 10]

        apply_filter(reads_tmp + "_tmp.tsv", filter_reads_file_tmp + "_tmp.tsv", masked, filters=applied_filters)

        with open(filter_reads_file, "wb") as f_out:
            with open(filter_reads_file_tmp + "_tmp.tsv", "rb") as f_in:
                f_out.write(f_in.read())

        for i in filters_suffixes:
            report_file_loc = reads_tmp + '_tmp.tsv_' + i + '.tsv'
            logger.info(report_file_loc)
            if os.path.isfile(report_file_loc) is True:
                logger.info("- Present {0}".format(os.path.getsize(report_file_loc)))
                with open(report_file_loc, "rb") as f_in:
                    if i == 'dangling-end':
                        logger.info("- Saving to:" + output_de)
                        with open(output_de, "wb") as f_out:
                            f_out.write(f_in.read())
                    elif i == 'duplicated':
                        logger.info("- Saving to:" + output_d)
                        with open(output_d, "wb") as f_out:
                            f_out.write(f_in.read())
                    elif i == 'error':
                        logger.info("- Saving to:" + output_e)
                        with open(output_e, "wb") as f_out:
                            f_out.write(f_in.read())
                    elif i == 'extra_dangling-end':
                        logger.info("- Saving to:" + output_ed)
                        with open(output_ed, "wb") as f_out:
                            f_out.write(f_in.read())
                    elif i == 'over-represented':
                        logger.info("- Saving to:" + output_or)
                        with open(output_or, "wb") as f_out:
                            f_out.write(f_in.read())
                    elif i == 'random_breaks':
                        logger.info("- Saving to:" + output_rb)
                        with open(output_rb, "wb") as f_out:
                            f_out.write(f_in.read())
                    elif i == 'self-circle':
                        logger.info("- Saving to:" + output_sc)
                        with open(output_sc, "wb") as f_out:
                            f_out.write(f_in.read())
                    elif i == 'too_close_from_RES':
                        logger.info("- Saving to:" + output_tc)
                        with open(output_tc, "wb") as f_out:
                            f_out.write(f_in.read())
                    elif i == 'too_large':
                        logger.info("- Saving to:" + output_tl)
                        with open(output_tl, "wb") as f_out:
                            f_out.write(f_in.read())
                    elif i == 'too_short':
                        logger.info("- Saving to:" + output_ts)
                        with open(output_ts, "wb") as f_out:
                            f_out.write(f_in.read())

        return masked

    def run(self, input_files, output_files, metadata=None):
        """
        The main function to filter the reads to remove experimental artifacts

        Parameters
        ----------
        input_files : list
            reads : str
                Location of the reads thats that has a matching location at both
                ends of the paired reads
        metadata : dict
            conservative : bool
                Level of filtering to apply [DEFAULT : True]


        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects

        """

        reads = input_files[0]

        conservative = True
        custom_filter = None
        if 'custom_filter' in metadata:
            custom_filter = metadata['filters']

        elif 'conservative' in metadata:
            conservative = metadata['conservative']

        min_dist_RE = 915
        max_fragment_size = 100000
        min_fragment_size = 100
        if 'min_dist_RE' in metadata:
            min_dist_RE = int(metadata['min_dist_RE'])
        if 'min_fragment_size' in metadata:
            min_fragment_size = int(metadata['min_fragment_size'])
        if 'max_fragment_size' in metadata:
            max_fragment_size = int(metadata['max_fragment_size'])

        root_name = reads.split("/")

        filtered_reads_file = "/".join(root_name[0:-1]) + '/' + metadata['expt_name'] + '_filtered_map.tsv'

        output_de = filtered_reads_file + '_dangling-end.tsv'
        output_d = filtered_reads_file + '_duplicated.tsv'
        output_e = filtered_reads_file + '_error.tsv'
        output_ed = filtered_reads_file + '_extra_dangling-end.tsv'
        output_or = filtered_reads_file + '_over-represented.tsv'
        output_rb = filtered_reads_file + '_random_breaks.tsv'
        output_sc = filtered_reads_file + '_self-circle.tsv'
        output_tc = filtered_reads_file + '_too_close_from_RES.tsv'
        output_tl = filtered_reads_file + '_too_large.tsv'
        output_ts = filtered_reads_file + '_too_short.tsv'

        # input and output share most metadata
        output_metadata = {}

        # handle error
        results = self.tb_filter(
            reads, filtered_reads_file, custom_filter, min_dist_RE,
            min_fragment_size, max_fragment_size, conservative,
            output_de, output_d, output_e, output_ed, output_or, output_rb,
            output_sc, output_tc, output_tl, output_ts)
        results = compss_wait_on(results)

        if 'outbam' in metadata:
            outbam = metadata['root_dir'] + '/' + metadata['outbam']
            bed2D_to_BAMhic(filtered_reads_file, True, 32, outbam, 'mid', results)
            filtered_reads_file = outbam

        hist_path = ''
        if 'histogram' in metadata:
            hist_path = "/".join(root_name[0:-1]) + '/histogram_fragment_sizes_.png'
            log_path = "/".join(root_name[0:-1]) + '/filter_log.txt'

            median, max_f, mad = insert_sizes(
                reads, nreads=1000000, stats=('median', 'first_decay', 'MAD'),
                savefig=hist_path)

            orig_stdout = sys.stdout
            f = open(log_path, "w")
            sys.stdout = f

            #insert size
            logger.info ('Insert size\n')

            logger.info ('  - median insert size = {0}'.format(median))
            logger.info ('  - double median absolution of insert size = {0}'.format(mad))
            logger.info ('  - max insert size (when a gap in continuity of > 10 bp is found in fragment lengths) = {0}'.format(max_f))

            max_mole = max_f # pseudo DEs
            min_dist = max_f + mad # random breaks
            logger.info ('   Using the maximum continuous fragment size'
                   '('+str(max_mole)+' bp) to check '
                   'for pseudo-dangling ends')
            logger.info ('   Using maximum continuous fragment size plus the MAD '
                   '('+str(min_dist)+' bp) to check for random breaks')

            sys.stdout = orig_stdout
            f.close()


        return ([filtered_reads_file, log_path, hist_path], output_metadata)

# ------------------------------------------------------------------------------
