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

from basic_modules.workflow import Workflow
from utils import logger

from mg_process_fastq.tool.biobambam_filter import biobambam


# ------------------------------------------------------------------------------

class process_biobambam(Workflow):  # pylint disable=too-few-public-methods, invalid-name
    """
    Functions for filtering FastQ alignments with BioBamBam.
    """

    def __init__(self, configuration=None):
        """
        Initialise the class

        Parameters
        ----------
        configuration : dict
            a dictionary containing parameters that define how the operation
            should be carried out, which are specific to each Tool.
        """
        logger.info("Processing BioBamBam Filtering")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        Main run function for filtering FastQ aligned reads using BioBamBam.

        Parameters
        ----------
        input_files : dict
            Location of the initial input files required by the workflow

            bam : str
                Location of BAM file

        metadata : dict
            Input file meta data associated with their roles

            bam : str

        output_files : dict
            Output file locations

            filtered : str

        Returns
        -------
        output_files : dict
            Output file locations associated with their roles, for the output

            filtered : str
                Filtered version of the bam file

        output_metadata : dict
            Output metadata for the associated files in output_files

            filtered : Metadata
        """
        output_files_generated = {}
        output_metadata = {}

        # Filter the bam
        b3f = biobambam(self.configuration)

        logger.progress("BioBamBam Filter", status="RUNNING")
        b3f_files, b3f_meta = b3f.run(
            {"input": input_files['bam']},
            {"input": metadata['bam']},
            {"output": output_files["filtered"]}
        )
        logger.progress("BioBamBam Filter", status="DONE")

        try:
            output_files_generated["filtered"] = b3f_files["bam"]
            output_metadata["filtered"] = b3f_meta["bam"]

            tool_name = output_metadata['filtered'].meta_data['tool']
            output_metadata['filtered'].meta_data['tool_description'] = tool_name
            output_metadata['filtered'].meta_data['tool'] = "process_biobambam"

            output_files_generated["filtered_bai"] = b3f_files["bai"]
            output_metadata["filtered_bai"] = b3f_meta["bai"]

            tool_name = output_metadata['filtered_bai'].meta_data['tool']
            output_metadata['filtered_bai'].meta_data['tool_description'] = tool_name
            output_metadata['filtered_bai'].meta_data['tool'] = "process_biobambam"
        except KeyError:
            logger.fatal("BioBamBam filtering failed")

        return output_files_generated, output_metadata
