"""
.. Copyright 2017 EMBL-European Bioinformatics Institute

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
# Indexers
from . import bowtie_indexer
#import bs_seeker_indexer
from . import bwa_indexer
from . import gem_indexer
from . import kallisto_indexer

# Aligners
#import bs_seeker_aligner
from . import bwa_aligner

# Filters
from . import biobambam_filter
#import bs_seeker_filter

# Analysis
#import bs_seeker_methylation_caller
from . import inps
from . import kallisto_quant
from . import macs2

__author__  = 'Mark McDowall'
__version__ = '0.0'
__license__ = 'Apache 2.0'
