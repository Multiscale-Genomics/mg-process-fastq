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

Testing
=======

To ensure that the code written to process the data within the MuG VRE works it is essential to have a testing framework to ensure that as the code evolves it does not break the funcitonality of the pipelines and tools. There are 2 key parts to this:

1. Sample data
2. Runnable scripts

For each tool within the mg-process-fastq repository there is an initial set of sample data and matching tests for each tool and pipeline.

Pipelines
---------
There is a test for each of the tools. This uses the "process" scripts to run each of the tools. This is to ensure that the pipeline scripts are able to call out to each of the tools and the correct parameters are handed to each one.

.. automodule:: tests

   ChIP-Seq
   ========
   To run the pipeline test:

   .. code-block:: none

      pytest tests/test_pipeline_chipseq.py


   Methods
   -------
   .. automodule:: tests.test_pipeline_chipseq
      :members:

   Sample Data
   -----------


   Genome Indexing
   ===============
   To run the pipeline test:

   .. code-block:: none

      pytest tests/test_pipeline_genome.py


   Methods
   -------
   .. automodule:: tests.test_pipeline_genome
      :members:

   Sample Data
   -----------


   Hi-C
   ====
   To run the pipeline test:

   .. code-block:: none

      pytest tests/test_pipeline_tb.py


   Methods
   -------
   .. automodule:: tests.test_pipeline_tb
      :members:

   Sample Data
   -----------
   :doc:`tests_hic`


   MNase-Seq
   =========
   To run the pipeline test:

   .. code-block:: none

      pytest tests/test_pipeline_mnaseseq.py


   Methods
   -------
   .. automodule:: tests.test_pipeline_mnaseseq
      :members:

   Sample Data
   -----------


   RNA-Seq
   =======
   To run the pipeline test:

   .. code-block:: none

      pytest tests/test_pipeline_rnaseq.py


   Methods
   -------
   .. automodule:: tests.test_pipeline_rnaseq
      :members:

   Sample Data
   -----------


   Whole Genome Bisulfate Sequencing (WGBS)
   ========================================
   To run the pipeline test:

   .. code-block:: none

      pytest tests/test_pipeline_wgbs.py


   Methods
   -------
   .. automodule:: tests.test_pipeline_wgbs
      :members:

   Sample Data
   -----------


Tools
-----
As the data stored is only the raw data, each of the sets of tools has been packaged up into a tool chain to run each of the tools without failing. This has been done with the tests.test_toolchains.py script.

.. code-block:: none

   python tests/test_toolchains.py --pipeline [genome | chipseq | hic | mnaseseq | rnaseq | wgbs]

This script automates the running of each of the tools that are required for a given pipeline.


Methods
^^^^^^^
.. automodule:: tests

   .. automodule:: tests.test_toolchains
      :members: