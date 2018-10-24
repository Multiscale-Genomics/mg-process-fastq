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

Workflows
=========

Workflows are the scripts that define the pipelines and launc each of the tools

.. automodule:: mg_process_fastq.workflow

   Multitool Workflows
   ===================

   ChIP-seq
   --------
   Used for the aalysis of FastQ files. Handles paired or single end FastQ inputs.
   It uses the BWA MEM alginer, filters out repeats and technical artifacts using
   BioBamBam2 and the peak calling uses MACS2.

   Methods
   ^^^^^^^

      .. autoclass:: mg_process_fastq.workflow.chipseq
         :members:

   Single Tool Workflows
   =====================