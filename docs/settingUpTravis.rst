Continuous Integration with Travis
==================================

This document summarizes the steps involved in setting up a continuous integration suite with Travis for our pipelines. 

Getting Started
----------------

Login with your gitHub account details on travis.ci. Add the "Multiscale Genomics" to your organization to gain access to its repositories. (you would have to send a request to be added to this organization). 

Follow the onscreen instructions on travis.ci.

.. code-block:: none
   :linenos:

   Flick the repository switch to "on".
   Add .travis.yml to your repository.
   Sync your travis account to your GitHub account.
   Pushing to Git will trigger the first Travis build after adding the above file.
   
  
Making .travis.yml File
-----------------------

To your added .travis.yml file in your GitHub repository include :

   - "python" in "language". With version/s specified in "python: " 
   - All packages required for running the pipelines in "addons: apt: packages: ". (For more information on the packages please see Full Installation from : http://multiscale-genomics.readthedocs.io/projects/mg-process-fastq/en/latest/full_installation.html#setup-the-system-environment)
   - Docker in "services" to tell Travis it needs to have docker installed.

.. code-block:: none

   services:
       - docker
   
   
   - All tools to be installed within "install:" section (For more information on the packages please see Full Installation from : http://multiscale-genomics.readthedocs.io/projects/mg-process-fastq/en/latest/full_installation.html#setup-the-system-environment )
   .. note:: libtbb did not seem to be installing correctly when put in "apt: packages: ". It has therefore been done with sudo in "install:"
   
   - Setup all symlinks in "before_script: "
   - Change execution permissions on shims folder and harness.sh
   - Add the shims folder path to your $PATH
   - Include harness.sh to your "script" to be executed
   
   

Making harness.sh File
-----------------------
   
Include all test scripts to be tested with pytest to your harness.sh file. As iNPS does not work with Python < 3, a conditional check is present to ensure that it does not run unless Travis is running Python 3
   
   
Running Docker container
-------------------------

To run docker within Travis, it has been included within the "services" in the .yml file. Every time travis runs docker it will pull the latest publicly available image from DockerHub and run it. This gets done from within the shims files. The pre-built container can also be pulled from : https://hub.docker.com/r/multiscalegenomics/mgprocessfastq/

using command : 

.. code-block:: none

   $ docker pull multiscalegenomics/mgprocessfastq:testdocker
   
For more details on the docker container, please refer to : https://github.com/Multiscale-Genomics/mg-process-fastq/blob/master/docs/docker.rst


Setting up Shims 
-----------------   

Libmaus2 and Biobambam2 have had to be installed within a docker container, as they were causing Travis to time out. As the container is non-interactable while on Travis and is not hosting live server. An intermediate layer to access the contents within docker from travis has been introduced in the form of shims. To make these files, take the list of all biobambam2 modules and construct bash script files with the following command for each of the modules : 

.. code-block:: none

   exec docker run -it  multiscalegenomics/mgprocessfastq:biobambamimage ~/lib/biobambam2/bin/biobambam_module_name $@ 
   
   
      
The .travis.yml file used for testing the mg-process-fastq pipelines can be found at : https://github.com/Multiscale-Genomics/mg-process-fastq/blob/master/.travis.yml