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

Setting up and using a Docker Container
=======================================


Our reason for using a container
--------------------------------

While working with Travis, the installation of libmaus2 and biobambam2 took up more than 45 minutes (accumulative with the rest of the installations), which caused Travis to time out. We therefore resorted to putting both tools in a container and accessing the commands from there. This document summarizes the steps involved in making a docker image for the above two tools and running a container from that image. As well as uploading your image to docker hub to make it publicly accessible.

This document has been prepared keeping macOS Sierra in mind, although many of the commands are cross platform (\*nix) compliant.


Getting Started
---------------

To be able to build a docker container you must have :

a) Docker installed on your machine
b) An account on one of the docker repositories (Docker Hub or Quay). We have used Docker Hub as this was free access.

a) Installing docker to your machine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For this work I had installed a command line based docker, along with the Virtual machine boot2docker. There is however a GUI distribution available for MAC as well. You may install boot2docker using :

.. code-block:: none
   :linenos:

   brew install boot2docker

b) Setting up account on Docker Hub
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Go to https://hub.docker.com and setup an account with Docker Hub. Add Multiscale Genomics to your organizations and create a repository : mgprocessfastq. You will be uploading your docker images to this repository later.

Constructing a docker container
-------------------------------

Run the following preliminary commands to get your boot2docker running:

.. code-block:: none
   :linenos:

   boot2docker up
   eval "$(boot2docker shellinit)"


You would also need to have :

   - docker-machine
   - Virtual Box


installed on your Mac. After these, execute the following commands:

.. code-block:: none
   :linenos:

   docker-machine create -d virtualbox dev
   eval $(docker-machine env dev)


To ensure your docker is running:

.. code-block:: none
   :linenos:

   docker

Making the Dockerfile for libmaus2 and biobambam2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You may want to create a new folder for this purpose, as the docker command
compiles the Dockerfile with the given path to the folder. Create a new file
with the name of "Dockerfile". Include the following lines within this file:

.. code-block:: none

   FROM ubuntu:14.04

   RUN apt-get update && \
         apt-get -y install sudo

   RUN  sudo apt-get install -y make build-essential libssl-dev zlib1g-dev \
   libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev \
   libncursesw5-dev xz-utils tk-dev unzip mcl libgtk2.0-dev r-base-core     \
   libcurl4-gnutls-dev python-rpy2 git

   RUN mkdir Mug  \
    && cd Mug  \
    && apt-get -y install git \
    && git config --global user.name "your_username"	\
    && git config --global user.email "your_emailId"	\
    && pwd  	\
    && mkdir bin lib code 	\
    && cd lib	\
    && git clone https://github.com/gt1/libmaus2.git
    && cd libmaus2  \
    && sudo apt-get -y install libtool m4 automake \
    && libtoolize \
    && aclocal 	\
    && autoheader 	\
    && automake --force-missing --add-missing 	\
    && autoconf \
    && ./configure --prefix=/Mug/lib/libmaus2 	\

    && make  \
    && make install \
    && cd /Mug/lib 	\


    && git clone https://github.com/gt1/biobambam2.git 	&& cd biobambam2 	\
    && autoreconf -i -f	\
    && ./configure --with-libmaus2=/Mug/lib/libmaus2 --prefix=/Mug/lib/biobambam2	\
    && make install

Making the docker image
^^^^^^^^^^^^^^^^^^^^^^^

Build a docker image from this file using:

.. code-block:: none
   :linenos:

   cd /path/to/your/dockerfile
   docker build –t multiscalegenomics/mgprocessfastq/biobambamimage.

Login with your docker hub account details :

.. code-block:: none
   :linenos:

   docker login

Push the above image to your docker hub repository

.. code-block:: none
   :linenos:

   docker push multiscalegenomics/mgprocessfastq:biobambamimage


Running a docker container
^^^^^^^^^^^^^^^^^^^^^^^^^^

You should be able to run the above image locally on your machine as well as pulling it elsewhere (on a system which has docker):

.. code-block:: none
   :linenos:

   docker pull multiscalegenomics/mgprocessfastq:biobambamimage

and then running a container via :

.. code-block:: none
   :linenos:

   docker run --name name_you_want multiscalegenomics/mgprocessfastq:biobambamimage


Our Travis build pulls the image from our mgprocessfastq repository from within the shims files, and runs the containers using the commands within.




