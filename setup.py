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

from setuptools import setup, find_packages
# from setuptools.command.install import install

# class Install_DamIDSeq(install):
#     """
#     Custom install steps to prepare the python environment for running iDEAR
#     """

#     def run(self):
#         import rpy2.robjects as robjects
#         from rpy2.robjects.packages import importr
#         import rpy2.robjects.packages as rpackages

#         base = importr("base")
#         utils = importr("utils")
#         utils.chooseCRANmirror(ind=1)

#         # evaluate locally a remote R script
#         base.source("http://www.bioconductor.org/biocLite.R")
#         biocinstaller = importr("BiocInstaller")
#         biocinstaller.biocLite("BSgenome")
#         biocinstaller.biocLite("DESeq2")

#         utils.install_packages("devtools")
#         dt = importr("devtools")
#         dt.install_bitbucket("juanlmateo/idear")

setup(
    name='mg_process_fastq',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy', 'h5py', 'pytest'
    ],
    setup_requires=[
        'pytest-runner',
    ],
    tests_require=[
        'pytest',
    ],
    # cmdclass={
    #     'install' : Install_DamIDSeq,
    # },
)
