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

import os.path
import shlex
import shutil
import subprocess
import tarfile

from utils import logger


class cd(object):  # pylint: disable=too-few-public-methods, invalid-name
    """
    Context manager for changing the current working directory
    """

    def __init__(self, newpath):
        self.savedpath = os.getcwd()
        self.newpath = os.path.expanduser(newpath)

    def __enter__(self):
        self.savedpath = os.getcwd()
        os.chdir(self.newpath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedpath)


class common(object):  # pylint: disable=too-few-public-methods, invalid-name
    """
    Common functions that can be used generically across tools and pipelines
    """

    @staticmethod
    def tar_folder(folder, tar_file, archive_name="tmp", keep_folder=False):
        """
        Archive a folder as a tar file.

        This function will overwrite an existing file of the same name.

        The function adds each of the files to the archive individually and
        removing the original once it has been added to save space on in the
        storage area

        Parameters
        ----------
        folder : str
            Location of the folder that is to be archived.
        tar_file : str
            Location of the archive tar file
        archive_name : str
            Name of the dir that the files in the archive should be stored in.
            Default: tmp
        keep_folder : bool
            By default (False) the files and folder are deleted once they are added to
            the archive.
        """
        onlyfiles = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]

        tar = tarfile.open(tar_file, "w")
        for tmp_file in onlyfiles:
            tar.add(
                os.path.join(folder, tmp_file),
                arcname=os.path.join(archive_name, tmp_file)
            )
            if keep_folder is False:
                os.remove(os.path.join(folder, tmp_file))
        tar.close()

        if keep_folder is False:
            shutil.rmtree(folder)

    @staticmethod
    def zip_file(location, cpu=1):
        """
        Use pigz (gzip as a fallback) to compress a file

        Parameters
        ----------
        location : str
            Location of the file to be zipped
        """
        try:
            command_line = "pigz -p {} {}".format(cpu, location)
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()
        except OSError:
            logger.warn("OSERROR: pigz not installed, using gzip")
            command_line = 'gzip ' + location
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()
