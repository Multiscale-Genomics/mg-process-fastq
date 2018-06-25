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
import subprocess

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
    def zip_file(location):
        """
        Use pigz (gzip as a fallback) to compress a file

        Parameters
        ----------
        location : str
            Location of the file to be zipped
        """
        try:
            command_line = 'pigz -p 2 ' + location
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()
        except OSError:
            logger.warn("OSERROR: pigz not installed, using gzip")
            command_line = 'gzip ' + location
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()
