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
    def to_output_file(input_file, output_file, empty=True):
        """
        When handling the output of files within the @task function copying the
        results into the correct output files should be done by reading from and
        writing to rather than renaming.

        In cases where there are a known set of output files, if the input file
        is missing then a blank file should be created and handled by the run()
        function of the tool. If an empty file should not be created then the
        empty parameter should be set to False.

        Parameters
        ----------
        input_file : str
            Location of the input file
        output_file : str
            Location of the output file
        empty : bool
            In cases where the input_file is missing an empty output_file is
            created. Should be set to False if no file shold be created.
        """
        logger.info(input_file + ' - ' + str(os.path.isfile(input_file)))
        if os.path.isfile(input_file) is True and os.path.getsize(input_file) > 0:
            with open(output_file, "wb") as f_out:
                with open(input_file, "rb") as f_in:
                    f_out.write(f_in.read())
            return True

        if empty:
            with open(output_file, "w") as f_out:
                f_out.write("")

        return True

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
