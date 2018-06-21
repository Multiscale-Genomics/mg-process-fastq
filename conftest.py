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

Config setting for pytest
"""


def pytest_configure(config):  # pylint: disable=unused-argument
    """
    Additional settings for pytest
    """
    import sys
    sys._run_from_cmdl = True  # pylint: disable=protected-access


def pytest_unconfigure(config):  # pylint: disable=unused-argument
    """
    Remove additional settings for pytest
    """
    import sys
    del sys._run_from_cmdl
