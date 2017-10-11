"""
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
