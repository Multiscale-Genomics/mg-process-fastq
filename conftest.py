"""
Config setting for pytest
"""

def pytest_configure(config):
    """
    Additional settings for pytest
    """
    import sys
    sys._run_from_cmdl = True

def pytest_unconfigure(config):
    """
    Remove additional settings for pytest
    """
    import sys
    del sys._run_from_cmdl
