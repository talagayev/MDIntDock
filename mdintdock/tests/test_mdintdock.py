"""
Unit and regression test for the mdintdock package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import mdintdock


def test_mdintdock_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "mdintdock" in sys.modules
