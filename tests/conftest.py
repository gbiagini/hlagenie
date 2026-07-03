"""Shared fixtures and configuration for the HLAGenie test suite.

The tests are written as a *characterization* (golden-master) suite: they lock in
the behaviour of the current, py-ard-based implementation so that a py-ard upgrade
(or any other refactor) can be validated against a known-good baseline.

All golden values below are pinned to IMGT/HLA database version 3510. Changing
``IMGT_VERSION`` will invalidate the exact-value assertions in the integration
tests (``test_genie.py`` and ``test_data_repository.py``).
"""

import pytest
import hlagenie

# IMGT/HLA database version every golden value in this suite is pinned to.
IMGT_VERSION = "3510"


@pytest.fixture(scope="session")
def genie():
    """Session-scoped ungapped GENIE object (the default configuration)."""
    return hlagenie.init(IMGT_VERSION)


@pytest.fixture(scope="session")
def genie_gapped():
    """Session-scoped gapped GENIE object (``ungap=False``)."""
    return hlagenie.init(IMGT_VERSION, ungap=False)
