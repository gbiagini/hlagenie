import pytest
import hlagenie

# fixtures for test functions


@pytest.fixture(autouse=True)
def aa_mm():
    # create an object for testing
    aa_mm = hlagenie.init("3510")

    # set known alleles for testing
    allele1 = "A*02:01"
    allele2 = "A*01:01"

    # set known position for testing
    position1 = 8

    # run tests here
    yield
