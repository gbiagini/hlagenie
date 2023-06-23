import pytest


# TODO - add more tests to this case
# test handling of deleted alleles
def test_deleted_alleles():
    assert aa_mm.getAAposition("B*38:158", 180) == "Q"
    assert aa_mm.getAAposition("B*38:158", 178) == "T"
