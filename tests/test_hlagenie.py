import pytest
import hlagenie

# create an object for testing
aa_mm = hlagenie.init("3510")

# set known alleles for testing
allele1 = "A*02:01"
allele2 = "A*01:01"

# set known position for testing
position1 = 8


# test getAA function
def test_getAA():
    assert aa_mm.getAA(allele1, position1) == "F"
    assert aa_mm.getAA(allele2, position1) == "F"
    assert aa_mm.getAA(allele1, 44) == "R"
    assert aa_mm.getAA(allele2, 44) == "K"
    assert aa_mm.getAA(allele1, 45) == "M"
    assert aa_mm.getAA(allele2, 45) == "M"


# test getPeptide function
def test_getPeptide():
    # check first 10 residues
    assert aa_mm.getPeptide(allele1, 1, 10) == "GSHSMRYFFT"


# test getEpitope function
def test_getEpitope():
    assert aa_mm.getEpitope(allele1, [1, 3, 5, 8, 10]) == "1G_3H_5M_8F_10T"


# test isPositionMismatched function
def test_isPositionMismatched():
    assert aa_mm.isPositionMismatched(allele1, allele2, 44) == True
    assert aa_mm.isPositionMismatched(allele1, allele2, 45) == False


# test countAAMismatchesAllele function
def test_countAAMismatchesAllele():
    assert aa_mm.countAAMismatchesAllele(allele1, allele1, allele2, allele2, 44) == 1


# test countAAMismatches function
def test_countAAMismatches():
    # this is directional, so if the recipient has either AA from the donor, it's not a mismatch
    assert aa_mm.countAAMismatches("Y", "Y", "D", "D") == 2
    assert aa_mm.countAAMismatches("Y", "Y", "Y", "D") == 0


# test defining of ARD
def test_getARD():
    assert len(aa_mm.getARD("A*01:01")) == 183
    assert aa_mm.getARD("A*01:01")[0:5] == "GSHSM"
    assert aa_mm.getARD("A*01:01")[-5:] == "LQRTD"


# test defining of XRD
def test_getXRD():
    assert len(aa_mm.getXRD("A*01:01")) == 275
    assert aa_mm.getXRD("A*01:01")[0:5] == "GSHSM"
    assert aa_mm.getXRD("A*01:01")[-5:] == "TLRWE"


# test handling of deleted alleles
# def test_deleted_alleles():
#     assert aa_mm.getAA("B*38:158", 180) == "Q"
#     assert aa_mm.getAA("B*38:158", 178) == "T"
