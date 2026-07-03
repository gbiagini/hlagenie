"""Characterization tests for hlagenie.data_repository table generation.

These exercise the ARD/XRD boundary tables and the sequence dictionaries that
``GENIE.__init__`` builds. Pinned to IMGT/HLA version 3510 and marked ``network``.
"""

import pytest

from hlagenie import data_repository as dr

pytestmark = pytest.mark.network


UNGAPPED_ARD_ENDS = {
    "A": 183, "B": 183, "C": 183, "DRB1": 95, "DRB3": 95, "DRB4": 95,
    "DRB5": 95, "DQA1": 88, "DQB1": 95, "DPA1": 85, "DPB1": 93,
}
UNGAPPED_XRD_ENDS = {
    "A": 275, "B": 275, "C": 275, "DRB1": 189, "DRB3": 189, "DRB4": 189,
    "DRB5": 189, "DQA1": 182, "DQB1": 189, "DPA1": 179, "DPB1": 187,
}
GAPPED_ARD_ENDS = {
    "A": 202, "B": 211, "C": 211, "DRB1": 98, "DRB3": 98, "DRB4": 98,
    "DRB5": 98, "DQA1": 88, "DQB1": 95, "DPA1": 86, "DPB1": 101,
}
GAPPED_XRD_ENDS = {
    "A": 295, "B": 303, "C": 303, "DRB1": 192, "DRB3": 192, "DRB4": 192,
    "DRB5": 192, "DQA1": 182, "DQB1": 192, "DPA1": 180, "DPB1": 196,
}


class TestUngappedBoundaryTables:
    def test_ard_ends(self, genie):
        assert genie.ards == UNGAPPED_ARD_ENDS

    def test_xrd_ends(self, genie):
        assert genie.xrds == UNGAPPED_XRD_ENDS

    def test_ard_within_xrd(self, genie):
        for locus in genie.ards:
            assert genie.ards[locus] < genie.xrds[locus]


class TestGappedBoundaryTables:
    def test_ard_ends(self, genie_gapped):
        assert genie_gapped.ards == GAPPED_ARD_ENDS

    def test_xrd_ends(self, genie_gapped):
        assert genie_gapped.xrds == GAPPED_XRD_ENDS


class TestSequenceTableConsistency:
    def test_mature_is_prefix_slice_of_full(self, genie):
        # Mature sequences are produced by trimming the leader from the full
        # sequence, so every mature sequence is a suffix of its full sequence.
        for allele in ["A*01:01", "B*07:02", "DRB1*01:01"]:
            assert genie.full_seqs[allele].endswith(genie.seqs[allele])

    def test_every_mature_allele_has_full_sequence(self, genie):
        assert set(genie.seqs) == set(genie.full_seqs)


class TestSetDbVersion:
    def test_returns_existing_version(self, genie):
        # The version was stamped during init; re-stamping is a no-op that
        # returns the stored value.
        assert dr.set_db_version(genie.db_connection, "3510") == 3510
