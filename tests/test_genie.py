"""Characterization tests for the hlagenie.genie.GENIE public API.

Every golden value is pinned to IMGT/HLA version 3510 (see conftest.IMGT_VERSION).
The first run downloads MSF files and builds the SQLite cache; subsequent runs
are offline and fast. Marked ``network`` for that reason.
"""

import pytest

pytestmark = pytest.mark.network


# (allele, mature_len, first_ten, ard_len, xrd_len, aa_pos1, aa_pos9)
LOCUS_REPRESENTATIVES = [
    ("A*01:01", 341, "GSHSMRYFFT", 183, 275, "G", "F"),
    ("A*02:01", 341, "GSHSMRYFFT", 183, 275, "G", "F"),
    ("B*07:02", 338, "GSHSMRYFYT", 183, 275, "G", "Y"),
    ("C*01:02", 342, "CSHSMKYFFT", 183, 275, "C", "F"),
    ("DRB1*01:01", 237, "GDTRPRFLWQ", 95, 189, "G", "W"),
    ("DRB3*01:01", 237, "GDTRPRFLEL", 95, 189, "G", "E"),
    ("DRB4*01:01", 237, "GDTQPRFLEQ", 95, 189, "G", "E"),
    ("DRB5*01:01", 237, "GDTRPRFLQQ", 95, 189, "G", "Q"),
    ("DQA1*01:01", 232, "EDIVADHVAS", 88, 182, "E", "A"),
    ("DQB1*05:01", 229, "RDSPEDFVYQ", 95, 189, "R", "Y"),
    ("DPA1*01:03", 229, "IKADHVSTYA", 85, 179, "I", "Y"),
    ("DPB1*01:01", 229, "RATPENYVYQ", 93, 187, "R", "Y"),
]


class TestGetAA:
    def test_known_positions(self, genie):
        assert genie.getAA("A*02:01", 8) == "F"
        assert genie.getAA("A*01:01", 8) == "F"
        assert genie.getAA("A*02:01", 44) == "R"
        assert genie.getAA("A*01:01", 44) == "K"
        assert genie.getAA("A*02:01", 45) == "M"
        assert genie.getAA("A*01:01", 45) == "M"

    def test_one_indexed(self, genie):
        # Position 1 is the first residue of the mature protein.
        assert genie.getAA("A*01:01", 1) == "G"

    def test_three_field_input_is_reduced(self, genie):
        assert genie.getAA("A*01:01:01", 1) == "G"

    def test_four_field_input_is_reduced(self, genie):
        assert genie.getAA("A*02:01:01:01", 44) == "R"

    @pytest.mark.parametrize(
        "allele,_len,_first,_ard,_xrd,pos1,pos9", LOCUS_REPRESENTATIVES
    )
    def test_positions_across_loci(
        self, genie, allele, _len, _first, _ard, _xrd, pos1, pos9
    ):
        assert genie.getAA(allele, 1) == pos1
        assert genie.getAA(allele, 9) == pos9


class TestGetNuc:
    def test_two_field_only_allele(self, genie):
        # A*01:06 exists in the nucleotide alignment under its 2-field name.
        assert "".join(genie.getNuc("A*01:06", i) for i in range(1, 7)) == "ATGGCC"

    def test_known_limitation_common_allele_raises(self, genie):
        # KNOWN QUIRK: nuc_seqs is keyed by full (4-field) allele IDs, but getNuc
        # reduces multi-colon input to the 2-field name via py-ard, which is then
        # absent from nuc_seqs. Documented here so an intentional fix will flag it.
        with pytest.raises(KeyError):
            genie.getNuc("A*01:01:01:01", 1)


class TestGetPeptide:
    def test_first_ten_residues(self, genie):
        assert genie.getPeptide("A*01:01", 1, 10) == "GSHSMRYFFT"

    def test_inclusive_of_stop(self, genie):
        # start-1 : stop slice -> positions [start, stop] inclusive.
        assert genie.getPeptide("A*01:01", 1, 1) == "G"
        assert genie.getPeptide("A*01:01", 2, 3) == "SH"


class TestGetEpitope:
    def test_formatting(self, genie):
        assert genie.getEpitope("A*02:01", [1, 3, 5, 8, 10]) == "1G_3H_5M_8F_10T"

    def test_single_position(self, genie):
        assert genie.getEpitope("A*01:01", [44]) == "44K"


class TestIsPositionMismatched:
    def test_mismatch(self, genie):
        assert genie.isPositionMismatched("A*02:01", "A*01:01", 44) is True

    def test_match(self, genie):
        assert genie.isPositionMismatched("A*02:01", "A*01:01", 45) is False


class TestCountAAMismatches:
    def test_both_donor_aas_absent_from_recipient(self, genie):
        assert genie.countAAMismatches("Y", "Y", "D", "D") == 2

    def test_directional_no_mismatch_when_recipient_shares(self, genie):
        assert genie.countAAMismatches("Y", "Y", "Y", "D") == 0

    def test_recipient_has_both(self, genie):
        assert genie.countAAMismatches("Y", "D", "Y", "D") == 0
        assert genie.countAAMismatches("D", "Y", "Y", "D") == 0


class TestCountAAMismatchesAllele:
    def test_homozygous_donor_collapses_to_one(self, genie):
        # Donor homozygous with two mismatches is counted as a single mismatch.
        assert (
            genie.countAAMismatchesAllele(
                "A*02:01", "A*02:01", "A*01:01", "A*01:01", 44
            )
            == 1
        )

    def test_heterozygous_donor(self, genie):
        assert (
            genie.countAAMismatchesAllele(
                "A*01:01", "A*02:01", "A*01:01", "A*01:01", 44
            )
            == 1
        )


class TestGetARD:
    def test_length_and_bounds(self, genie):
        ard = genie.getARD("A*01:01")
        assert len(ard) == 183
        assert ard[:5] == "GSHSM"
        assert ard[-5:] == "LQRTD"

    @pytest.mark.parametrize(
        "allele,_len,_first,ard_len,_xrd,_p1,_p9", LOCUS_REPRESENTATIVES
    )
    def test_ard_length_across_loci(
        self, genie, allele, _len, _first, ard_len, _xrd, _p1, _p9
    ):
        assert len(genie.getARD(allele)) == ard_len


class TestGetXRD:
    def test_length_and_bounds(self, genie):
        xrd = genie.getXRD("A*01:01")
        assert len(xrd) == 275
        assert xrd[:5] == "GSHSM"
        assert xrd[-5:] == "TLRWE"

    @pytest.mark.parametrize(
        "allele,_len,_first,_ard,xrd_len,_p1,_p9", LOCUS_REPRESENTATIVES
    )
    def test_xrd_length_across_loci(
        self, genie, allele, _len, _first, _ard, xrd_len, _p1, _p9
    ):
        assert len(genie.getXRD(allele)) == xrd_len


class TestSequenceDictionaries:
    def test_mature_sequence_metadata(self, genie):
        for allele, mature_len, first_ten, _ard, _xrd, _p1, _p9 in LOCUS_REPRESENTATIVES:
            seq = genie.seqs[allele]
            assert len(seq) == mature_len
            assert seq[:10] == first_ten

    def test_dictionary_sizes(self, genie):
        assert len(genie.seqs) == 21875
        assert len(genie.full_seqs) == 21875
        assert len(genie.nuc_seqs) == 34714

    def test_full_seq_contains_leader(self, genie):
        # The full sequence includes the leader peptide, so it is longer than
        # the mature sequence and does not start at the mature motif.
        assert len(genie.full_seqs["A*01:01"]) > len(genie.seqs["A*01:01"])


class TestListFunctions:
    # (locus, completes, incompletes, extendeds) for protein sequences, v3510.
    LOCUS_COUNTS = [
        ("A", 2164, 2248, 7),
        ("B", 2606, 2840, 15),
        ("C", 2258, 1918, 58),
        ("DRB1", 244, 1983, 4),
        ("DRB3", 17, 320, 1),
        ("DRB4", 5, 140, 0),
        ("DRB5", 6, 137, 0),
        ("DQA1", 67, 195, 0),
        ("DQB1", 348, 1123, 36),
        ("DPA1", 94, 142, 2),
        ("DPB1", 317, 1007, 4),
    ]

    @pytest.mark.parametrize("locus,completes,incompletes,extendeds", LOCUS_COUNTS)
    def test_counts(self, genie, locus, completes, incompletes, extendeds):
        assert len(genie.listCompletes(locus)) == completes
        assert len(genie.listIncompletes(locus)) == incompletes
        assert len(genie.listExtendeds(locus)) == extendeds

    def test_reference_allele_is_complete(self, genie):
        assert "A*01:01" in genie.listCompletes("A")

    def test_invalid_seqtype_returns_none(self, genie):
        assert genie.listCompletes("A", seqtype="bogus") is None
        assert genie.listIncompletes("A", seqtype="bogus") is None
        assert genie.listExtendeds("A", seqtype="bogus") is None


class TestGappedConfiguration:
    def test_gapped_domain_ends_differ_from_ungapped(self, genie_gapped):
        # Gapped domain boundaries include alignment gaps, so they are longer.
        assert genie_gapped.ards["A"] == 202
        assert genie_gapped.xrds["A"] == 295

    def test_gapped_ard_is_at_least_ungapped(self, genie, genie_gapped):
        for locus in genie.ards:
            assert genie_gapped.ards[locus] >= genie.ards[locus]
