"""Tests for hlagenie.redux, the py-ard-free two-field (U2) reducer.

The pure tests need no network. ``test_parity_with_pyard`` is the accuracy
guarantee: it rebuilds the reduction data from the IMGT files and asserts the
reducer matches py-ard's ``redux(allele, "U2")`` for *every* allele ID in the
IMGT/HLA database. It is skipped if py-ard is not installed (py-ard is only a
test dependency now — it is the oracle, not a runtime requirement).
"""

import pytest

from hlagenie import redux

# Pinned to match the rest of the golden-master suite (see conftest.IMGT_VERSION).
IMGT_VERSION = "3510"


class TestFieldHelpers:
    def test_get_n_field_allele(self):
        assert redux.get_n_field_allele("A*01:01:01:01", 2) == "A*01:01"
        assert redux.get_n_field_allele("A*01:01:01:01", 3) == "A*01:01:01"

    def test_preserve_expression(self):
        assert (
            redux.get_n_field_allele("A*01:01:38L", 2, preserve_expression=True)
            == "A*01:01L"
        )
        # No expression char -> nothing appended.
        assert (
            redux.get_n_field_allele("A*01:01:01:01", 2, preserve_expression=True)
            == "A*01:01"
        )

    def test_get_2field_and_3field_strip_pg(self):
        assert redux.get_2field_allele("A*01:01:01:01G") == "A*01:01"
        assert redux.get_3field_allele("A*01:01:01G") == "A*01:01:01"


class TestReduceToTwoField:
    def test_plain_reduction(self):
        assert redux.reduce_to_two_field("A*01:01:01:01", {"A*01:01"}, {}) == "A*01:01"

    def test_already_two_field(self):
        assert redux.reduce_to_two_field("A*01:01", {"A*01:01"}, {}) == "A*01:01"

    def test_p_not_g_override_applies_to_two_field(self):
        # A P-not-G override reduces even an already-two-field allele.
        assert (
            redux.reduce_to_two_field("A*01:335", {"A*01:335"}, {"A*01:335": "A*01:01"})
            == "A*01:01"
        )

    def test_expression_suffix_preserved_when_valid(self):
        assert (
            redux.reduce_to_two_field("DPA1*03:05:01:01Q", {"DPA1*03:05Q"}, {})
            == "DPA1*03:05Q"
        )

    def test_expression_suffix_dropped_when_invalid(self):
        # The suffixed two-field is not valid -> fall back to plain two-field.
        assert (
            redux.reduce_to_two_field("B*15:01:01:02N", {"B*15:01"}, {}) == "B*15:01"
        )

    def test_non_g_group_locus_unchanged(self):
        assert redux.reduce_to_two_field("MICA*008:04:05", set(), {}) == "MICA*008:04:05"


class TestBuilders:
    def test_build_valid_alleles(self):
        names = ["A*01:01:01:01", "A*01:04:01:01N", "A*01:04:01:02N"]
        valid = redux.build_valid_alleles(names)
        assert "A*01:01:01:01" in valid       # full name
        assert "A*01:01" in valid             # 2-field
        assert "A*01:01:01" in valid          # 3-field
        # Whole A*01:04 group shares 'N' -> A*01:04N is valid.
        assert "A*01:04N" in valid

    def test_build_valid_alleles_mixed_expression_not_promoted(self):
        # A*01:07 group has both an N and an L -> no 2-field expression form.
        names = ["A*01:07:01:01N", "A*01:07:02L"]
        valid = redux.build_valid_alleles(names)
        assert "A*01:07N" not in valid
        assert "A*01:07L" not in valid

    def test_build_p_not_g(self):
        g_alleles = ["A*01:01:01:01"]                       # 2d: A*01:01
        p_pairs = [("A*01:335", "A*01:01P")]                # 2d: A*01:335 (not in G)
        p_not_g = redux.build_p_not_g(g_alleles, p_pairs)
        assert p_not_g == {"A*01:335": "A*01:01"}


@pytest.mark.network
class TestParityWithPyard:
    def test_parity_with_pyard(self):
        """The reducer must match py-ard for every allele ID in the DB."""
        pyard = pytest.importorskip("pyard")
        from hlagenie.load import (
            load_allele_names,
            load_g_group_alleles,
            load_p_group_pairs,
        )

        valid = redux.build_valid_alleles(load_allele_names(IMGT_VERSION))
        p_not_g = redux.build_p_not_g(
            load_g_group_alleles(IMGT_VERSION), load_p_group_pairs(IMGT_VERSION)
        )
        ard = pyard.init(IMGT_VERSION, load_mac=False)

        names = load_allele_names(IMGT_VERSION)
        mismatches = [
            (n, redux.reduce_to_two_field(n, valid, p_not_g), ard.redux(n, "U2"))
            for n in names
        ]
        mismatches = [t for t in mismatches if t[1] != t[2]]
        assert mismatches == [], f"{len(mismatches)} mismatches, e.g. {mismatches[:5]}"
