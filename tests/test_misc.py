"""Tests for hlagenie.misc helper functions.

All tests here are pure (no network) except ``test_get_imgt_db_versions``,
which is marked ``network``.
"""

import pathlib
import tempfile

import pytest

from hlagenie.misc import (
    coordinate,
    coordinate_end,
    find_gaps,
    get_data_dir,
    get_default_db_directory,
    get_imgt_db_versions,
    get_imgt_version,
    regex_gen,
)


class TestGetImgtVersion:
    def test_dotted_version_is_normalised(self):
        assert get_imgt_version("3.51.0") == "3510"

    def test_plain_digits_pass_through(self):
        assert get_imgt_version("3510") == "3510"


class TestGetDataDir:
    def test_none_returns_default(self):
        assert get_data_dir(None) == get_default_db_directory()

    def test_valid_directory_is_returned_as_path(self):
        with tempfile.TemporaryDirectory() as td:
            assert get_data_dir(td) == pathlib.Path(td)

    def test_invalid_directory_raises(self):
        with pytest.raises(RuntimeError):
            get_data_dir("/no/such/dir/xyz123")


class TestGetDefaultDbDirectory:
    def test_ends_with_hlagenie(self):
        d = get_default_db_directory()
        assert isinstance(d, pathlib.Path)
        assert d.name == "hlagenie"


class TestFindGaps:
    def test_returns_indices_of_dash_characters(self):
        assert find_gaps("AB-CD--E") == [2, 5, 6]

    def test_no_gaps_returns_empty_list(self):
        assert find_gaps("ABCDE") == []


class TestRegexGen:
    def test_builds_expected_pattern(self):
        expected = (
            "G[^GSHSMRYFFT]*?S[^GSHSMRYFFT]*?H[^GSHSMRYFFT]*?S[^GSHSMRYFFT]*?"
            "M[^GSHSMRYFFT]*?R[^GSHSMRYFFT]*?Y[^GSHSMRYFFT]*?F[^GSHSMRYFFT]*?"
            "F[^GSHSMRYFFT]*?T[^GSHSMRYFFT]*?"
        )
        assert regex_gen("GSHSMRYFFT") == expected


class TestCoordinate:
    # A synthetic sequence: a leader peptide followed by the mature start motif.
    SEQ = "MABCDEFXYZGSHSMRYFFTQQQ"
    RGX = regex_gen("GSHSMRYFFT")

    def test_start_coordinate_of_mature_protein(self):
        # The motif "GSHSMRYFFT" begins at index 10.
        assert coordinate(self.SEQ, self.RGX) == 10

    def test_end_coordinate_of_domain(self):
        # The match spans through the end of the motif at index 20.
        assert coordinate_end(self.SEQ, self.RGX) == 20


@pytest.mark.network
class TestGetImgtDbVersions:
    def test_returns_known_versions(self):
        versions = get_imgt_db_versions()
        assert isinstance(versions, list)
        assert "3510" in versions
        assert "Latest" in versions
