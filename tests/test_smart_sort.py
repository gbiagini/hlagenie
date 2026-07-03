"""Tests for hlagenie.smart_sort.smart_sort_comparator.

Pure function, no network or database required.
"""

from hlagenie.smart_sort import smart_sort_comparator as ss


def test_identical_alleles_return_zero():
    assert ss("A*01:01", "A*01:01") == 0


def test_first_field_ordering():
    assert ss("A*01:01", "A*02:01") == -1
    assert ss("A*02:01", "A*01:01") == 1


def test_second_field_ordering():
    assert ss("A*01:01", "A*01:02") == -1
    assert ss("A*01:02", "A*01:01") == 1


def test_numeric_not_lexicographic():
    # Lexicographically "A*01:100" < "A*01:9"; numerically it must be greater.
    assert ss("A*01:100", "A*01:09") == 1


def test_third_and_fourth_fields():
    assert ss("A*01:01:01", "A*01:01:02") == -1
    assert ss("A*01:01:01:01", "A*01:01:01:02") == -1
    # Everything past the fourth field is considered equal.
    assert ss("A*01:01:01:01", "A*01:01:01:01") == 0


def test_expression_characters_are_stripped_before_compare():
    # The expression suffix (here 'N') is removed, so these are equal alleles.
    assert ss("A*01:01N", "A*01:01") == 0


def test_serology_compared_lexicographically():
    assert ss("A1", "A2") == -1
    assert ss("A2", "A1") == 1


def test_glstring_compared_lexicographically():
    # Presence of a GL-string character short-circuits to string comparison.
    assert ss("A*01:01/A*01:02", "A*02:01") == -1
    assert ss("A*02:01", "A*01:01/A*01:02") == 1


def test_expression_suffix_with_extra_field():
    assert ss("B*44:02:01:02S", "B*44:02") == 1
