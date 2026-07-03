"""Tests for hlagenie.db low-level SQLite helpers.

These run entirely against an in-memory SQLite database, so no network or
IMGT data is required.
"""

import sqlite3

import pytest

from hlagenie import db


@pytest.fixture
def conn():
    connection = sqlite3.connect(":memory:")
    yield connection
    connection.close()


class TestTableExistence:
    def test_missing_table(self, conn):
        assert db.table_exists(conn, "nope") is False

    def test_present_table(self, conn):
        db.save_dict(conn, "foo", {"a": "1"}, ("k", "v"))
        assert db.table_exists(conn, "foo") is True

    def test_tables_exist_all_present(self, conn):
        db.save_dict(conn, "foo", {"a": "1"}, ("k", "v"))
        db.save_dict(conn, "bar", {"b": "2"}, ("k", "v"))
        assert db.tables_exist(conn, ["foo", "bar"]) is True

    def test_tables_exist_one_missing(self, conn):
        db.save_dict(conn, "foo", {"a": "1"}, ("k", "v"))
        assert db.tables_exist(conn, ["foo", "bar"]) is False


class TestSaveAndLoadDict:
    def test_round_trip(self, conn):
        data = {"a": "1", "b": "2", "c": "3"}
        db.save_dict(conn, "t", data, ("k", "v"))
        assert db.load_dict(conn, "t", ("k", "v")) == data

    def test_count_rows(self, conn):
        db.save_dict(conn, "t", {"a": "1", "b": "2"}, ("k", "v"))
        assert db.count_rows(conn, "t") == 2

    def test_save_dict_overwrites_existing_table(self, conn):
        db.save_dict(conn, "t", {"a": "1", "b": "2"}, ("k", "v"))
        db.save_dict(conn, "t", {"z": "9"}, ("k", "v"))
        assert db.load_dict(conn, "t", ("k", "v")) == {"z": "9"}


class TestSaveAndLoadSet:
    def test_round_trip(self, conn):
        data = {"x", "y", "z"}
        db.save_set(conn, "s", data, ("col"))
        assert db.load_set(conn, "s", ("col")) == data


class TestUserVersion:
    def test_default_is_none(self, conn):
        assert db.get_user_version(conn) is None

    def test_set_then_get(self, conn):
        db.set_user_version(conn, 3510)
        assert db.get_user_version(conn) == 3510
