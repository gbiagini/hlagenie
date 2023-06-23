#!/usr/bin/env python3

# aa_matching.py - module for amino acid matching functions

# import necessary modules
from pathlib import Path  # for path manipulation
import pyard  # for ARD reduction
import yaml  # for reading YAML files


class genie:
    """
    Matching functions for HLA alleles
    Allows for matching of HLA alleles at the amino acid level
    """

    def __init__(
        self,
        imgt_version: str = "Latest",
        data_dir: str = None,
        max_cache_size: int = DEFAULT_CACHE_SIZE,
        config: dict = None,
        ungap=True,
    ):
        # set values for needed variables
        self._data_dir = data_dir

        # create database connection to SQLite database
        self.db_connection = db.create_db_connection(data_dir, imgt_version)
