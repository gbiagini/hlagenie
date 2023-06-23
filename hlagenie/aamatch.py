#!/usr/bin/env python3

# aa_matching.py - module for amino acid matching functions

# import necessary modules
from pathlib import Path  # for path manipulation
import pyard  # for HLA nomenclature
import yaml  # for reading YAML files
import db_ops  # for database operations


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
        self.db_connection = db_ops.create_db_connection(data_dir, imgt_version)

    def getAAposition(self, allele, position):
        """
        Get the amino acid at a specific position in an allele

        :param allele: The allele to get the amino acid from
        :param position: The position to get the amino acid from
        """

        # get the amino acid at the specified position
        return self.seqs[allele][position - 1]

    def getAAsubstring(self, allele, start, stop):
        """
        Get the amino acid substring from a specified position to another specified position

        :param allele: The allele to get the amino acid substring from
        :param start: The position to start the substring
        :param stop: The position to end the substring
        """

        # get the amino acid substring
        return self.seqs[allele][start - 1 : stop]

    def getEpitope(self, allele, positions):
        """
        Get the epitope string from a list of positions

        :param allele: The allele to get the epitope from
        :param positions: A list of positions to retrieve the epitope from
        """

        # get the epitope string
        return "_".join(
            [f"{position}{self.seqs[allele][position-1]}" for position in positions]
        )

    def isPositionMismatched(self, allele1, allele2, position):
        """
        Check if two alleles have a mismatch at a specified position

        :param allele1: The first allele to check
        :param allele2: The second allele to check
        :param position: The position to check
        """

        # get the amino acid at the specified position for each allele
        aa1 = self.getAAposition(allele1, position)
        aa2 = self.getAAposition(allele2, position)

        # check if the amino acids are the same
        return aa1 == aa2

    def countAAMismatchesAllele(
        self, allele1donor, allele2donor, allele1recip, allele2recip, position
    ):
        """
        Count the number of amino acid mismatches between two alleles at a specified position, adjusting for donor homozygosity

        :param allele1donor: The first allele of the donor to check
        :param allele2donor: The second allele of the donor to check
        :param allele1recip: The first allele of the recipient to check
        :param allele2recip: The second allele of the recipient to check
        :param position: The position to check
        """

        # check if donor is homozygous
        donor_homozygous = False
        if allele1donor == allele2donor:
            donor_homozygous = True

        # get amino acids at specified position
        aa1_donor = self.getAAposition(allele1donor, position)
        aa2_donor = self.getAAposition(allele2donor, position)
        aa1_recip = self.getAAposition(allele1recip, position)
        aa2_recip = self.getAAposition(allele2recip, position)

        # count mismatches between donor and recipient
        mm_count = self.countAAMismatches(aa1_donor, aa2_donor, aa1_recip, aa2_recip)

        # adjust if donor is homozygous, due to mismatch being same AA
        if (mm_count == 2) and (donor_homozygous):
            mm_count = 1

        return mm_count
