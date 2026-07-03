#!/usr/bin/env python3

# aa_matching.py - module for amino acid matching functions

# import necessary modules
import pyard  # for HLA nomenclature
from . import db  # for database operations
from . import data_repository as dr  # for data repository operations
from .load import load_latest_version  # get most updated version of IMGT database
from .configs import config  # for configurations


class GENIE:
    """
    Matching functions for HLA alleles
    Allows for matching of HLA alleles at the amino acid level
    """

    def __init__(
        self,
        imgt_version: str = "Latest",
        data_dir: str = None,
        load_mac: bool = True,
        cache_size: int = config["DEFAULT_CACHE_SIZE"],
        ungap: bool = True,
        imputed: bool = False,
        imputation_method: str = "nearest",
    ):
        # set values for needed variables
        self._data_dir = data_dir
        self.ungap = ungap
        self.load_mac = load_mac

        # py-ard object, initialized lazily on first multi-field allele lookup
        # so the cached fast path never pays py-ard's initialization cost
        self.ard = None

        # if database version is "Latest", get the latest version
        if imgt_version == "Latest":
            imgt_version = load_latest_version()

        self.imgt_version = imgt_version

        # make sure imputation method is valid
        if imputation_method not in ["nearest", "nearest10", "nearest10_weighted"]:
            raise ValueError(
                "Invalid imputation method specified, must be nearest or nearest10"
            )

        # create database connection to SQLite database
        self.db_connection = db.create_db_connection(
            data_dir, imgt_version, imputed, imputation_method
        )

        # save the IMGT version
        dr.set_db_version(self.db_connection, imgt_version)

        # load sequence data from database
        if self.ungap:
            self.full_seqs = dr.generate_ungapped_tables(
                self.db_connection, imgt_version, imputed, imputation_method, self.load_mac
            )
            self.nuc_seqs = dr.generate_ungapped_nuc_tables(
                self.db_connection, imgt_version, imputed, imputation_method
            )
            self.seqs = dr.generate_ungapped_mature_tables(self.db_connection)
            self.ards = dr.generate_ungapped_ard_table(self.db_connection, self.seqs)
            self.xrds = dr.generate_ungapped_xrd_table(self.db_connection, self.seqs)
        else:
            self.full_seqs = dr.generate_gapped_tables(
                self.db_connection, imgt_version, imputed, imputation_method, self.load_mac
            )
            self.nuc_seqs = dr.generate_gapped_nuc_tables(
                self.db_connection, imgt_version, imputed, imputation_method
            )
            self.seqs = dr.generate_gapped_mature_tables(self.db_connection)
            self.ards = dr.generate_gapped_ard_table(self.db_connection, self.seqs)
            self.xrds = dr.generate_gapped_xrd_table(self.db_connection, self.seqs)

    def __del__(self):
        """Close the db connection, when HLAGenie instance goes away

        :return:
        """
        if hasattr(self, "db_connection") and self.db_connection:
            self.db_connection.close()

    def _reduce(self, allele: str) -> str:
        """Reduce a >2-field allele name to its two-field form via py-ard.

        Names that are already two-field (or coarser) are returned unchanged.
        The py-ard object is created on first use and reused thereafter.

        :param allele: allele name to reduce
        :return: two-field allele name
        """

        if allele.count(":") <= 1:
            return allele

        if self.ard is None:
            self.ard = pyard.init(self.imgt_version, load_mac=self.load_mac)

        return self.ard.redux(allele, "U2")

    def getAA(self, allele: str, position: int):
        """
        Get the amino acid at a specific position in an allele

        :param allele: The allele to get the amino acid from
        :param position: The position to get the amino acid from
        :return: The amino acid at the specified position
        """

        allele = self._reduce(allele)

        # get the amino acid at the specified position
        return self.seqs[allele][position - 1]

    def getNuc(self, allele: str, position: int):
        """Get the nucleotide at a specific position in an allele

        Unlike the amino-acid lookups, a full (3- or 4-field) allele name is
        required: a single two-field name maps to many distinct nucleotide
        sequences (synonymous differences at the third field), so the name is
        not reduced. Alleles whose only name is two-field are looked up directly.

        :param allele: The full allele name to get the nucleotide from
        :type allele: str
        :param position: The position to get the nucleotide from
        :type position: int
        :return: The nucleotide at the specified position
        :raises ValueError: if the allele has no nucleotide sequence (e.g. a
            bare two-field name that maps to multiple sequences, or an unknown allele)
        """

        if allele not in self.nuc_seqs:
            raise ValueError(
                f"No nucleotide sequence for {allele!r}. getNuc requires a full "
                "(3- or 4-field) allele name, because one two-field name maps to "
                "multiple distinct nucleotide sequences."
            )

        # get the nucleotide at the specified position
        return self.nuc_seqs[allele][position - 1]

    def getPeptide(self, allele: str, start: int, stop: int):
        """
        Get the amino acid substring from a specified position to another specified position

        :param allele: The allele to get the amino acid substring from
        :param start: The position to start the substring
        :param stop: The position to end the substring
        :return: The amino acid substring from the specified positions
        """

        allele = self._reduce(allele)

        # get the amino acid substring
        return self.seqs[allele][start - 1 : stop]

    def getEpitope(self, allele: str, positions: list[int]):
        """
        Get the epitope string from a list of positions

        :param allele: The allele to get the epitope from
        :param positions: A list of positions to retrieve the epitope from
        :return: The epitope string from the specified positions
        """

        allele = self._reduce(allele)

        # get the epitope string
        return "_".join(
            [f"{position}{self.seqs[allele][position-1]}" for position in positions]
        )

    def isPositionMismatched(self, allele1: str, allele2: str, position: int):
        """
        Check if two alleles have a mismatch at a specified position

        :param allele1: The first allele to check
        :param allele2: The second allele to check
        :param position: The position to check
        :return: True if the alleles have a mismatch at the specified position, False otherwise
        """

        allele1 = self._reduce(allele1)
        allele2 = self._reduce(allele2)

        # get the amino acid at the specified position for each allele
        aa1 = self.getAA(allele1, position)
        aa2 = self.getAA(allele2, position)

        # check if the amino acids are the same
        return not (aa1 == aa2)

    def countAAMismatchesAllele(
        self,
        allele1donor: str,
        allele2donor: str,
        allele1recip: str,
        allele2recip: str,
        position: int,
    ):
        """
        Count the number of amino acid mismatches between two alleles at a specified position, adjusting for donor homozygosity

        :param allele1donor: The first allele of the donor to check
        :param allele2donor: The second allele of the donor to check
        :param allele1recip: The first allele of the recipient to check
        :param allele2recip: The second allele of the recipient to check
        :param position: The position to check
        :return: The number of amino acid mismatches between the two alleles at the specified position
        """
        allele1donor = self._reduce(allele1donor)
        allele2donor = self._reduce(allele2donor)
        allele1recip = self._reduce(allele1recip)
        allele2recip = self._reduce(allele2recip)

        # check if donor is homozygous
        donor_homozygous = False
        if allele1donor == allele2donor:
            donor_homozygous = True

        # get amino acids at specified position
        aa1_donor = self.getAA(allele1donor, position)
        aa2_donor = self.getAA(allele2donor, position)
        aa1_recip = self.getAA(allele1recip, position)
        aa2_recip = self.getAA(allele2recip, position)

        # count mismatches between donor and recipient
        mm_count = self.countAAMismatches(aa1_donor, aa2_donor, aa1_recip, aa2_recip)

        # adjust if donor is homozygous, due to mismatch being same AA
        if (mm_count == 2) and (donor_homozygous):
            mm_count = 1

        return mm_count

    def countAAMismatches(
        self, aa1_donor: str, aa2_donor: str, aa1_recip: str, aa2_recip: str
    ):
        """
        Count the number of amino acid mismatches at a position bewteen donor and recipient

        :param aa1_donor: The first amino acid of the donor to check
        :param aa2_donor: The second amino acid of the donor to check
        :param aa1_recip: The first amino acid of the recipient to check
        :param aa2_recip: The second amino acid of the recipient to check
        :return: The number of amino acid mismatches
        """

        # count mismatches between donor and recipient
        mm_count = 0
        if (aa1_donor != aa1_recip) and (aa1_donor != aa2_recip):
            mm_count += 1
        if (aa2_donor != aa2_recip) and (aa2_donor != aa1_recip):
            mm_count += 1

        return mm_count

    def getARD(self, allele: str):
        """
        Get the ARD sequence of an allele

        :param allele: The allele to get the ARD sequence from
        :return: The ARD sequence
        """

        # reduce to two field if greater than two field
        allele = self._reduce(allele)

        # get locus
        locus = allele.split("*")[0]

        # get the ARD sequence
        return self.seqs[allele][: self.ards[locus]]

    def getXRD(self, allele: str):
        """
        Get the XRD sequence of an allele

        :param allele: The allele to get the XRD sequence from
        :return: The XRD sequence
        """

        # reduce to two field if greater than two field
        allele = self._reduce(allele)

        # get locus
        locus = allele.split("*")[0]

        # get the ARD sequence
        return self.seqs[allele][: self.xrds[locus]]

    def listIncompletes(self, locus: str, seqtype: str = "prot"):
        """
        List the incomplete alleles in the database

        :param locus: The locus to list the incomplete alleles from
        :param seqtype: The sequence type to list the incomplete alleles for (prot or nuc)
        :return: A list of incomplete alleles
        """

        if seqtype == "prot":
            return dr.generate_incomplete_table(
                self.db_connection, locus, seqtype, self.seqs
            )
        elif seqtype == "nuc":
            return dr.generate_incomplete_table(
                self.db_connection, locus, seqtype, self.nuc_seqs
            )
        else:
            print("Invalid sequence type specified")

    def listCompletes(self, locus: str, seqtype: str = "prot"):
        """
        List the complete alleles in the database

        :param locus: The locus to list the complete alleles from
        :param seqtype: The sequence type to list the complete alleles for (prot or nuc)
        :return: A list of complete alleles
        """

        if seqtype == "prot":
            return dr.generate_completed_table(
                self.db_connection, locus, seqtype, self.seqs
            )
        elif seqtype == "nuc":
            return dr.generate_completed_table(
                self.db_connection, locus, seqtype, self.nuc_seqs
            )
        else:
            print("Invalid sequence type specified")

    def listExtendeds(self, locus: str, seqtype: str = "prot"):
        """
        List the extended alleles in the database

        :param locus: The locus to list the extended alleles from
        :param seqtype: The sequence type to list the extended alleles for (prot or nuc)
        :return: A list of extended alleles
        """

        if seqtype == "prot":
            return dr.generate_extended_table(
                self.db_connection, locus, seqtype, self.ungap, self.seqs
            )
        elif seqtype == "nuc":
            return dr.generate_extended_table(
                self.db_connection, locus, seqtype, self.ungap, self.nuc_seqs
            )
        else:
            print("Invalid sequence type specified")
