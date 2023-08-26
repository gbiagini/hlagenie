import copy
import functools
import sqlite3
import hlagenie.load
from .load import (
    load_sequence_alignment,
    load_nucleotide_alignment,
)
from hlagenie.configs import config
from hlagenie.smart_sort import smart_sort_comparator
from . import db
from .misc import find_gaps, regex_gen, coordinate, coordinate_end
import pyard  # for HLA nomenclature


def generate_gapped_tables(db_conn: sqlite3.Connection, imgt_version, imputed):
    """
    Create tables with gapped sequences for every allele in the IMGT/HLA database for each locus

    :param db_conn: The database connection object
    :return: dictionary of gapped sequences
    """

    # check if the tables exist so as to not rebuild if unnecessary
    if db.tables_exist(db_conn, config["gapped_tables"]):
        return db.load_gapped_tables(db_conn)

    # initialize pyard object
    ard = pyard.init(imgt_version)

    # initialize the dictionary to store all sequences
    gapped_seqs = {}

    # retrieve multiple sequence alignment for each locus
    for locus in config["loci"]:
        # load the sequence alignment
        multi_seq = load_sequence_alignment(imgt_version, locus, imputed)

        # turn the sequence alignment into a dictionary

        ## initialize a per-locus dictionary
        loc_seqs = {}

        ## iterate through the sequence alignment
        for record in multi_seq:
            ### use py-ard to get two-field allele
            allele = ard.redux(record.id, "U2")

            # for DRB345, skip if not specific locus
            if allele.split("*")[0] != locus:
                continue

            ### only add allele if not already present to avoid overwriting with less complete sequence
            if allele not in loc_seqs.keys():
                loc_seqs[allele] = str(record.seq)

        # save the sequence alignment to the database
        db.save_dict(db_conn, f"{locus}_gapped", loc_seqs, ("allele", "seq"))

        # update overall dictionary
        gapped_seqs.update(loc_seqs)

    return gapped_seqs


def generate_gapped_mature_tables(db_conn: sqlite3.Connection):
    """
    Create tables with gapped mature sequences for every allele in the IMGT/HLA database for each locus

    :return: dictionary of gapped mature sequences
    """

    # check if the tables exist so as to not rebuild if unnecessary
    if db.tables_exist(db_conn, config["gapped_mature_tables"]):
        return db.load_gapped_mature_tables(db_conn)

    # initialize dictionary to store all sequences
    mature_seqs = {}

    # retrieve gapped sequences for each locus
    for locus in config["loci"]:
        # get the gapped sequences
        locus_seqs = db.load_dict(db_conn, f"{locus}_gapped", ("allele", "seq"))

        # get the reference sequence
        ref_allele = config["refseq"][locus]
        for record, seq in locus_seqs.items():
            if record == ref_allele:
                ref_seq = seq
                break

        # get the regex based on the first ten residues
        regex = regex_gen(config["first_ten"][ref_allele])

        # get the start coordinates based on reference sequence
        start_coords = coordinate(ref_seq, regex)

        # initialize a per-locus dictionary
        loc_dict = {}

        # iterate through the sequences, removing the leader sequence
        for key, value in locus_seqs.items():
            loc_dict[key] = value[start_coords:]

        # save the dictionary to the database
        db.save_dict(db_conn, f"{locus}_gapped_mature", loc_dict, ("allele", "seq"))

        # update overall dictionary
        mature_seqs.update(loc_dict)

    return mature_seqs


# TODO - consider if this should leave positions which are simply unknown (current) or also remove these
def generate_ungapped_tables(db_conn: sqlite3.Connection, imgt_version, imputed):
    """
    Create tables with ungapped sequences for every allele in the IMGT/HLA database for each locus

    :param db_conn: The database connection object
    :return: dictionary of ungapped sequences
    """

    # check if the tables exist so as to not rebuild if unnecessary
    if db.tables_exist(db_conn, config["ungapped_tables"]):
        return db.load_ungapped_tables(db_conn)

    # initialize pyard object
    ard = pyard.init(imgt_version)

    # initialize the dictionary to store all sequences
    ungapped_seqs = {}

    # retrieve multiple sequence alignment for each locus
    for locus in config["loci"]:
        # load the sequence alignment
        multi_seq = load_sequence_alignment(imgt_version, locus, imputed)

        # get the reference sequence
        ref_allele = config["refseq"][locus]
        for record in multi_seq:
            if ard.redux(record.id, "U2") == ref_allele:
                ref_seq = str(record.seq)
                break

        # get the gaps in the reference sequence
        gaps = find_gaps(ref_seq)

        # turn the sequence alignment into a dictionary

        ## initialize a per-locus dictionary
        loc_seqs = {}

        ## iterate through the sequence alignment
        for record in multi_seq:
            ### use py-ard to get two-field allele
            allele = ard.redux(record.id, "U2")

            # for DRB345, skip if not specific locus
            if allele.split("*")[0] != locus:
                continue

            ### only add allele if not already present to avoid overwriting with less complete sequence
            if allele not in loc_seqs.keys():
                #### remove gaps from the sequence (if actually a gap)
                sequence = "".join(
                    [
                        char
                        for i, char in enumerate(str(record.seq))
                        if ((i not in gaps) or (char != "-"))
                    ]
                )

                loc_seqs[allele] = sequence

        # save the sequence alignment to the database
        db.save_dict(db_conn, f"{locus}_ungapped", loc_seqs, ("allele", "seq"))

        # update overall dictionary
        ungapped_seqs.update(loc_seqs)

    return ungapped_seqs


def generate_completed_table(
    db_conn: sqlite3.Connection, locus: str, seqtype: str, seqs: dict
):
    """
    Create table with list of allele with completed sequences

    :param db_conn: SQLite3 database connection object to HLAGenie database
    :type db_conn: sqlite3.Connection
    :param locus: HLA locus to generate completed table for
    :type locus: str
    :param seqtype: sequence type to generate completed table for
    :type seqtype: str
    :param seqs: dictionary of sequences
    :type seqs: dict
    :return: list of alleles with completed sequences
    """

    if db.table_exists(db_conn, f"{locus}_completed_{seqtype}"):
        return db.load_set(db_conn, f"{locus}_completed_{seqtype}", ("allele"))

    if seqtype == "prot":
        # get the reference sequence
        ref_allele = config["refseq"][locus]
    else:
        # get the reference sequence
        ref_allele = config["refseq_full"][locus]

    # get the reference sequence
    ref_seq = seqs[ref_allele]

    # get the gaps in the reference sequence
    gaps = set(find_gaps(ref_seq))

    # get the sequences for the locus
    loc_seqs = {
        allele: seq for allele, seq in seqs.items() if allele.split("*")[0] == locus
    }

    # initialize list to store completed alleles
    completed = [
        allele for allele, seq in loc_seqs.items() if set(find_gaps(seq)) == gaps
    ]

    # remove null alleles
    completed = [allele for allele in completed if not allele[-1].isalpha()]

    # save the list to the database
    db.save_set(db_conn, f"{locus}_completed_{seqtype}", set(completed), ("allele"))

    return completed


def generate_incomplete_table(
    db_conn: sqlite3.Connection, locus: str, seqtype: str, seqs: dict
):
    """
    Create table with list of allele with incomplete sequences

    :param db_conn: SQLite3 database connection object to HLAGenie database
    :type db_conn: sqlite3.Connection
    :param locus: HLA locus to generate incomplete table for
    :type locus: str
    :param seqtype: sequence type to generate incomplete table for
    :type seqtype: str
    :param seqs: dictionary of sequences
    :type seqs: dict
    :return: list of alleles with incomplete sequences
    """

    if db.table_exists(db_conn, f"{locus}_incomplete_{seqtype}"):
        return db.load_set(db_conn, f"{locus}_incomplete_{seqtype}", ("allele"))

    if seqtype == "prot":
        # get the reference sequence
        ref_allele = config["refseq"][locus]
    else:
        # get the reference sequence
        ref_allele = config["refseq_full"][locus]

    # get the reference sequence
    ref_seq = seqs[ref_allele]

    # get the gaps in the reference sequence
    gaps = len(find_gaps(ref_seq))

    # limit to alleles in locus
    loc_seqs = {
        allele: seq for allele, seq in seqs.items() if allele.split("*")[0] == locus
    }

    # initialize list to store incomplete alleles
    incomplete = [
        allele for allele, seq in loc_seqs.items() if len(find_gaps(seq)) > gaps
    ]

    # remove null alleles
    incomplete = [allele for allele in incomplete if not allele[-1].isalpha()]

    # save the list to the database
    db.save_set(db_conn, f"{locus}_incomplete_{seqtype}", set(incomplete), ("allele"))

    return incomplete


def generate_extended_table(
    db_conn: sqlite3.Connection, locus: str, seqtype: str, ungap: bool, seqs: dict
):
    """
    Create table with list of allele with extended sequences

    :param db_conn: SQLite3 database connection object to HLAGenie database
    :type db_conn: sqlite3.Connection
    :param locus: HLA locus to generate extended table for
    :type locus: str
    :param seqtype: sequence type to generate extended table for
    :type seqtype: str
    :param ungap: whether ungapped sequences were used
    :type ungap: bool
    :param seqs: dictionary of sequences
    :type seqs: dict
    :return: list of alleles with extended sequences
    """

    if db.table_exists(db_conn, f"{locus}_extended_{seqtype}"):
        return db.load_set(db_conn, f"{locus}_extended_{seqtype}", ("allele"))

    if seqtype == "prot":
        # get the reference sequence
        ref_allele = config["refseq"][locus]
    else:
        # get the reference sequence
        ref_allele = config["refseq_full"][locus]

    # get the reference sequence
    ref_seq = seqs[ref_allele]

    # limit to alleles in locus
    loc_seqs = {
        allele: seq for allele, seq in seqs.items() if allele.split("*")[0] == locus
    }

    # check if ungap=True was used
    if ungap:
        # get the length of the reference sequence
        ref_len = len(ref_seq)

        # initialize list to store extended alleles
        extended = [allele for allele, seq in loc_seqs.items() if len(seq) > ref_len]
    else:
        # get the gaps in the reference sequence
        gaps = len(find_gaps(ref_seq))
        # initialize list to store extended alleles
        extended = [
            allele for allele, seq in loc_seqs.items() if len(find_gaps(seq)) < gaps
        ]

    # save the list to the database
    db.save_set(db_conn, f"{locus}_extended_{seqtype}", set(extended), ("allele"))

    return extended


def generate_ungapped_mature_tables(db_conn: sqlite3.Connection):
    """
    Create tables with ungapped mature sequences for every allele in the IMGT/HLA database for each locus

    :return: dictionary of ungapped mature sequences
    """

    # check if the tables exist so as to not rebuild if unnecessary
    if db.tables_exist(db_conn, config["ungapped_mature_tables"]):
        return db.load_ungapped_mature_tables(db_conn)

    # initialize dictionary to store all sequences
    mature_seqs = {}

    # retrieve gapped sequences for each locus
    for locus in config["loci"]:
        # get the gapped sequences
        locus_seqs = db.load_dict(db_conn, f"{locus}_ungapped", ("allele", "seq"))

        # get the reference sequence
        ref_allele = config["refseq"][locus]
        for record, seq in locus_seqs.items():
            if record == ref_allele:
                ref_seq = seq
                break

        # get the regex based on the first ten residues
        regex = regex_gen(config["first_ten"][ref_allele])

        # get the start coordinates based on reference sequence
        start_coords = coordinate(ref_seq, regex)

        # initialize a per-locus dictionary
        loc_dict = {}

        # iterate through the sequences, removing the leader sequence
        for key, value in locus_seqs.items():
            loc_dict[key] = value[start_coords:]

        # save the dictionary to the database
        db.save_dict(db_conn, f"{locus}_ungapped_mature", loc_dict, ("allele", "seq"))

        # update overall dictionary
        mature_seqs.update(loc_dict)

    return mature_seqs


def generate_gapped_ard_table(db_conn: sqlite3.Connection, seqs: dict):
    """
    Generate a table with the ARD ending position for each locus

    :param db_conn: The database connection object
    :param seqs: dictionary of sequences to get reference sequence
    :return: dictionary of ARD ending positions
    """

    # check if the table exists so as to not rebuild if unnecessary
    if db.table_exists(db_conn, "gapped_ard"):
        ard_ends = db.load_dict(db_conn, "gapped_ard", ("locus", "ard_end"))
        # convert the values to integers
        ard_ends = dict(map(lambda x: (x[0], int(x[1])), ard_ends.items()))
        return ard_ends

    # initialize dictionary to store ard positions
    ard_ends = {}

    # iterate through loci
    for locus in config["loci"]:
        # get the reference allele
        ref_allele = config["refseq"][locus]

        # get the regex based on the last ten residues of the ARD
        regex = regex_gen(config["ard_last_ten"][ref_allele])

        # get the reference sequence
        ref_seq = seqs[ref_allele]

        # get the end coordinates based on reference sequence
        end_coords = coordinate_end(ref_seq, regex)

        # add to dictionary
        ard_ends[locus] = int(end_coords)

    # save the dictionary to the database
    db.save_dict(db_conn, "gapped_ard", ard_ends, ("locus", "ard_end"))

    return ard_ends


def generate_gapped_xrd_table(db_conn: sqlite3.Connection, seqs: dict):
    """
    Generate a table with the XRD ending position for each locus

    :param db_conn: The database connection object
    :param seqs: dictionary of sequences to get reference sequence
    :return: dictionary of XRD ending positions
    """

    # check if the table exists so as to not rebuild if unnecessary
    if db.table_exists(db_conn, "gapped_xrd"):
        xrd_ends = db.load_dict(db_conn, "gapped_xrd", ("locus", "xrd_end"))
        # convert the values to integers
        xrd_ends = dict(map(lambda x: (x[0], int(x[1])), xrd_ends.items()))
        return xrd_ends

    # initialize dictionary to store ard positions
    xrd_ends = {}

    # iterate through loci
    for locus in config["loci"]:
        # get the reference allele
        ref_allele = config["refseq"][locus]

        # get the regex based on the last ten residues of the XRD
        regex = regex_gen(config["xrd_last_ten"][ref_allele])

        # get the reference sequence
        ref_seq = seqs[ref_allele]

        # get the end coordinates based on reference sequence
        end_coords = coordinate_end(ref_seq, regex)

        # add to dictionary
        xrd_ends[locus] = int(end_coords)

    # save the dictionary to the database
    db.save_dict(db_conn, "gapped_xrd", xrd_ends, ("locus", "xrd_end"))

    return xrd_ends


def generate_ungapped_ard_table(db_conn: sqlite3.Connection, seqs: dict):
    """
    Generate a table with the ARD ending position for each locus

    :param db_conn: The database connection object
    :param seqs: dictionary of sequences to get reference sequence
    :return: dictionary of ARD ending positions
    """

    # check if the table exists so as to not rebuild if unnecessary
    if db.table_exists(db_conn, "ungapped_ard"):
        ard_ends = db.load_dict(db_conn, "ungapped_ard", ("locus", "ard_end"))
        # convert the values to integers
        ard_ends = dict(map(lambda x: (x[0], int(x[1])), ard_ends.items()))
        return ard_ends

    # initialize dictionary to store ard positions
    ard_ends = {}

    # iterate through loci
    for locus in config["loci"]:
        # get the reference allele
        ref_allele = config["refseq"][locus]

        # get the regex based on the last ten residues of the ARD
        regex = regex_gen(config["ard_last_ten"][ref_allele])

        # get the reference sequence
        ref_seq = seqs[ref_allele]

        # get the end coordinates based on reference sequence
        end_coords = coordinate_end(ref_seq, regex)

        # add to dictionary
        ard_ends[locus] = int(end_coords)

    # save the dictionary to the database
    db.save_dict(db_conn, "ungapped_ard", ard_ends, ("locus", "ard_end"))

    return ard_ends


def generate_ungapped_xrd_table(db_conn: sqlite3.Connection, seqs: dict):
    """
    Generate a table with the XRD ending position for each locus

    :param db_conn: The database connection object
    :param seqs: dictionary of sequences to get reference sequence
    :return: dictionary of XRD ending positions
    """

    # check if the table exists so as to not rebuild if unnecessary
    if db.table_exists(db_conn, "ungapped_xrd"):
        xrd_ends = db.load_dict(db_conn, "ungapped_xrd", ("locus", "xrd_end"))
        # convert the values to integers
        xrd_ends = dict(map(lambda x: (x[0], int(x[1])), xrd_ends.items()))
        return xrd_ends

    # initialize dictionary to store ard positions
    xrd_ends = {}

    # iterate through loci
    for locus in config["loci"]:
        # get the reference allele
        ref_allele = config["refseq"][locus]

        # get the regex based on the last ten residues of the XRD
        regex = regex_gen(config["xrd_last_ten"][ref_allele])

        # get the reference sequence
        ref_seq = seqs[ref_allele]

        # get the end coordinates based on reference sequence
        end_coords = coordinate_end(ref_seq, regex)

        # add to dictionary
        xrd_ends[locus] = int(end_coords)

    # save the dictionary to the database
    db.save_dict(db_conn, "ungapped_xrd", xrd_ends, ("locus", "xrd_end"))

    return xrd_ends


def generate_ungapped_nuc_tables(db_conn: sqlite3.Connection, imgt_version, imputed):
    """Generate a table with ungapped nucleotide sequences for every locus

    :param db_conn: The database connection object
    :type db_conn: sqlite3.Connection
    :return: dictionary of ungapped nucleotide sequences
    """
    # check if the tables exist so as to not rebuild if unnecessary
    if db.tables_exist(db_conn, config["ungapped_nuc_tables"]):
        return db.load_ungapped_nuc_tables(db_conn)

    # initialize pyard object
    ard = pyard.init(imgt_version)

    # initialize the dictionary to store all sequences
    ungapped_seqs = {}

    # retrieve multiple sequence alignment for each locus
    for locus in config["loci"]:
        # load the nucleotide sequence alignment
        multi_seq = load_nucleotide_alignment(imgt_version, locus, imputed)

        # get the reference sequence
        ref_allele = config["refseq_full"][locus]
        for record in multi_seq:
            if record.id == ref_allele:
                ref_seq = str(record.seq)
                break

        # get the gaps in the reference sequence
        gaps = find_gaps(ref_seq)

        # turn the sequence alignment into a dictionary

        ## initialize a per-locus dictionary
        loc_seqs = {}

        ## iterate through the sequence alignment
        for record in multi_seq:
            allele = record.id

            # for DRB345, skip if not specific locus
            if allele.split("*")[0] != locus:
                continue

            ### only add allele if not already present to avoid overwriting with less complete sequence
            if allele not in loc_seqs.keys():
                #### remove gaps from the sequence (if actually a gap)
                sequence = "".join(
                    [
                        char
                        for i, char in enumerate(str(record.seq))
                        if ((i not in gaps) or (char != "-"))
                    ]
                )

                loc_seqs[allele] = sequence

        # save the sequence alignment to the database
        db.save_dict(db_conn, f"{locus}_ungapped_nuc", loc_seqs, ("allele", "seq"))

        # update overall dictionary
        ungapped_seqs.update(loc_seqs)

    return ungapped_seqs


def generate_gapped_nuc_tables(db_conn: sqlite3.Connection, imgt_version, imputed):
    """
    Create tables with gapped nucleotide sequences for every allele in the IMGT/HLA database for each locus

    :param db_conn: The database connection object
    :return: dictionary of gapped sequences
    """

    # check if the tables exist so as to not rebuild if unnecessary
    if db.tables_exist(db_conn, config["gapped_nuc_tables"]):
        return db.load_gapped_nuc_tables(db_conn)

    # initialize pyard object
    ard = pyard.init(imgt_version)

    # initialize the dictionary to store all sequences
    gapped_seqs = {}

    # retrieve multiple sequence alignment for each locus
    for locus in config["loci"]:
        # load the sequence alignment
        multi_seq = load_nucleotide_alignment(imgt_version, locus, imputed)

        # turn the sequence alignment into a dictionary

        ## initialize a per-locus dictionary
        loc_seqs = {}

        ## iterate through the sequence alignment
        for record in multi_seq:
            allele = record.id

            # for DRB345, skip if not specific locus
            if allele.split("*")[0] != locus:
                continue

            ### only add allele if not already present to avoid overwriting with less complete sequence
            if allele not in loc_seqs.keys():
                loc_seqs[allele] = str(record.seq)

        # save the sequence alignment to the database
        db.save_dict(db_conn, f"{locus}_gapped_nuc", loc_seqs, ("allele", "seq"))

        # update overall dictionary
        gapped_seqs.update(loc_seqs)

    return gapped_seqs


def set_db_version(db_connection: sqlite3.Connection, imgt_version):
    """
    Set the IMGT database version number as a user_version string in
    the database itself.

    :param db_connection: Active SQLite Connection
    :param imgt_version: current imgt_version
    """
    # If version already exists, don't reset
    version = db.get_user_version(db_connection)
    if version:
        return version

    db.set_user_version(db_connection, int(imgt_version))
    print("Version:", imgt_version)
    return imgt_version
