import copy
import functools
import sqlite3
import hlagenie.load
from .load import load_sequence_alignment
from hlagenie.configs import config
from hlagenie.smart_sort import smart_sort_comparator
from . import db
from .misc import find_gaps, regex_gen, coordinate, coordinate_end
import pyard  # for HLA nomenclature


def generate_gapped_tables(db_conn: sqlite3.Connection, imgt_version):
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
        multi_seq = load_sequence_alignment(imgt_version, locus)

        # turn the sequence alignment into a dictionary

        ## initialize a per-locus dictionary
        loc_seqs = {}

        ## iterate through the sequence alignment
        for record in multi_seq:
            ### use py-ard to get two-field allele
            allele = ard.redux(record.id, "U2")

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
def generate_ungapped_tables(db_conn: sqlite3.Connection, imgt_version):
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
        multi_seq = load_sequence_alignment(imgt_version, locus)

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
        return db.load_dict(db_conn, "gapped_ard", ("locus", "ard_end"))

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
        ard_ends[locus] = end_coords

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
        return db.load_dict(db_conn, "gapped_xrd", ("locus", "xrd_end"))

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
        xrd_ends[locus] = end_coords

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
        return db.load_dict(db_conn, "ungapped_ard", ("locus", "ard_end"))

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
        ard_ends[locus] = end_coords

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
        return db.load_dict(db_conn, "ungapped_xrd", ("locus", "xrd_end"))

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
        xrd_ends[locus] = end_coords

    # save the dictionary to the database
    db.save_dict(db_conn, "ungapped_xrd", xrd_ends, ("locus", "xrd_end"))

    return xrd_ends
