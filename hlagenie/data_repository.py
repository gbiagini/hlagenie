import copy
import functools
import sqlite3
import hlagenie.load
from .load import load_sequence_alignment
from hlagenie.configs import config
from hlagenie.smart_sort import smart_sort_comparator
from . import db
from .misc import find_gaps
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
            if allele not in gapped_seqs.keys():
                loc_seqs[allele] = str(record.seq)

        # save the sequence alignment to the database
        db.save_dict(db_conn, f"{locus}_gapped", loc_seqs, ("allele", "seq"))

        # update overall dictionary
        gapped_seqs.update(loc_seqs)

    return gapped_seqs


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
            if record.id == ref_allele:
                ref_seq = str(record.seq)

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
            if allele not in ungapped_seqs.keys():
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
