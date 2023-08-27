import sqlite3  # for database operations
import pathlib  # for path manipulation
from .misc import (
    get_imgt_db_versions,
    get_default_db_directory,
)  # for getting default database directory
from .load import load_latest_version
from hlagenie.configs import config  # configurations


def create_db_connection(data_dir, imgt_version, imputed, imputation_method):
    """
    Create connection to SQLite database

    :param data_dir: The directory where the database is stored
    :param imgt_version: The version of the IMGT/HLA database to use
    :return: The database connection
    """

    # set data directory
    if data_dir is None:
        data_dir = get_default_db_directory()

    # set database filename
    if imputed:
        db_filename = (
            f"{data_dir}/hlagenie-{imgt_version}-imputed-{imputation_method}.db"
        )
    else:
        db_filename = f"{data_dir}/hlagenie-{imgt_version}.db"

    # Check if imgt_version is valid
    # if not pathlib.Path(db_filename).exists():
    #    all_versions = get_imgt_db_versions()
    #    if str(imgt_version) not in all_versions:
    #        raise ValueError(f"{imgt_version} is not a valid IMGT version")

    # Create the data directory if it doesn't exist
    if not pathlib.Path(data_dir).exists():
        pathlib.Path(data_dir).mkdir(parents=True, exist_ok=True)

    # Create the database file and alert user
    if not pathlib.Path(db_filename).exists():
        print(f"Creating database file {db_filename} as cache")

    # Open the database connection
    file_uri = f"file:{db_filename}"
    return sqlite3.connect(file_uri, uri=True)


def table_exists(connection: sqlite3.Connection, table_name: str) -> bool:
    """
    Does the table exist in the database ?

    :param connection: db connection of type sqlite.Connection
    :param table_name: table in the sqlite db
    :return: bool indicating whether table_name exists as a table
    """
    query = """SELECT count(*) from sqlite_master where type = 'table' and name = ?"""
    cursor = connection.execute(query, (table_name,))
    result = cursor.fetchone()
    cursor.close()
    return result[0] > 0


def tables_exist(connection: sqlite3.Connection, table_names: list[str]):
    """
    Do all the given tables exist in the database ?

    :param connection: db connection of type sqlite.Connection
    :param table_names: names of tables in the sqlite db
    :return: bool indicating whether all table_names exists
    """
    return all([table_exists(connection, table_name) for table_name in table_names])


def count_rows(connection: sqlite3.Connection, table_name: str) -> int:
    """
    Count number of rows in the table.

    :param connection: db connection of type sqlite.Connection
    :param table_name: table in the sqlite db
    :return: bool indicating whether table_name exists as a table
    """
    query = f"SELECT count(*) from '{table_name}'"
    cursor = connection.execute(query)
    result = cursor.fetchone()
    cursor.close()
    return result[0]


def save_dict(
    connection: sqlite3.Connection,
    table_name: str,
    dictionary: dict[str, str],
    columns: tuple[str, str],
) -> bool:
    """
    Save the dictionary as a table with column names from columns Tuple.

    :param connection: db connection of type sqlite.Connection
    :param table_name: name of the table to create
    :param dictionary: the dictionary which to take the key and values from
    :param columns: column names in the table
    :return: success status
    """
    cursor = connection.cursor()

    # Drop the table first
    drop_table_sql = f"DROP TABLE IF EXISTS {table_name}"
    cursor.execute(drop_table_sql)

    # Create table
    create_table_sql = f"""CREATE TABLE {table_name} (
                            {columns[0]} TEXT PRIMARY KEY,
                            {columns[1]} TEXT NOT NULL
                    )"""
    cursor.execute(create_table_sql)

    # insert
    cursor.executemany(f"INSERT INTO {table_name} VALUES (?, ?)", dictionary.items())

    # commit transaction - writes to the db
    connection.commit()
    # close the cursor
    cursor.close()

    return True


def save_set(
    connection: sqlite3.Connection, table_name: str, rows: set, column: str
) -> bool:
    """
    Save the set rows to the table table_name in the column

    :param connection: db connection of type sqlite.Connection
    :param table_name: name of the table to create
    :param rows: set which will become the the column in the table
    :param column: name of the column in the table
    :return: success status
    """
    cursor = connection.cursor()

    # Drop the table first
    drop_table_sql = f"DROP TABLE IF EXISTS {table_name}"
    cursor.execute(drop_table_sql)

    # Create table
    create_table_sql = f"""CREATE TABLE {table_name} (
                            {column} TEXT PRIMARY KEY
                    )"""
    cursor.execute(create_table_sql)

    # insert
    cursor.executemany(
        f"INSERT INTO {table_name} VALUES (?)",
        zip(
            rows,
        ),
    )

    # commit transaction - writes to the db
    connection.commit()
    # close the cursor
    cursor.close()

    return True


def load_set(connection: sqlite3.Connection, table_name: str, column: str) -> set:
    """
    Retrieve the first column of the table as a set

    :param connection: db connection of type sqlite.Connection
    :param table_name: name of the table to query
    :param column: name of the column in the table to query
    :return: set containing values from the column
    """
    cursor = connection.cursor()
    select_all_query = f"SELECT {column} FROM {table_name}"
    cursor.execute(select_all_query)
    table_as_set = set(map(lambda t: t[0], cursor.fetchall()))
    cursor.close()
    return table_as_set


def load_dict(
    connection: sqlite3.Connection, table_name: str, columns: tuple[str, str]
) -> dict[str, str]:
    """
    Retrieve the values in columns as a name, value pair and create a dict.

    :param connection: db connection of type sqlite.Connection
    :param table_name: name of the table to query
    :param columns: column names in the table
    :return: a dict containing key,value pairs from the columns
    """
    cursor = connection.cursor()
    select_all_query = f"SELECT {columns[0]}, {columns[1]} FROM {table_name}"
    cursor.execute(select_all_query)
    table_as_dict = {k: v for k, v in cursor.fetchall()}
    cursor.close()
    return table_as_dict


def get_user_version(connection: sqlite3.Connection) -> int:
    """
    Retrieve user_version from db

    :connection: sqlite3.Connection: SQLite DB Connection
    """
    query = "PRAGMA user_version"
    cursor = connection.execute(query)
    result = cursor.fetchone()
    version = result[0]
    cursor.close()

    if version:
        return version
    return None


def set_user_version(connection: sqlite3.Connection, version: int):
    """
    Save the version number as user_version in db

    :connection: sqlite3.Connection:
    :version: int: version number to store
    @return:
    """
    query = f"PRAGMA user_version={version}"
    cursor = connection.execute(query)
    # commit transaction - writes to the db
    connection.commit()
    # close the cursor
    cursor.close()


# TODO: combine load_gapped_tables and load_gapped_nuc_tables
def load_gapped_tables(connection: sqlite3.Connection):
    """
    Return a dictionary of gapped sequence tables

    :param connection: db connection of type sqlite.Connection
    :return: dict of gapped sequence tables
    """

    # initialize dictionary to store allele:seq key:value pairs
    gapped_seqs = {}

    for loc in config["loci"]:
        # extend the dictionary with each locus
        gapped_seqs.update(
            load_dict(connection, table_name=f"{loc}_gapped", columns=("allele", "seq"))
        )

    return gapped_seqs


def load_gapped_nuc_tables(connection: sqlite3.Connection):
    """
    Return a dictionary of gapped nucleotide tables

    :param connection: db connection of type sqlite.Connection
    :return: dict of gapped sequence tables
    """

    # initialize dictionary to store allele:seq key:value pairs
    gapped_seqs = {}

    for loc in config["loci"]:
        # extend the dictionary with each locus
        gapped_seqs.update(
            load_dict(
                connection, table_name=f"{loc}_gapped_nuc", columns=("allele", "seq")
            )
        )

    return gapped_seqs


def load_gapped_mature_tables(connection: sqlite3.Connection):
    """
    Return a dictionary of gapped mature sequence tables

    :param connection: db connection of type sqlite.Connection
    :return: dict of gapped mature sequence tables
    """

    # initialize dictionary to store allele:seq key:value pairs
    gapped_seqs = {}

    for loc in config["loci"]:
        # extend the dictionary with each locus
        gapped_seqs.update(
            load_dict(
                connection, table_name=f"{loc}_gapped_mature", columns=("allele", "seq")
            )
        )

    return gapped_seqs


# TODO: combine load_ungapped_tables and load_ungapped_nuc_tables
def load_ungapped_tables(connection: sqlite3.Connection):
    """
    Return a dictionary of ungapped sequence tables

    :param connection: db connection of type sqlite.Connection
    :return: dict of ungapped sequence tables
    """

    # initialize dictionary to store allele:seq key:value pairs
    ungapped_seqs = {}

    for loc in config["loci"]:
        # extend the dictionary with each locus
        ungapped_seqs.update(
            load_dict(
                connection, table_name=f"{loc}_ungapped", columns=("allele", "seq")
            )
        )

    return ungapped_seqs


def load_ungapped_nuc_tables(connection: sqlite3.Connection):
    """
    Return a dictionary of ungapped nucleotide tables

    :param connection: db connection of type sqlite.Connection
    :return: dict of ungapped sequence tables
    """

    # initialize dictionary to store allele:seq key:value pairs
    ungapped_seqs = {}

    for loc in config["loci"]:
        # extend the dictionary with each locus
        ungapped_seqs.update(
            load_dict(
                connection, table_name=f"{loc}_ungapped_nuc", columns=("allele", "seq")
            )
        )

    return ungapped_seqs


def load_ungapped_mature_tables(connection: sqlite3.Connection):
    """
    Return a dictionary of ungapped mature sequence tables

    :param connection: db connection of type sqlite.Connection
    :return: dict of ungapped mature sequence tables
    """

    # initialize dictionary to store allele:seq key:value pairs
    ungapped_seqs = {}

    for loc in config["loci"]:
        # extend the dictionary with each locus
        ungapped_seqs.update(
            load_dict(
                connection,
                table_name=f"{loc}_ungapped_mature",
                columns=("allele", "seq"),
            )
        )

    return ungapped_seqs
