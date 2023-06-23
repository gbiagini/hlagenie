import sqlite3  # for database operations
import pathlib  # for path manipulation
import misc  # for getting IMGT database versions


def create_db_connection(data_dir, imgt_version):
    """
    Create connection to SQLite database

    :param data_dir: The directory where the database is stored
    :param imgt_version: The version of the IMGT/HLA database to use
    :return: The database connection
    """

    # set database filename
    db_filename = f"{data_dir}/hlagenie-{imgt_version}.db"

    # Check if imgt_version is valid
    if imgt_version != "Latest":
        if not pathlib.Path(db_filename).exists():
            all_versions = misc.get_imgt_db_versions()
            if str(imgt_version) not in all_versions:
                raise ValueError(f"{imgt_version} is not a valid IMGT version")

    # Create the data directory if it doesn't exist
    if not pathlib.Path(data_dir).exists():
        pathlib.Path(data_dir).mkdir(parents=True, exist_ok=True)

    # Create the database file and alert user
    if not pathlib.Path(db_filename).exists():
        print(f"Creating database file {db_filename} as cache")

    # Open the database connection
    file_uri = f"file:{db_filename}"
    return sqlite3.connect(file_uri, uri=True)
