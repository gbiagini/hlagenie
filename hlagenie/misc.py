# adapted from py-ard (github.com/nmdp-bioinformatics/py-ard)

import pathlib  # for path manipulation
import tempfile  # for temporary directory access
from .load import load_latest_version  # to get latest IMGTHLA database version


# manipulate input of imgt_version
def get_imgt_version(imgt_version):
    if imgt_version:
        version = imgt_version.replace(".", "")
        if version.isdigit():
            return version
    # if no version specified, use latest
    return load_latest_version()


# define directory to store IMGT database
def get_data_dir(data_dir=None):
    # validate user-input directory
    if data_dir:
        path = pathlib.Path(data_dir)
        if not path.exists() or not path.is_dir():
            raise RuntimeError(f"{data_dir} is not a valid directory")
        data_dir = path
    # use default directory in /tmp if none specified
    else:
        data_dir = pathlib.Path(tempfile.gettempdir()) / "hlagenie"
    return data_dir


# directly from py-ard
def get_imgt_db_versions() -> list[str]:
    """
    Get a list of all available IMGT/HLA database versions

    :return: list of available IMGT/HLA database versions
    """
    import urllib.request
    import json

    req = urllib.request.Request(
        url="https://api.github.com/repos/ANHIG/IMGTHLA/branches?per_page=100"
    )
    res = urllib.request.urlopen(req, timeout=5)
    if res.status == 200:
        json_body = json.loads(res.read())
        versions = list(map(lambda x: x["name"], json_body))
        return versions


# find gap characters in sequence
def find_gaps(sequence: str):
    """
    Identify gap characters in the reference sequence, and return their indices

    :param sequence: reference sequence
    :return: list of gap indices
    """

    # create list of gap indices
    gaps = [i for i, x in enumerate(sequence) if x == "-"]

    return gaps


# get default database directory
def get_default_db_directory():
    return pathlib.Path(tempfile.gettempdir()) / "hlagenie"


def regex_gen(first_ten: str):
    """
    Generate a regex for adjusting to the mature protein sequence

    :param first_ten: first ten amino acids of the mature protein sequence
    :return: regex for adjusting to the mature protein sequence
    """

    # TODO - Comment this code
    regex = ""
    for i in range(10):
        regex += f"{first_ten[i]}[^{first_ten}]*?"
    return regex


def coordinate(sequence: str, regex: str):
    """
    Find the start coordinates of the mature protein sequence

    :param regex: regex for adjusting to the mature protein sequence
    :param sequence: reference sequence
    :return: start coordinates of the mature protein sequence
    """

    import re  # for regex matching

    # find match
    match = re.search(regex, sequence)

    # if match is found, return start coordinate
    start = match.start()

    return start


def coordinate_end(sequence: str, regex: str):
    """
    Find the end coordinates of either the ARD or XRD based on a regular expression

    :param regex: regex for identifying the end of a domain
    :param sequence: reference sequence
    :return: end coordinates of the domain
    """

    import re  # for regex matching

    # find match
    match = re.search(regex, sequence)

    # if match is found, return end coordinate
    end = match.span()[1]

    return end
