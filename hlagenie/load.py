import sys
import os
import csv
import tempfile
import requests
from Bio import AlignIO
from urllib.error import URLError

# Base URL for the IMGT/HLA GitHub repository (same source py-ard reads).
IMGT_HLA_URL = "https://raw.githubusercontent.com/ANHIG/IMGTHLA"


def load_sequence_alignment(
    imgt_version: str, loc: str, imputed: bool, imputation_method: str
):
    """Retrieve sequence alignment from the IMGTHLA GitHub repository and store in a Named Temporary File until processed into the database

    :param imgt_version: The version of the IMGT/HLA database to use
    :param loc: The HLA locus to retrieve the sequence alignment for
    :return: Bio.Align.MultipleSeqAlignment object
    """
    # for DRB3, DRB4, and DRB5, use DRB345 alignment
    if loc in ["DRB3", "DRB4", "DRB5"]:
        loc = "DRB345"

    # if imputed is True, use the imputed sequence alignment
    if imputed:
        # GitHub URL for imputed sequences
        IMGT_HLA_URL = (
            "https://raw.githubusercontent.com/gbiagini/hla-imputed-sequences"
        )
        msf_p_url = (
            f"{IMGT_HLA_URL}/{imgt_version}/{imputation_method}/msf/{loc}_prot.msf"
        )
    else:
        # GitHub URL for IMGT/HLA
        IMGT_HLA_URL = "https://raw.githubusercontent.com/ANHIG/IMGTHLA"
        msf_p_url = f"{IMGT_HLA_URL}/{imgt_version}/msf/{loc}_prot.msf"

    try:
        # download the file data from the IMGT_HLA GitHub repository
        request = requests.get(msf_p_url, timeout=15)

        # create a Named Temporary File to store the file data
        tf = tempfile.NamedTemporaryFile(delete=False)

        # write the file data to the Named Temporary File
        tf.write(request.content)

        # read the data into a multiple sequence alignment object
        multi_seq = AlignIO.read(tf.name, "msf")

        # close the Named Temporary File, deleting it
        tf.close()
        os.unlink(tf.name)

    except URLError as e:
        print(f"Error downloading {msf_p_url}", e, file=sys.stderr)
        sys.exit(1)

    # return the multiple sequence alignment object
    return multi_seq


def load_nucleotide_alignment(
    imgt_version: str, loc: str, imputed: bool, imputation_method: str
):
    """Retrieve nucleotide alignment from the IMGTHLA GitHub repository and store in a Named Temporary File until processed into the database

    :param imgt_version: The version of the IMGT/HLA database to use
    :param loc: The HLA locus to retrieve the sequence alignment for
    :return: Bio.Align.MultipleSeqAlignment object
    """

    # for DRB3, DRB4, and DRB5, use DRB345 alignment
    if loc in ["DRB3", "DRB4", "DRB5"]:
        loc = "DRB345"

    # if imputed is True, use the imputed sequence alignment
    if imputed:
        # GitHub URL for imputed sequences
        IMGT_HLA_URL = (
            "https://raw.githubusercontent.com/gbiagini/hla-imputed-sequences"
        )
        msf_n_url = (
            f"{IMGT_HLA_URL}/{imgt_version}/{imputation_method}/msf/{loc}_nuc.msf"
        )
    else:
        # GitHub URL for IMGT/HLA
        IMGT_HLA_URL = "https://raw.githubusercontent.com/ANHIG/IMGTHLA"
        msf_n_url = f"{IMGT_HLA_URL}/{imgt_version}/msf/{loc}_nuc.msf"

    try:
        # download the file data from the IMGT_HLA GitHub repository
        request = requests.get(msf_n_url, timeout=15)

        # create a Named Temporary File to store the file data
        tf = tempfile.NamedTemporaryFile(delete=False)

        # write the file data to the Named Temporary File
        tf.write(request.content)

        # read the data into a multiple sequence alignment object
        multi_seq = AlignIO.read(tf.name, "msf")

        # close the Named Temporary File, deleting it
        tf.close()
        os.unlink(tf.name)

    except URLError as e:
        print(f"Error downloading {msf_n_url}", e, file=sys.stderr)
        sys.exit(1)

    # return the multiple sequence alignment object
    return multi_seq


def load_latest_version():
    """From py-ard. Get latest version of the IMGT/HLA database

    :return: latest version of the IMGT/HLA database
    :rtype: str
    """
    from urllib.request import urlopen

    version_txt = (
        "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/release_version.txt"
    )
    try:
        response = urlopen(version_txt, timeout=15)
    except URLError as e:
        print(f"Error downloading {version_txt}", e, file=sys.stderr)
        sys.exit(1)

    version = 0
    for line in response:
        l = line.decode("utf-8")
        if l.find("version:") != -1:
            # Version line looks like
            # # version: IPD-IMGT/HLA 3.51.0
            version = l.split()[-1].replace(".", "")
    return version


# The following loaders read the same IMGT/HLA nomenclature files that py-ard
# uses to build its allele reduction tables. They let HLAGenie reduce allele
# names to two-field form without depending on py-ard at runtime.


def _fetch_lines(url: str):
    """Download a text file and return its stripped lines."""
    try:
        response = requests.get(url, timeout=15)
    except URLError as e:
        print(f"Error downloading {url}", e, file=sys.stderr)
        sys.exit(1)
    return [line.strip() for line in response.text.splitlines()]


def load_allele_names(imgt_version: str):
    """Return the list of allele names from ``Allelelist.<version>.txt``.

    :param imgt_version: IMGT/HLA database version
    :return: list of allele names (e.g. ``A*01:01:01:01``)
    """
    if imgt_version == "3130":
        # 3130 was renamed to 3131 for the Allelelist file only (per py-ard).
        imgt_version = "3131"
    url = f"{IMGT_HLA_URL}/Latest/allelelist/Allelelist.{imgt_version}.txt"
    lines = _fetch_lines(url)
    # Skip the first 6 header lines; the 7th line is the CSV header.
    reader = csv.DictReader(lines[6:])
    return [row["Allele"] for row in reader]


def load_g_group_alleles(imgt_version: str):
    """Return full allele names appearing in ``hla_nom_g.txt``.

    :param imgt_version: IMGT/HLA database version
    :return: list of full allele names (locus + allele)
    """
    url = f"{IMGT_HLA_URL}/{imgt_version}/wmda/hla_nom_g.txt"
    alleles = []
    for line in _fetch_lines(url)[6:]:
        if not line:
            continue
        fields = line.split(";")
        if len(fields) >= 3 and fields[1] and fields[2]:
            locus, allele_list = fields[0], fields[1]
            alleles.extend(locus + a for a in allele_list.split("/"))
    return alleles


def load_p_group_pairs(imgt_version: str):
    """Return ``(full_allele, full_p_group_name)`` pairs from ``hla_nom_p.txt``.

    :param imgt_version: IMGT/HLA database version
    :return: list of ``(allele, p_group)`` tuples
    """
    url = f"{IMGT_HLA_URL}/{imgt_version}/wmda/hla_nom_p.txt"
    pairs = []
    for line in _fetch_lines(url)[6:]:
        if not line:
            continue
        fields = line.split(";")
        if len(fields) >= 3 and fields[1] and fields[2]:
            locus, allele_list, p_group = fields[0], fields[1], fields[2]
            full_p = locus + p_group
            pairs.extend((locus + a, full_p) for a in allele_list.split("/"))
    return pairs
