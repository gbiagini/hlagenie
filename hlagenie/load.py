import sys
import tempfile
import requests
from Bio import AlignIO
from urllib.error import URLError


def load_sequence_alignment(imgt_version: str, loc: str, imputed: bool):
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
    else:
        # GitHub URL for IMGT/HLA
        IMGT_HLA_URL = "https://raw.githubusercontent.com/ANHIG/IMGTHLA"

    msf_p_url = f"{IMGT_HLA_URL}/{imgt_version}/msf/{loc}_prot.msf"
    try:
        # download the file data from the IMGT_HLA GitHub repository
        request = requests.get(msf_p_url, timeout=15)

        # create a Named Temporary File to store the file data
        tf = tempfile.NamedTemporaryFile()

        # write the file data to the Named Temporary File
        tf.write(request.content)

        # read the data into a multiple sequence alignment object
        multi_seq = AlignIO.read(tf.name, "msf")

        # close the Named Temporary File, deleting it
        tf.close()

    except URLError as e:
        print(f"Error downloading {msf_p_url}", e, file=sys.stderr)
        sys.exit(1)

    # return the multiple sequence alignment object
    return multi_seq


def load_nucleotide_alignment(imgt_version: str, loc: str, imputed: bool):
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
    else:
        # GitHub URL for IMGT/HLA
        IMGT_HLA_URL = "https://raw.githubusercontent.com/ANHIG/IMGTHLA"

    msf_n_url = f"{IMGT_HLA_URL}/{imgt_version}/msf/{loc}_nuc.msf"
    try:
        # download the file data from the IMGT_HLA GitHub repository
        request = requests.get(msf_n_url, timeout=15)

        # create a Named Temporary File to store the file data
        tf = tempfile.NamedTemporaryFile()

        # write the file data to the Named Temporary File
        tf.write(request.content)

        # read the data into a multiple sequence alignment object
        multi_seq = AlignIO.read(tf.name, "msf")

        # close the Named Temporary File, deleting it
        tf.close()

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
