import sys
import tempfile
from urllib.error import URLError

# GitHub URL for IMGT/HLA
IMGT_HLA_URL = "https://raw.githubusercontent.com/ANHIG/IMGTHLA"


def load_sequence_alignment(imgt_version: str):
    """Retrieve sequence alignment from the IMGTHLA GitHub repository and store in a Named Temporary File until processed into the database

    :param imgt_version: The version of the IMGT/HLA database to use
    :return: The Named Temporary File containing the sequence alignment
    """
    msf_p_url = f"{IMGT_HLA_URL}/{imgt_version}/msf/{loc}_prot.msf"
    try:
        
    except URLError as e:
        print(f"Error downloading {msf_p_url}", e, file=sys.stderr)
        sys.exit(1)