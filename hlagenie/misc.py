# adapted from py-ard (github.com/nmdp-bioinformatics/py-ard)

import pathlib  # for path manipulation
import tempfile  # for temporary directory access


# manipulate input of imgt_version
def get_imgt_version(imgt_version):
    if imgt_version:
        version = imgt_version.replace(".", "")
        if version.isdigit():
            return version
    # if no version specified, use latest
    return "Latest"


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
def get_imgt_db_versions() -> List[str]:
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
