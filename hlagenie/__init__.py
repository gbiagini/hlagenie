from .configs import config  # for configurations

__author__ = "Giovanni Biagini"
__version__ = "0.1.0"


def init(
    imgt_version: str = "Latest",
    data_dir: str = None,
    cache_size: int = config["DEFAULT_CACHE_SIZE"],
    ungap: bool = True,
):
    from .genie import GENIE

    genie = GENIE(
        imgt_version=imgt_version,
        data_dir=data_dir,
        cache_size=cache_size,
        ungap=ungap,
    )

    return genie
