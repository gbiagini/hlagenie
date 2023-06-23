

def init(
    imgt_version: str = "Latest",
    data_dir: str = None,
    cache_size: int = DEFAULT_CACHE_SIZE,
    config: dict = None,
):
    from .aamatch import AAMatch
    
    aamatch = AAMatch(
        dbversion = imgt_version,
        data_dir = data_dir,
        cache_size = cache_size,
        config = config
    )
    
    return aamatch