import pandas as pd
from importlib import resources

_CONSTELLATION_FILES = {
    "starlink_march25": "starlink_march25.csv",
    "starlink_scaled40000": "starlink_scaled40000.csv",
    "oneweb": "oneweb.csv",
    "qianfan": "qianfan.csv",
    "guowang": "guowang.csv",
    "starlink_filing1": "starlink_filing1.csv",
    "starlink_filing2": "starlink_filing2.csv",
    "leo": "leo.csv",
}

def _load_csv(filename: str) -> pd.DataFrame:
    with resources.files("analytical_satsky.shells").joinpath(filename).open("rb") as f:
        return pd.read_csv(f, index_col=0)

def list_constellations() -> list[str]:
    """Return the list of available constellation tables."""
    return sorted(_CONSTELLATION_FILES)

def load_constellation(name: str) -> pd.DataFrame:
    """
    Load a packaged constellation table by name.
    Parameters
    ----------
    name : str
        Name of the constellation table. Use `list_constellations()` to see available options.

    Returns
    -------
    pandas.DataFrame
        Constellation table as a dataframe.
    """
    try:
        filename = _CONSTELLATION_FILES[name]
    except KeyError as e:
        available = ", ".join(list_constellations())
        raise ValueError(
            f"Unknown constellation table '{name}'. Available options are: {available}"
        ) from e

    return _load_csv(filename)

