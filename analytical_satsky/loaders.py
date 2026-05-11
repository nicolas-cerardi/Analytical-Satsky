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

def load_constellation(*names: str) -> pd.DataFrame:
    """
    Load one or several packaged constellation tables by name.

    Parameters
    ----------
    *names : str
        Name(s) of the constellation table(s) to load. Names are case-insensitive.
        Use `list_constellations()` to see available options.

    Returns
    -------
    pandas.DataFrame
        Constellation table as a dataframe. If several names are provided, the
        corresponding tables are concatenated in the order of the input names.
    """
    if len(names) == 0:
        available = ", ".join(list_constellations())
        raise ValueError(
            f"At least one constellation name must be provided. "
            f"Available options are: {available}"
        )

    normalized_files = {
        key.lower(): filename for key, filename in _CONSTELLATION_FILES.items()
    }

    tables = []

    for name in names:
        name_lower = name.lower()

        try:
            filename = normalized_files[name_lower]
        except KeyError as e:
            available = ", ".join(list_constellations())
            raise ValueError(
                f"Unknown constellation table '{name}'. "
                f"Available options are: {available}"
            ) from e

        tables.append(_load_csv(filename))

    return pd.concat(tables, ignore_index=True)

