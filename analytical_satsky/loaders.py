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
    """
    Return the list of predefined satellite constellation tables bundled with
    the package.

    This is the recommended way to discover which built-in constellation
    datasets are available on the current installed version.

    Returns
    -------
    list[str]
        A list of constellation names that can be passed to
        ``load_constellation()``.

    Examples
    --------
    >>> from analytical_satsky import list_constellations
    >>> list_constellations()
    ['guowang', 'leo', 'oneweb', 'qianfan', 'starlink_filing1', 'starlink_filing2', 'starlink_march25', 'starlink_scaled40000']

    Notes
    -----
    The returned names depend on the installed package version. Use
    ``load_constellation(name)`` to load one of these datasets as a
    ``pandas.DataFrame``.
    """
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
        Each row describes one orbital shell, with columns:

        - ``i`` : orbital inclination, in degrees
        - ``h`` : orbital altitude, in km
        - ``n`` : number of satellites in the shell, dimensionless

    Raises
    ------
    ValueError
        If no name is provided, or if one of the input names does not match
        an available packaged constellation.

    Examples
    --------
    Load a single constellation:

    >>> from analytical_satsky import load_constellation
    >>> shells = load_constellation("oneweb")

    Load several constellations and combine them:

    >>> shells = load_constellation("leo", "qianfan", "oneweb")

    Notes
    -----
    Returned values are plain numeric columns, not Astropy quantities. Users
    may also create custom constellation tables manually using a compatible
    ``pandas.DataFrame``.
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

