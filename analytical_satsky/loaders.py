import pandas as pd
from importlib import resources

def _load_csv(name: str) -> pd.DataFrame:
    """Helper function to load a CSV file from the data folder."""
    with resources.files("analytical_satsky.shells").joinpath(name).open("rb") as f:
        return pd.read_csv(f, index_col=0)

def load_table_starlink_march25():
    return _load_csv("starlink_march25.csv")
