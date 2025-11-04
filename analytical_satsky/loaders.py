import pandas as pd
from importlib import resources

def _load_csv(name: str) -> pd.DataFrame:
    """Helper function to load a CSV file from the data folder."""
    with resources.files("analytical_satsky.shells").joinpath(name).open("rb") as f:
        return pd.read_csv(f, index_col=0)

def load_table_starlink_march25():
    return _load_csv("starlink_march25.csv")

def load_table_starlink_scaled40000():
    return _load_csv("starlink_scaled40000.csv")

def load_table_oneweb():
    return _load_csv("oneweb.csv")

def load_table_qianfan():
    return _load_csv("qianfan.csv")

def load_table_guowang():
    return _load_csv("guowang.csv")
