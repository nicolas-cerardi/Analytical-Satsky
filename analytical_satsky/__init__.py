"""
Analytical-Satsky: analytical tools to estimate satellite sky occupancy and exposure.
"""

from .loaders import (
    load_table_starlink_march25,
    load_table_oneweb_march25,
    load_table_leo,
)

from .model import (
    compute_total_satellite_density,
    compute_shell_density,
    compute_exposed_time,
)

from .plot_utils import (
    plot_sky_map,
)

__all__ = [
    "load_table_starlink_march25",
    "load_table_oneweb_march25",
    "load_table_leo",
    "compute_total_satellite_density",
    "compute_shell_density",
    "compute_exposed_time",
    "plot_sky_map",
]

__version__ = "0.1.0"