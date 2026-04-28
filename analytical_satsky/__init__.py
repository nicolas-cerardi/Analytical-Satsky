"""
Analytical-Satsky: analytical tools to estimate satellite sky occupancy in astronomical observations.
"""

from .loaders import (
    list_constellations,
    load_constellation,
)

from .model import (
    compute_total_satellite_density,
    compute_shell_satellite_density,
    simulate_exposed_time,
    compute_exposure_fraction,
)

from .plot_utils import (
    plot_sky_map,
)

__all__ = [
    "list_constellations",
    "load_constellation",
    "compute_total_satellite_density",
    "compute_shell_satellite_density",
    "simulate_exposed_time",
    "compute_exposure_fraction",
    "plot_sky_map",
]

__version__ = "0.1.0"