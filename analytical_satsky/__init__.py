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
    compute_occupancy_fraction,
)

from .plot_utils import (
    plot_sky_map,
)

from .shell_obs import (
    SingleShellObs,
    MultiShellObs,
)

__all__ = [
    "list_constellations",
    "load_constellation",
    "compute_total_satellite_density",
    "compute_shell_satellite_density",
    "simulate_exposed_time",
    "compute_exposure_fraction",
    "compute_occupancy_fraction",
    "plot_sky_map",
    "SingleShellObs",
    "MultiShellObs",
]

__version__ = "0.2.0"