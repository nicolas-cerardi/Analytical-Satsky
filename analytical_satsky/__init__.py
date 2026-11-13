"""
Analytical-Satsky: analytical tools to estimate satellite sky occupancy in astronomical observations.
"""

from .loaders import (
    list_constellations,
    load_constellation,
)

from .model import (
    compute_occupancy_fraction,
)

from .plot_utils import (
    plot_sky_map,
)

from .shell_obs import (
    SingleShellObs,
    MultiShellObs,
    BasicObservation,
    SingleShellFoV,
    MultiShellFoV,
    SingleShellFlux,
    MultiShellFlux,
    IntegralObsModel
)

__all__ = [
    "list_constellations",
    "load_constellation",
    "compute_occupancy_fraction",
    "plot_sky_map",
    "SingleShellObs",
    "MultiShellObs",
    "BasicObservation",
    "SingleShellFoV",
    "MultiShellFoV",
    "SingleShellFlux",
    "MultiShellFlux",
    "IntegralObsModel"
]

__version__ = "0.3.0"