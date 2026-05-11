# API Reference

This page documents the public functions exposed by Analytical-Satsky.

```python
from analytical_satsky import *
```

## Overview

The public API is organized into four categories:

### Constellation management

* [`list_constellations`](#list_constellations)
* [`load_constellation`](#load_constellation)

### Satellite density and occupancy modelling

* [`SingleShellObs`](#singleshellobs)
* [`MultiShellObs`](#multishellobs)
* [`compute_occupancy_fraction`](#compute_occupancy_fraction)

### Visualisation

* [`plot_sky_map`](#plot_sky_map)

---

# Constellation management

## `list_constellations`

```python
list_constellations() -> list[str]
```

Return the list of predefined satellite constellation tables bundled with the package.

This is the recommended way to discover which built-in constellation datasets are available on the current installed version.

### Returns

* `list[str]`
  A list of constellation names that can be passed to `load_constellation()`.

### Example

```python
from analytical_satsky import list_constellations

names = list_constellations()
print(names)
```

Example output:

```python
[
    "starlink_march25",
    "starlink_scaled40000",
    "oneweb",
    "qianfan",
    "guowang",
    "starlink_filing1",
    "starlink_filing2",
    "leo",
]
```

### Notes

* The returned names depend on the installed package version.
* Use `load_constellation(name)` to load one of these datasets as a `pandas.DataFrame`.
* For more on constellations see the [constellation page](constellations.md).

## `load_constellation`

```python
load_constellation(*names: str) -> pandas.DataFrame
```

Load one or several predefined satellite constellation tables bundled with the package.

The returned object is a `pandas.DataFrame` where each row represents one orbital shell of a constellation.

If several names are provided, the corresponding tables are concatenated in the order of the input names.

### Parameters

* `*names` (`str`)
  One or several constellation names to load.

Names are case-insensitive.

Use `list_constellations()` to see all available built-in options.

### Returns

* `pandas.DataFrame`
  A table containing one row per orbital shell, with the columns:

| Column | Meaning                           | Unit          |
| ------ | --------------------------------- | ------------- |
| `i`    | Orbital inclination               | degrees       |
| `h`    | Orbital altitude                  | km            |
| `n`    | Number of satellites in the shell | dimensionless |

### Examples

Load a single constellation:

```python
from analytical_satsky import load_constellation

shells = load_constellation("oneweb")
```

Load several constellations and combine them:

```python
shells = load_constellation("leo", "qianfan", "oneweb")
```

### Raises

* `ValueError`
  If no name is provided, or if one of the input names does not match an available packaged constellation.

### Notes

* Use `list_constellations()` to show valid names.
* Returned values are plain numeric columns (not Astropy quantities).
* Users may also create custom constellation tables manually using a compatible `pandas.DataFrame`.
* When several names are given, rows are appended in the same order as the input arguments.
* For more on constellations see the [constellation page](constellations.md).

---
## `SingleShellObs`

Cached model for the contribution of a single satellite shell.

Use this class to evaluate the geometry, projected satellite density, apparent motion, and expected number of satellite crossings for one orbital shell.

```python
SingleShellObs(
    obsloc,
    shell,
    target_dec,
    target_lha,
    Lfov,
    tobs,
)
```

### Parameters

* `obsloc` (`astropy.coordinates.EarthLocation`)
  Location of the observer.

* `shell` (`pandas.Series`)
  Row describing one orbital shell.

  Required fields:

| Field | Meaning                           | Unit          |
| ----- | --------------------------------- | ------------- |
| `i`   | Orbital inclination               | degrees       |
| `h`   | Orbital altitude                  | km            |
| `n`   | Number of satellites in the shell | dimensionless |

* `target_dec` (`astropy.units.Quantity`)
  Declination of the target line of sight. Must be angular.

* `target_lha` (`astropy.units.Quantity`)
  Local hour angle of the target line of sight, defined in the observer frame. Must be angular.

* `Lfov` (`astropy.units.Quantity`)
  Effective angular diameter of the telescope field of view.

* `tobs` (`astropy.units.Quantity`)
  Observation duration. Must be convertible to seconds.

### Useful attributes

* `rho_sat`
  Projected satellite density for the shell.

* `wsat`
  Apparent satellite angular velocity.

* `nsats`
  Expected number of satellites crossing the field of view during the observation.

### Example

```python
import astropy.units as u
from astropy.coordinates import EarthLocation
from analytical_satsky import load_constellation, SingleShellObs

obsloc = EarthLocation.of_site("SKA-Mid")
shells = load_constellation("starlink_filing1")

obs = SingleShellObs(
    obsloc=obsloc,
    shell=shells.iloc[0],
    target_dec=-20.0 * u.deg,
    target_lha=0.0 * u.deg,
    Lfov=1.0 * u.deg,
    tobs=1.0 * u.hour,
)

print(obs.nsats)
```

### Notes

* Shell tables use plain numeric values (`deg`, `km`, counts), not Astropy quantities.
* Derived quantities are evaluated lazily and cached after first access.
* `target_lha` is the local hour angle of the target, not its right ascension.
* For information on the model assumptions see the [model page](model.md).

---

## `MultiShellObs`

Cached model for a full constellation made of several orbital shells.

This class builds one `SingleShellObs` model per shell and combines their contributions.

```python
MultiShellObs(
    obsloc,
    shells_df,
    target_dec,
    target_lha,
    Lfov,
    tobs,
)
```

### Parameters

* `obsloc` (`astropy.coordinates.EarthLocation`)
  Location of the observer.

* `shells_df` (`pandas.DataFrame`)
  Table describing the orbital shells of the constellation.

  Required columns:

| Column | Meaning                           | Unit          |
| ------ | --------------------------------- | ------------- |
| `i`    | Orbital inclination               | degrees       |
| `h`    | Orbital altitude                  | km            |
| `n`    | Number of satellites in the shell | dimensionless |

* `target_dec` (`astropy.units.Quantity`)
  Declination of the target line of sight. Must be angular.

* `target_lha` (`astropy.units.Quantity`)
  Local hour angle of the target line of sight, defined in the observer frame. Must be angular.

* `Lfov` (`astropy.units.Quantity`)
  Effective angular diameter of the telescope field of view.

* `tobs` (`astropy.units.Quantity`)
  Observation duration. Must be convertible to seconds.

### Useful attributes and methods

* `shell_models`
  List of `SingleShellObs` models, one per shell.

* `nsats_per_shell`
  Expected number of satellite crossings for each shell.

* `total_satellite_density`
  Sum of projected satellite densities over all shells.

* `total_nsats`
  Total expected number of satellites crossing the field of view.

* `sample_passes(nstat)`
  Draw stochastic satellite passes from the model.

### Example

```python
import astropy.units as u
from astropy.coordinates import EarthLocation
from analytical_satsky import load_constellation, MultiShellObs

obsloc = EarthLocation.of_site("SKA-Mid")
shells = load_constellation("starlink_filing1")

obs = MultiShellObs(
    obsloc=obsloc,
    shells_df=shells,
    target_dec=-20.0 * u.deg,
    target_lha=0.0 * u.deg,
    Lfov=1.0 * u.deg,
    tobs=1.0 * u.hour,
)

print(obs.total_nsats)
print(obs.nsats_per_shell)
```

### Notes

* Use `SingleShellObs` for one shell and `MultiShellObs` for a full constellation.
* `total_nsats` replaces the older high-level density-count interface.
* `target_lha` is the local hour angle of the target, not its right ascension.
* The returned values are expected satellite counts, not integer counts from a deterministic simulation.


## `compute_occupancy_fraction`

```python
compute_occupancy_fraction(
    tobs,
    ntimestep,
    all_inits,
    all_ts,
)
```

Compute the fraction of observing time during which at least one satellite is present in the effective beam.

This function converts sampled satellite crossing events into time-occupancy fractions for each statistical realisation.

It is typically used after generating satellite passes with `MultiShellObs.sample_passes(...)`.

### Parameters

* `tobs` (`astropy.units.Quantity`)
  Observation duration. Must be convertible to seconds.

* `ntimestep` (`int`)
  Number of time samples used to discretize the interval `[0, tobs]`.

* `all_inits` (`numpy.ndarray`)
  Ingress times into the effective beam for each statistical realisation.

  Expected shape: `(nstat, nevents)`.

  Values are assumed to be expressed in seconds.

* `all_ts` (`numpy.ndarray`)
  Fly-through durations across the effective beam for each statistical realisation.

  Must have the same shape as `all_inits`.

  Values are assumed to be expressed in seconds.

  Entries with non-positive duration are ignored. This allows the use of
  placeholder events when the number of sampled events varies between
  statistical realisations.

### Returns

* `numpy.ndarray`
  Exposure fraction for each statistical realisation.

  Returned shape: `(nstat,)`.

### Example

```python
import astropy.units as u
import numpy as np
import pandas as pd

from astropy.coordinates import EarthLocation

from analytical_satsky import load_constellations, MultiShellObs, compute_occupancy_fraction

obsloc = EarthLocation.of_site("SKA-Mid")
shells = load_constellation("starlink_filing1")
target_dec = np.array([30.0]) * u.deg
target_lha = np.array([0.0]) * u.deg

obs = MultiShellObs(obsloc, shells, target_dec, target_lha, 10.0 * u.deg, 3600 * u.s)
all_inits, all_ts = obs.sample_passes(nstat=100)

fractions = compute_occupancy_fraction(
    tobs=tobs,
    ntimestep=3600,
    all_inits=all_inits,
    all_ts=all_ts,
)

print(fractions.mean())
```

### Notes

* This function estimates occupancy numerically using a discrete time grid.
* Larger `ntimestep` values provide finer time resolution at the cost of additional runtime.
* The implementation uses a difference-array / cumulative-sum approach to avoid explicitly constructing a boolean occupancy mask for every satellite pass and every timestep.
* `compute_occupancy_fraction()` replaces the older `compute_exposure_fraction()` interface.

---

# Visualisation

## `plot_sky_map`

```python
plot_sky_map(
    sky_map,
    obsloc,
    target_lha,
    target_dec,
    cmap="magma",
    vmin=10.0,
    vmax=10000.0,
    return_fig=False,
)
```

Plot a two-dimensional sky map in local horizon coordinates using a zenith-centred polar projection.

The input map is defined on a local hour angle / declination grid and is transformed to altitude / azimuth coordinates for the selected observatory. Points below the horizon are automatically masked.

### Parameters

* `sky_map` (`astropy.units.Quantity`)
  Two-dimensional map to display. Must have the same shape as the grid defined by `target_lha` and `target_dec`.

* `obsloc` (`astropy.coordinates.EarthLocation`)
  Location of the observer.

* `target_lha` (`astropy.units.Quantity`)
  One-dimensional local hour angle grid. Must be convertible to radians.

* `target_dec` (`astropy.units.Quantity`)
  One-dimensional declination grid. Must be convertible to radians.

* `cmap` (`str`, optional)
  Matplotlib colormap used for the plot. Default is `"magma"`.

* `vmin` (`float`, optional)
  Minimum value of the logarithmic color scale. Default is `10.0`.

* `vmax` (`float`, optional)
  Maximum value of the logarithmic color scale. Default is `10000.0`.

* `return_fig` (`bool`, optional)
  If `True`, return the Matplotlib figure and axes objects. Default is `False`.

### Returns

* `tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]` or `None`
  Returns `(fig, ax)` if `return_fig=True`, otherwise returns `None`.

### Example

```python
import astropy.units as u
from astropy.coordinates import EarthLocation
from analytical_satsky import plot_sky_map

plot_sky_map(
    sky_map=sky_map,
    obsloc=EarthLocation.of_site("greenwich"),
    target_lha=lha_grid * u.deg,
    target_dec=dec_grid * u.deg,
)
```

### Notes

* The plot uses a zenith-centred polar projection.
* Cardinal directions are shown as `N`, `E`, `S`, `W`.
* Values below the horizon are masked before plotting.
