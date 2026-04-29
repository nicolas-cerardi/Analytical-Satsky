# API Reference

This page documents the public functions exposed by Analytical-Satsky.

```python
from analytical_satsky import *
```

## Overview

The public API is organized into four categories:

### Constellation management

* `list_constellations()`
* `load_constellation(name)`

### Satellite density modelling

* `compute_total_satellite_density(...)`
* `compute_shell_satellite_density(...)`

### Exposure and observing statistics

* `simulate_exposed_time(...)`
* `compute_exposure_fraction(...)`

### Visualisation

* `plot_sky_map(...)`

---

# Constellation management

## `list_constellations()`

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

## `load_constellation(*names)`

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

---

## `compute_total_satellite_density(...)`

```python
compute_total_satellite_density(
    obsloc,
    shells,
    target_dec,
    target_ra,
    Lfov,
    tobs,
)
```

Compute the expected number of satellites crossing the target field of view during an observation, summed over all input orbital shells.

This is the main high-level modelling function of the package. It evaluates each shell independently and adds their contributions.

### Parameters

* `obsloc` (`astropy.coordinates.EarthLocation`)
  Location of the observer.

* `shells` (`pandas.DataFrame`)
  Table describing the orbital shells of the constellation.

  Required columns:

| Column | Meaning                           | Unit          |
| ------ | --------------------------------- | ------------- |
| `i`    | Orbital inclination               | degrees       |
| `h`    | Orbital altitude                  | km            |
| `n`    | Number of satellites in the shell | dimensionless |

* `target_dec` (`astropy.units.Quantity`)
  Declination of the target. Must be convertible to radians.

* `target_lha` (`astropy.units.Quantity`)
  Local hour angle of the target, defined in the observer frame. Must be convertible to radians.

* `Lfov` (`astropy.units.Quantity`)
  Angular diameter of the telescope field of view. Must be convertible to radians.

* `tobs` (`astropy.units.Quantity`)
  Observation duration. Must be convertible to seconds.

### Returns

* `astropy.units.Quantity`
  Expected number of satellites crossing the target field of view during the observation, summed over all shells.

### Example

```python
import astropy.units as u
from astropy.coordinates import EarthLocation
from analytical_satsky import load_constellation
from analytical_satsky import compute_total_satellite_density

obsloc = EarthLocation.of_site("SKA-Mid")
shells = load_constellation("starlink_filing1")

nsats = compute_total_satellite_density(
    obsloc=obsloc,
    shells=shells,
    target_dec=-20.0 * u.deg,
    target_ra=0.0 * u.deg,
    Lfov=1.0 * u.deg,
    tobs=1.0 * u.hour,
)

print(nsats)
```

### Notes

* Shell tables use plain numeric values (`deg`, `km`, counts), not Astropy quantities.
* For the contribution of a single shell only, use `compute_shell_satellite_density()`.

## `compute_shell_satellite_density(...)`

```python
compute_shell_satellite_density(
    obsloc,
    i,
    nsat,
    hsat,
    target_dec,
    target_lha,
    Lfov,
    tobs,
)
```

Compute the expected number of satellites crossing the target field of view during an observation, for a single orbital shell.

This is a lower-level function than `compute_total_satellite_density()`.

### Parameters

* `obsloc` (`astropy.coordinates.EarthLocation`)
  Location of the observer.

* `i` (`astropy.units.Quantity`)
  Orbital inclination of the shell. Must be convertible to radians.

* `nsat` (`float`)
  Number of satellites in the shell.

* `hsat` (`astropy.units.Quantity`)
  Altitude of the shell above the Earth surface. Must be convertible to metres.

* `target_dec` (`astropy.units.Quantity`)
  Declination of the target. Must be convertible to radians.

* `target_lha` (`astropy.units.Quantity`)
  Local hour angle of the target in the observer frame. Must be convertible to radians.

* `Lfov` (`astropy.units.Quantity`)
  Angular diameter of the telescope field of view. Must be convertible to radians.

* `tobs` (`astropy.units.Quantity`)
  Observation duration. Must be convertible to seconds.

### Returns

* `astropy.units.Quantity`
  Expected number of satellites crossing the target field of view during the observation.

### Example

```python
import astropy.units as u
from astropy.coordinates import EarthLocation
from analytical_satsky import compute_shell_satellite_density

obsloc = EarthLocation.of_site("SKA-Mid")

nsats = compute_shell_satellite_density(
    obsloc=obsloc,
    i=53.0 * u.deg,
    nsat=1584,
    hsat=550.0 * u.km,
    target_dec=30.0 * u.deg,
    target_lha=0.0 * u.deg,
    Lfov=5.0 * u.deg,
    tobs=1.0 * u.hour,
)

print(nsats)
```

### Notes

* This function evaluates one circular orbital shell only.
* To model a full constellation made of several shells, use `compute_total_satellite_density()`.
* `target_lha` is the local hour angle of the target, not its right ascension.
* The returned value is an expected number of satellites, not an integer count from a simulation.

## `simulate_exposed_time(...)`

```python
simulate_exposed_time(
    shells,
    Lfov,
    obsloc,
    target_dec,
    target_lha,
    tobs,
    nstat=100,
)
```

Simulate satellite crossing events during an observation.

This function performs a statistical sampling of satellite passages through the telescope field of view and returns ingress times and crossing durations for all simulated events.

It is useful for estimating time occupancy, data loss fractions, or constructing synthetic timelines of satellite contamination.

### Parameters

* `shells` (`pandas.DataFrame`)
  Table describing the orbital shells of the constellation.

  Required columns:

| Column | Meaning                           | Unit          |
| ------ | --------------------------------- | ------------- |
| `i`    | Orbital inclination               | degrees       |
| `h`    | Orbital altitude                  | km            |
| `n`    | Number of satellites in the shell | dimensionless |

* `Lfov` (`astropy.units.Quantity`)
  Angular diameter of the telescope field of view. Must be convertible to radians.

* `obsloc` (`astropy.coordinates.EarthLocation`)
  Location of the observer.

* `target_dec` (`astropy.units.Quantity`)
  Declination of the target. Must be convertible to radians.

* `target_lha` (`astropy.units.Quantity`)
  Local hour angle of the target. Must be convertible to radians.

* `tobs` (`astropy.units.Quantity`)
  Observation duration. Must be convertible to seconds.

* `nstat` (`int`, optional)
  Number of statistical realisations used in the sampling. Default is `100`.

### Returns

* `all_inits` (`array-like`)
  Simulated ingress times into the effective beam for all sampled satellites from all shells.

* `all_ts` (`array-like`)
  Simulated fly-through durations across the effective beam for all sampled satellites from all shells.

### Example

```python
import astropy.units as u
from astropy.coordinates import EarthLocation
from analytical_satsky import load_constellation
from analytical_satsky import simulate_exposed_time

shells = load_constellation("starlink_march25")

all_inits, all_ts = simulate_exposed_time(
    shells=shells,
    Lfov=5.0 * u.deg,
    obsloc=EarthLocation.of_site("greenwich"),
    target_dec=20.0 * u.deg,
    target_lha=0.0 * u.deg,
    tobs=1.0 * u.hour,
    nstat=200,
)
```

### Notes

* This function returns simulated event times, not deterministic orbital predictions.
* Increasing `nstat` improves statistical robustness at the cost of runtime.
* Outputs can be post-processed to estimate occupancy fractions, contamination timelines, or dead-time statistics.
* For analytical expected counts instead of a sampled catalog, use `compute_total_satellite_density()`. 


## `compute_exposure_fraction(...)`

```python
compute_exposure_fraction(
    tobs,
    ntimestep,
    all_inits,
    all_ts,
)
```

Compute the fraction of observing time during which at least one satellite is present in the effective beam.

This function converts simulated satellite crossing events into time-occupancy fractions for each statistical realisation.

It is typically used after `simulate_exposed_time()`.

### Parameters

* `tobs` (`astropy.units.Quantity`)
  Observation duration. Must be convertible to seconds.

* `ntimestep` (`int`)
  Number of time samples used to discretize the interval `[0, tobs]`.

* `all_inits` (`numpy.ndarray`)
  Ingress times into the effective beam for each statistical realisation.

  Expected shape: `(nstat, nevents)`.

* `all_ts` (`numpy.ndarray`)
  Fly-through durations across the effective beam for each statistical realisation.

  Must have the same shape as `all_inits`.

### Returns

* `numpy.ndarray`
  Exposure fraction for each statistical realisation.

  Returned shape: `(nstat,)`.

Each value lies between `0` and `1`.

### Example

```python
from analytical_satsky import simulate_exposed_time
from analytical_satsky import compute_exposure_fraction

all_inits, all_ts = simulate_exposed_time(...)

fractions = compute_exposure_fraction(
    tobs=1.0 * u.hour,
    ntimestep=1000,
    all_inits=all_inits,
    all_ts=all_ts,
)
```

### Notes

* This function estimates occupancy numerically using a discrete time grid.
* Larger `ntimestep` values provide finer time resolution at the cost of additional runtime.

## `plot_sky_map(...)`

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
* Radial coordinate corresponds to zenith angle (`0°` at zenith, `90°` at horizon).
* Cardinal directions are shown as `N`, `E`, `S`, `W`.
* The color normalization is logarithmic.
* Values below the horizon are masked before plotting.

