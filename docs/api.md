# API Reference

This page documents the public functions exposed by Analytical-Satsky.

```python id="rmwdqv"
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
* `ra_to_lha(...)`
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

## `load_constellation(name)`

```python
load_constellation(name: str) -> pandas.DataFrame
```

Load a predefined satellite constellation table bundled with the package.

The returned object is a `pandas.DataFrame` where each row represents one orbital shell of the constellation.

### Parameters

* `name` (`str`)
  Name of the constellation to load.

Use `list_constellations()` to see all available built-in options.

### Returns

* `pandas.DataFrame`
  A table containing one row per orbital shell, with the columns:

| Column | Meaning                           | Unit          |
| ------ | --------------------------------- | ------------- |
| `i`    | Orbital inclination               | degrees       |
| `h`    | Orbital altitude                  | km            |
| `n`    | Number of satellites in the shell | dimensionless |

### Example

```python
from analytical_satsky import load_constellation

shells = load_constellation("oneweb")
print(shells.head())
```

Example output:

```python
      i      h      n
0  87.9  1200.0   720
1  40.0  1200.0   648
...
```

### Raises

* `ValueError`
  If `name` does not match an available packaged constellation.

### Notes

* Use `list_constellations()` to discover valid names.
* The exact available constellations may evolve between package versions.
* Returned values are plain numeric columns (not Astropy quantities).
* Users may also create custom constellation tables manually using a compatible `pandas.DataFrame`.

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

*To be documented.*

## `ra_to_lha(...)`

*To be documented.*

## `compute_exposure_fraction(...)`

*To be documented.*

## `plot_sky_map(...)`

*To be documented.*
