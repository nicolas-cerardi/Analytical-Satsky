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

### Typical workflow

```python
from analytical_satsky import load_constellation
from analytical_satsky import compute_total_satellite_density

shells = load_constellation("starlink_march25")

result = compute_total_satellite_density(
    obsloc=obsloc,
    shells=shells,
    target_ra=target_ra,
    target_dec=target_dec,
    beam_width=beam_width,
    exposure_time=exposure_time,
)
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

*To be documented.*

## `compute_shell_satellite_density(...)`

*To be documented.*

## `simulate_exposed_time(...)`

*To be documented.*

## `ra_to_lha(...)`

*To be documented.*

## `compute_exposure_fraction(...)`

*To be documented.*

## `plot_sky_map(...)`

*To be documented.*
