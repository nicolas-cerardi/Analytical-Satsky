# Constellations

Analytical-Satsky represents a satellite constellation as a collection of orbital shells.

Each shell is described by three parameters:

| Column | Meaning                           | Unit          |
| ------ | --------------------------------- | ------------- |
| `i`    | Orbital inclination               | degrees       |
| `h`    | Orbital altitude                  | km            |
| `n`    | Number of satellites in the shell | dimensionless |

The constellation data are stored as CSV files and loaded as `pandas.DataFrame` objects.

## Available constellations

The package currently includes the following predefined constellations:

| Name                   | File                       |
| ---------------------- | -------------------------- |
| `starlink_filing1`     | `starlink_filing1.csv`     |
| `starlink_filing2`     | `starlink_filing2.csv`     |
| `oneweb`               | `oneweb.csv`               |
| `qianfan`              | `qianfan.csv`              |
| `guowang`              | `guowang.csv`              |
| `leo`                  | `leo.csv`                  |
| `starlink_march25`     | `starlink_march25.csv`     |
| `starlink_scaled40000` | `starlink_scaled40000.csv` |

These can be loaded with:

```python
from analytical_satsky import load_constellation

shells = load_constellation("leo")
```

## Example: the `leo` constellation

A constellation is a table where each row corresponds to one orbital shell.

For example, the `leo` constellation can be inspected with:

```python
from analytical_satsky import load_constellation

leo = load_constellation("leo")
print(leo)
```

This returns a table with the following structure:

|   i |   h |   n |
| --: | --: | --: |
| 51.9 | 630 | 1156 |
| 42.0 | 610 | 1296 |
| 33.0 | 590 | 784 |

where:

* `i` is the inclination of orbits in the shell, in degrees
* `h` is the shell altitude, in kilometres
* `n` is the number of satellites in that shell

## Creating a custom constellation

Users can define their own constellation by creating a `pandas.DataFrame` with the same three columns:

```python
import pandas as pd

custom_constellation = pd.DataFrame({
    "i": [53.0, 70.0, 97.6],
    "h": [550.0, 570.0, 1200.0],
    "n": [1584, 720, 300]
})
```

This custom constellation can then be passed directly to the model routines:

```python
from analytical_satsky import compute_total_satellite_density

result = compute_total_satellite_density(
    obsloc=obsloc,
    shells=custom_constellation,
    target_ra=target_ra,
    target_dec=target_dec,
    beam_width=beam_width,
    exposure_time=exposure_time,
)
```

## Unit conventions

The constellation table should not contain Astropy quantities.

Expected units:

* `i`: degrees
* `h`: kilometres
* `n`: number of satellites

The package expects the shell parameters in these units and transforms them later in Astropy quantities.

## Notes

The model assumes that each shell is spherical (with polar caps cut off) and uniformly populated by `n` satellites. More complex orbital populations can be approximated by splitting a shell into several orbital planes with different longitude of the ascending node, $\Omega$.
