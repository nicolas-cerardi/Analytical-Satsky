# Analytical-Satsky

[![Tests](https://github.com/nicolas-cerardi/Analytical-Satsky/actions/workflows/tests.yml/badge.svg)](https://github.com/nicolas-cerardi/Analytical-Satsky/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/nicolas-cerardi/Analytical-Satsky/blob/main/LICENSE)

Modelling satellite constellations seen from astronomical observatories

## Overview

The number of artificial satellites in low earth orbit is increasing very fast, due to the deployment of large communication constellations, comprising from hundreds to tens of thousands of satellites.
This already has in impact on astronomical observations, as optical and radio telescopes respectively see more frequent trails and radio frequency interferences left by these satellites. Still, this is only a glimpse of what will happen in the next decade, when the number of satellites could be multiplied by ten.

Predicting the impact on science from these satellites is a difficult task, and requires to simulate the trajectories of up to 100,000 satellite. Instead of this discrete description of each individual satellite, this package implements an analytical model for predicting the number of satellites, as presented in [Bassa et al. 2022](https://www.aanda.org/articles/aa/abs/2022/01/aa42101-21/aa42101-21.html). It is here adapted to radio telescopes, providing dedicated features

## Features

 - Load predefined satellite shell tables from the major upcoming constellations
 - Compute sky maps of the instantaneous satellite density, from a specific observatory
 - Predict the number of satellites entering the field of view for a given observation

## Installation

Simply git clone and pip install locally

```
git clone https://github.com/nicolas-cerardi/Analytical-Satsky.git
cd Analytical-Satsky 
pip install -e .
```

## Quick start

A minimal working code:

```python
import numpy as np
import astropy.units as u
from analytical_satsky import compute_total_satellite_density, load_constellation

shells_df = load_constellation("starlink_filing2")

n_satellite_in_obs = compute_total_satellite_density(
    -30*u.deg, 
    116.*u.deg, 
    shells_df, 
    np.array([-31.])*u.deg,
    np.array([110.])*u.deg, 
    10.*u.deg, 
    3600*u.s
)

print(n_satellite_in_obs)
```

For more see the demo notebooks:
- Compute a full sky map of instantaneous sky density: [notebooks/satellite_density_maps.ipynb](https://github.com/nicolas-cerardi/Analytical-Satsky/blob/main/notebooks/satellite_density_maps.ipynb)
- Compute satellite occupancy, as a function of target declination, with error bars: [notebooks/compute_satellite_occupancy.ipynb](https://github.com/nicolas-cerardi/Analytical-Satsky/blob/main/notebooks/compute_satellite_occupancy.ipynb)

## Repository structure

analytical_satsky/
notebooks/

## Roadmap

 - In development: sampling realistic satellite trajectories and computing an integral model.

## Citation

Cerardi, Tolley, di Vruno et al., accepted in Astronomy and Astrophysics, 2026

## License
MIT


