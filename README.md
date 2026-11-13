# Analytical-Satsky

[![Tests](https://github.com/nicolas-cerardi/Analytical-Satsky/actions/workflows/tests.yml/badge.svg)](https://github.com/nicolas-cerardi/Analytical-Satsky/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/nicolas-cerardi/Analytical-Satsky/blob/main/LICENSE)
[![astro-ph.IM](https://img.shields.io/badge/astro--ph.IM-arXiv%3A2604.22694-B31B1B.svg)](https://arxiv.org/abs/2604.22694)

Modelling satellite constellations seen from astronomical observatories

<p align="center">
  <img src="notebooks/sky_map_skamid.png" alt="Sky map from SKA-Mid" width="400">
</p>

## Overview

The number of artificial satellites in low earth orbit is increasing very fast, due to the deployment of large communication constellations, comprising from hundreds to tens of thousands of satellites.
This already has an impact on astronomical observations, as optical and radio telescopes respectively see more frequent trails and radio frequency interferences left by these satellites. Still, this is only a glimpse of what will happen in the next decade, when the number of satellites could be multiplied by ten.

Predicting the impact on science from these satellites is a difficult task, and requires to simulate the trajectories of up to 100,000 satellites. Instead of this discrete description of each individual satellite, this package implements an analytical model for predicting the number of satellites, as presented in [Bassa et al. 2022](https://www.aanda.org/articles/aa/abs/2022/01/aa42101-21/aa42101-21.html). It is here adapted to radio telescopes, providing dedicated features.

## Features

- Load predefined satellite shell tables from the major upcoming constellations
- Implements a analytical modelling of satellite constellations, specifically adapted to radio telescopes' large fields of view.
- Compute sky maps of the instantaneous satellite density, from a specific observatory
- Predict the number of satellites entering the field of view for a given observation
- Estimate the fraction of time with at least one satellite in the field of view.
- Samples mock satellite orbits from the analytical density. 

## Installation

Simply git clone and pip install locally

```bash
git clone https://github.com/nicolas-cerardi/Analytical-Satsky.git
cd Analytical-Satsky 
pip install -e .
```

## Quick start

A minimal working code:

```python
import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from analytical_satsky import load_constellation, MultiShellObs

shells_df = load_constellation("starlink_filing2")
obsloc = EarthLocation.of_site("SKA-Low")

multi_shell_obs = MultiShellObs(
    obsloc, 
    shells_df, 
    np.array([-31.])*u.deg, 
    np.array([0.])*u.deg, 
    10.*u.deg, 
    3600*u.s
)

print(multi_shell_obs.total_nsats)
```

For large fields of view, the integral model samples a mock catalogue of the satellites crossing the field of view, from a true RA/Dec pointing and an explicit observation time:

```python
from astropy.time import Time
from analytical_satsky import BasicObservation, IntegralObsModel

obs = BasicObservation(obsloc, ra=27.0*u.deg, dec=-30.0*u.deg)

integral_obs_model = IntegralObsModel(
    obs,
    shells_df,
    Lfov=4.0*u.deg,
    t_exp_mjd=(10*u.min).to_value(u.day),
    t_init_mjd=Time("2026-01-01T11:16:00"),
    ndec=900,
    nra=1800,
    Npoints=1000,
    dt=5.0*u.s,
)

satellite_catalogue = integral_obs_model.sample_satellites(seed=1)
print(satellite_catalogue)
```

For more, see the demo notebooks:

- Compute a full sky map of number of satellite per hour of observation: [notebooks/satellite_density_maps.ipynb](https://github.com/nicolas-cerardi/Analytical-Satsky/blob/main/notebooks/satellite_density_maps.ipynb)
- Compute satellite occupancy, as a function of target declination, with error bars: [notebooks/compute_satellite_occupancy.ipynb](https://github.com/nicolas-cerardi/Analytical-Satsky/blob/main/notebooks/compute_satellite_occupancy.ipynb)
- Sample mock satellite catalogues: [notebooks/satellite_sampling.ipynb](https://github.com/nicolas-cerardi/Analytical-Satsky/blob/main/notebooks/satellite_sampling.ipynb)
- To explore various functionalities of the package, look at the walkthrough: [notebooks/walkthrough.ipynb](https://github.com/nicolas-cerardi/Analytical-Satsky/blob/main/notebooks/walkthrough.ipynb)

And the [online documentation](https://nicolas-cerardi.github.io/Analytical-Satsky/) (in progress):

- [API reference](https://nicolas-cerardi.github.io/Analytical-Satsky/api/)
- [Constellations](https://nicolas-cerardi.github.io/Analytical-Satsky/constellations/)
- [Models Description](https://nicolas-cerardi.github.io/Analytical-Satsky/model/)

## Repository structure

- `analytical_satsky/` — Python package containing the analytical model and utilities  
- `notebooks/` — Tutorials and figure reproduction examples  
- `tests/` — Unit and smoke tests executed automatically
- `docs/` — Documentation  

## Roadmap

- In development: cleanup of unused legacy functions.

## Citation

If you reuse this code in an academic work, please cite the accompanying paper:

- Cerardi, N. Tolley, E. and di Vruno, F., Forecasting the occupancy of satellite megaconstellations in SKA observations, Letter to the Editor, accepted in Astronomy and Astrophysics, 2026

```bibtex
@article{Cerardi_2026,
   title={Forecasting the occupancy of satellite megaconstellations in SKA observations},
   ISSN={1432-0746},
   url={http://dx.doi.org/10.1051/0004-6361/202558309},
   DOI={10.1051/0004-6361/202558309},
   journal={Astronomy &amp; Astrophysics},
   publisher={EDP Sciences},
   author={Cerardi, N. and Tolley, E. and di Vruno, F.},
   year={2026},
   month=Apr }
```
