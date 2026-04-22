"""
This testing file provides test for:
 - units in outputs of internal routines
 - check values in some simple cases (e.g., zero satellite density for high declinations or low inclinations)
"""

import numpy as np
import astropy.units as u
from analytical_satsky import compute_total_satellite_density   
from analytical_satsky.model import compute_d_phi, compute_cosalpha, compute_wsat, compute_nsats
import pandas as pd

def test_zero_density_low_i():
    #create fake Dataframe with 1 shell with low inclination
    shells_df = pd.DataFrame({
        'n': [1000],
        'h': [500],
        'i': [10]
    })

    n_satellite_in_obs = compute_total_satellite_density(
        -30*u.deg, 
        116.*u.deg, 
        shells_df, 
        np.array([-30.])*u.deg,
        np.linspace(100., 120., num=10)*u.deg, 
        10.*u.deg, 
        3600*u.s
    )
    assert np.all(n_satellite_in_obs.value == 0)

def test_zero_density_high_dec():
    #create fake Dataframe with 1 shell with high inclination
    shells_df = pd.DataFrame({
        'n': [1000],
        'h': [500],
        'i': [30]
    })

    n_satellite_in_obs = compute_total_satellite_density(
        -30*u.deg, 
        116.*u.deg, 
        shells_df, 
        np.array([-50.])*u.deg,
        np.linspace(100., 120., num=10)*u.deg, 
        10.*u.deg, 
        3600*u.s
    )
    assert np.all(n_satellite_in_obs.value == 0)

def test_unit_d_phi():
    obslat_rad = (30*u.deg).to(u.rad)
    hsat_m = 5e5*u.m
    target_dec_rad = (np.array([30])*u.deg).to(u.rad)
    target_lha_rad = (np.array([0])*u.deg).to(u.rad)

    d, lat, lon = compute_d_phi(obslat_rad, hsat_m, target_dec_rad, target_lha_rad)
    assert d.unit == u.m
    assert lat.unit == u.rad
    assert lon.unit == u.rad

def test_unit_cosalpha():
    d = 6e5*u.m
    hsat_m = 5e5*u.m

    cosalpha = compute_cosalpha(d, hsat_m)
    assert cosalpha.unit == u.dimensionless_unscaled

def test_unit_wsat():
    i_rad = (35*u.deg).to(u.rad)
    lat = 0*u.rad
    lon = 0*u.rad
    hsat_m = 5e5*u.m
    obslat_rad = (30*u.deg).to(u.rad)
    obslon_rad = 0*u.rad
    target_dec_rad = (np.array([30])*u.deg).to(u.rad)
    target_lha_rad = (np.array([0])*u.deg).to(u.rad)

    wsat = compute_wsat(i_rad, lat, lon, hsat_m, obslat_rad, obslon_rad, target_dec_rad, target_lha_rad)
    assert wsat.unit == u.m / u.s

def test_unit_nsats():
    rho_sat = 1/u.rad**2
    Lfov_rad = (10*u.deg).to(u.rad)
    wsat = 1*u.m/u.s
    tobs = 3600*u.s
    d= 6e5*u.m

    nsats = compute_nsats(rho_sat, Lfov_rad, wsat/d*u.rad, tobs)
    assert nsats.unit == u.dimensionless_unscaled


