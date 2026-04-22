import numpy as np
import astropy.units as u
from analytical_satsky import compute_total_satellite_density, load_constellation

def minimal_code(ndec, nra):
    shells_df = load_constellation("starlink_filing2")

    n_satellite_in_obs = compute_total_satellite_density(
        -30*u.deg, 
        116.*u.deg, 
        shells_df, 
        np.linspace(-40., -20., num=ndec)*u.deg,
        np.linspace(100., 120., num=nra)*u.deg, 
        10.*u.deg, 
        3600*u.s
    )
    return n_satellite_in_obs

def test_shape():
    ndec, nra = 10, 10
    n_satellite_in_obs = minimal_code(ndec, nra)
    assert n_satellite_in_obs.shape == (ndec, nra)
    
def test_units():
    #should be unitless
    ndec, nra = 10, 10
    n_satellite_in_obs = minimal_code(ndec, nra)
    assert n_satellite_in_obs.unit == u.dimensionless_unscaled

def test_positive():
    ndec, nra = 10, 10
    n_satellite_in_obs = minimal_code(ndec, nra)
    assert np.all(n_satellite_in_obs.value >= 0)

