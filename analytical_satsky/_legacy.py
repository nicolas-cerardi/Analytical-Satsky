import numpy as np
from astropy.coordinates import EarthLocation
import astropy.units as u
from astropy.units import Quantity

from .model import (
    compute_d_phi,
    compute_cosalpha,
    satellite_density,
    compute_wsat,
    compute_nsats,
    draw_passes
)

def run_analytical_ska(obslat, obslon, obstime, tobs, target_dec, target_ra, shells_i, shells_h, shells_n, Lfov, nstat=100):
    """
    Deprecated - early version of the code, not tested and not used in the current version of the package.
    """
    location = EarthLocation(lat=obslat, lon=obslon)
    LST = obstime.sidereal_time('apparent', longitude=location.lon)
    target_dec = target_dec.to(u.rad)
    target_ra = target_ra.to(u.rad) - LST.to(u.rad)
    all_inits, all_ts = [],[]
    nshells = len(shells_i)
    for idx in range(nshells):
        i = (shells_i[idx]*u.deg).to(u.rad)
        Nsat = shells_n[idx]
        hsat = shells_h[idx] * u.km
        d, lat, lon = compute_d_phi(obslat, hsat, target_dec, target_ra)
        cosalpha = compute_cosalpha(d, hsat)
        rho_sat = satellite_density(Nsat, lat, i, hsat, d, cosalpha)
        wsat = compute_wsat(i, lat, lon, hsat, obslat, target_dec, target_ra)
        
        nsats = np.nan_to_num(compute_nsats(rho_sat, Lfov, wsat/d.to(u.m)*u.rad, tobs))
        n_sample = np.random.poisson(nsats.value.squeeze(), size=(nstat,)).squeeze()
        if n_sample.any()>0:
            tmp_inits, tmp_ts = draw_passes(n_sample, Lfov, wsat, d, tobs)
            all_inits.append(tmp_inits)
            all_ts.append(tmp_ts.value)
    all_inits = np.concatenate(all_inits, axis=1)
    all_ts = np.concatenate(all_ts, axis=1)
    return all_inits, all_ts

def compute_exposure_fraction(
    tobs: Quantity,
    ntimestep: int,
    all_inits: np.ndarray,
    all_ts: np.ndarray,
) -> np.ndarray:
    """
    Compute the fraction of observing time during which at least one satellite
    is present in the effective beam.

    Parameters
    ----------
    tobs : astropy.units.Quantity
        Observation duration. Must be convertible to seconds.
    ntimestep : int
        Number of time samples used to discretize the interval ``[0, tobs]``.
    all_inits : numpy.ndarray
        Ingress times into the effective beam for each statistical realisation.
        Expected shape is ``(nstat, nevents)``, where each row contains the
        ingress times for one realisation. Values are assumed to be in seconds.
    all_ts : numpy.ndarray
        Fly-through durations across the effective beam for each statistical
        realisation. Must have the same shape as ``all_inits``. Values are
        assumed to be in seconds.

    Returns
    -------
    numpy.ndarray
        Exposure fraction for each statistical realisation. The returned array
        has shape ``(nstat,)``.

    Notes
    -----
    Deprecated implementation based on explicit boolean masking. 
    Please use ``compute_occupancy_fraction`` instead, which uses a more efficient approach.
    """
    
    timeframe = np.linspace(0,tobs.to(u.s).value,ntimestep)
    nsat_at_t = np.zeros((all_inits.shape[0], ntimestep))
    
    for j in range(all_inits.shape[0]):
        for init, tpass in zip(all_inits[j], all_ts[j]):
            nsat_at_t[j, (timeframe>init)&(timeframe<init+tpass)] += 1
    return np.mean(nsat_at_t>=1, axis=1)