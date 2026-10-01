import numpy as np
import astropy.units as u
from astropy.units import Quantity

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