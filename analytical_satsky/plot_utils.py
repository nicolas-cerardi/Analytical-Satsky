from typing import Optional

import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.units import Quantity
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes


def plot_sky_map(
    sky_map: Quantity,
    obsloc: EarthLocation,
    target_lha: Quantity,
    target_dec: Quantity,
    cmap: str = "magma",
    vmin: float = 10.0,
    vmax: float = 10000.0,
    return_fig: bool = False,
) -> Optional[tuple[Figure, Axes]]:
    """
    Plot a sky map in local horizon coordinates.

    The input map is defined on a local hour angle / declination grid and is
    projected to altitude / azimuth coordinates for the given observatory.
    Points below the horizon are masked.

    Parameters
    ----------
    sky_map : astropy.units.Quantity
        Two-dimensional sky map to plot. Must have the same shape as the grid
        defined by ``target_lha`` and ``target_dec``.
    obsloc : astropy.coordinates.EarthLocation
        Location of the observer.
    target_lha : astropy.units.Quantity
        One-dimensional local hour angle grid. Must be convertible to radians.
    target_dec : astropy.units.Quantity
        One-dimensional declination grid. Must be convertible to radians.
    cmap : str, optional
        Matplotlib colormap used for the plot. Default is ``"magma"``.
    vmin : float, optional
        Minimum value of the logarithmic color scale. Default is 10.0.
    vmax : float, optional
        Maximum value of the logarithmic color scale. Default is 10000.0.
    return_fig : bool, optional
        If True, return the Matplotlib figure and axes. Default is False.

    Returns
    -------
    tuple[matplotlib.figure.Figure, matplotlib.axes.Axes] or None
        Figure and polar axes if ``return_fig`` is True, otherwise None.

    Notes
    -----
    The color scale is logarithmic. Therefore, plotted values should be
    strictly positive after masking.
    """
   
    
    # Get Local Sidereal Time (RA at meridian)
    obs_time = Time('2026-01-01T00:00:00')  # UTC time #2025-03-18T08:00:00 #2026-01-01T00:00:00
    altaz_frame = AltAz(obstime=obs_time, location=obsloc)
    LST = obs_time.sidereal_time('apparent', longitude=obsloc.lon)
    lha_grid, dec_grid = np.meshgrid(target_lha+LST.to(u.rad), target_dec, indexing='ij')
    
    # Convert LHA/Dec to Alt/Az
    sky_coords = SkyCoord(ra=lha_grid, dec=dec_grid, frame='icrs')
    altaz_coords = sky_coords.transform_to(altaz_frame)

    alt = altaz_coords.alt.deg
    az = altaz_coords.az.deg

    # Mask points below the horizon
    mask = alt < 0
    plot_data = np.copy(sky_map)
    plot_data[mask] = np.nan

    # Step 3: Polar plot (zenith-centered)
    r = (90 - alt)
    theta = np.unwrap(np.deg2rad(az), axis=0)
    theta = np.unwrap(theta, axis=1)

    fig = plt.figure(figsize=(4.5, 4))
    ax = fig.add_subplot(111, polar=True)
    c = ax.pcolormesh(theta, r, plot_data.value, shading='auto', cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax))

    ax.set_ylim(0, 90)
    ax.set_theta_zero_location('N')  # Azimuth 0 at top (North)
    ax.set_theta_direction(1)       # Clockwise: N → E → S → W
    ax.set_yticks([30, 60, 70, 80])
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])
    ax.set_xticklabels(['N', 'E', 'S', 'W'])

    cbar_ax = fig.add_axes([0.95, 0.1, 0.04, 0.8])  # [left, bottom, width, height]
    cbar = fig.colorbar(c, cax=cbar_ax, label='nsats / h')
    
    plt.show()
    if return_fig:
        return fig, ax
    return
    