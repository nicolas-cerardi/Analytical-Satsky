import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

def plot_sky_map(sky_map, obslat, obslon, target_ra, target_dec, cmap='magma', vmin=10., vmax=10000., return_fig=False):
    '''
    Plots a sky map in ra-dec coordinates.

    Parameters
    ----------
    sky_map : 2D array
        Sky map data to plot (e.g., satellite density).
    obslat : float
        Observer latitude in degrees.
    obslon : float
        Observer longitude in degrees.
    target_ra : array
        Target right ascension grid in degrees.
    target_dec : array
        Target declination grid in degrees.
    cmap : str, optional
        Colormap to use for the plot. Default is 'magma'.
    vmin : float, optional
        Minimum value for color scale. Default is 10.0.
    vmax : float, optional
        Maximum value for color scale. Default is 10000.0.
    return_fig : bool, optional
        If True, returns the figure and axis objects. Default is False.
    
    Returns
    -------
    fig, ax : matplotlib Figure and Axes objects (if return_fig is True)
        The figure and axis of the plot.
    '''
   
    
    # Get Local Sidereal Time (RA at meridian)
    obs_time = Time('2025-03-18T08:00:00')  # UTC time
    location = EarthLocation(lat=obslat, lon=obslon)
    altaz_frame = AltAz(obstime=obs_time, location=location)
    LST = obs_time.sidereal_time('apparent', longitude=location.lon)
    ra_grid, dec_grid = np.meshgrid(target_ra+LST.to(u.rad), target_dec, indexing='ij') #+LST.to(u.rad)
    
    # Convert RA/Dec to Alt/Az
    sky_coords = SkyCoord(ra=ra_grid, dec=dec_grid, frame='icrs') #*u.deg
    altaz_coords = sky_coords.transform_to(altaz_frame)

    alt = altaz_coords.alt.deg
    az = altaz_coords.az.deg

    # Mask points below the horizon
    mask = alt < 0 #dec_grid<0 #
    plot_data = np.copy(sky_map)
    plot_data[mask] = np.nan

    # Step 3: Polar plot (zenith-centered)
    r = (90 - alt)
    theta = np.unwrap(np.deg2rad(az), axis=0)
    theta = np.unwrap(theta, axis=1)

    fig = plt.figure(figsize=(4.5, 4))
    ax = fig.add_subplot(111, polar=True)
    c = ax.pcolormesh(theta, r, plot_data.value, shading='auto', cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax)) #vmin=0.04, vmax=10000

    ax.set_ylim(0, 90)
    ax.set_theta_zero_location('N')  # Azimuth 0 at top (North)
    ax.set_theta_direction(-1)       # Clockwise: N → E → S → W
    ax.set_yticks([30, 60, 70, 80])
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])

    cbar_ax = fig.add_axes([0.95, 0.1, 0.015, 0.8])  # [left, bottom, width, height]
    cbar = fig.colorbar(c, cax=cbar_ax, label='nsats / h')

    plt.tight_layout()
    plt.show()
    if return_fig:
        return fig, ax
    return
    