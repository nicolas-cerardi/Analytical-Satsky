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

from analytical_satsky.large_fov_model import local_basis_at_radec, satellite_positions_gcrs, positions_to_lm
from analytical_satsky.shell_obs import MultiShellFlux


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
    projected to altitude / azimuth coordinates for the given observatory,
    using a zenith-centred polar projection with cardinal directions labelled
    ``N``, ``E``, ``S``, ``W``. Points below the horizon are masked.

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

    Examples
    --------
    >>> import astropy.units as u
    >>> from astropy.coordinates import EarthLocation
    >>> from analytical_satsky import plot_sky_map
    >>> plot_sky_map(
    ...     sky_map=sky_map,
    ...     obsloc=EarthLocation.of_site("greenwich"),
    ...     target_lha=lha_grid * u.deg,
    ...     target_dec=dec_grid * u.deg,
    ... )

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


def plot_satellite_tracks_lm(
    integral_obs_model,
    satellite_catalogue,
    n_t: int = 500,
    fov_lm_deg: float = 10.0,
    ax: Optional[Axes] = None,
    return_fig: bool = False,
) -> Optional[tuple[Figure, Axes]]:
    """
    Plot sampled satellite tracks in the l,m plane around the pointing center.

    Propagates each satellite in ``satellite_catalogue`` over the exposure
    covered by ``integral_obs_model``, and projects its position into l,m
    direction cosines (flat-sky approximation) as seen from the moving
    observer, together with the field-of-view boundary at the initial
    observation time.

    Parameters
    ----------
    integral_obs_model : analytical_satsky.shell_obs.IntegralObsModel
        Model used to generate ``satellite_catalogue``. Supplies the
        observer location, pointing, field of view, and exposure duration.
    satellite_catalogue : pandas.DataFrame
        Satellite catalogue, as returned by
        ``integral_obs_model.sample_satellites()``.
    n_t : int, optional
        Number of time samples used to propagate each satellite's track.
        Default is 500.
    fov_lm_deg : float, optional
        Full width/height of the plotted field, in degrees. Default is 10.0.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If not given, a new figure and axes are created.
    return_fig : bool, optional
        If True, return the Matplotlib figure and axes. Default is False.

    Returns
    -------
    tuple[matplotlib.figure.Figure, matplotlib.axes.Axes] or None
        Figure and axes if ``return_fig`` is True, otherwise None.

    Notes
    -----
    Points on the far/antipodal hemisphere from the pointing center (where
    the direction cosines ``l``/``m`` alone can't distinguish a genuine
    nearby crossing from one on the opposite side of the sky) are masked out
    of the plotted tracks.
    """
    nsat = len(satellite_catalogue)
    t_array_s = np.linspace(0, integral_obs_model.t_exp_mjd * 86400, n_t)

    pointing_ra_deg = integral_obs_model.initial_multi_shell_fov.pointing_ra.to_value(u.deg)
    pointing_dec_deg = integral_obs_model.initial_multi_shell_fov.pointing_dec.to_value(u.deg)
    e_east0, e_north0, e_center0 = local_basis_at_radec(pointing_ra_deg, pointing_dec_deg)

    # FoV boundary at the initial observation time, built the same way
    # IntegralObsModel.sample_satellites builds it internally.
    multi_shell_flux_t0 = MultiShellFlux(
        obs=integral_obs_model.obs,
        shells_df=integral_obs_model.shells_df,
        Npoints=integral_obs_model.Npoints,
        Lfov=integral_obs_model.Lfov,
        t_mjd=integral_obs_model.t_init_mjd,
        t_init_mjd=integral_obs_model.t_init_mjd,
        dt=integral_obs_model.dt,
    )
    circle_l_deg = np.rad2deg(multi_shell_flux_t0.circle_vecs @ e_east0)
    circle_m_deg = np.rad2deg(multi_shell_flux_t0.circle_vecs @ e_north0)
    circle_l_deg = np.append(circle_l_deg, circle_l_deg[0])
    circle_m_deg = np.append(circle_m_deg, circle_m_deg[0])

    positions_gcrs = satellite_positions_gcrs(satellite_catalogue, t_array_s)
    obs_times = integral_obs_model.t_init_mjd + t_array_s * u.s
    l_traj_deg, m_traj_deg, n_traj = positions_to_lm(
        positions_gcrs, integral_obs_model.obsloc, obs_times, e_east0, e_north0, e_center0
    )
    l_traj_deg = np.where(n_traj > 0, l_traj_deg, np.nan)
    m_traj_deg = np.where(n_traj > 0, m_traj_deg, np.nan)

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    else:
        fig = ax.figure

    ax.plot(circle_l_deg, circle_m_deg, linestyle='--', color='gray', label='FoV boundary')
    for isat in range(nsat):
        ax.plot(l_traj_deg[isat], m_traj_deg[isat], lw=1)
    ax.plot(0, 0, marker='+', color='red', markersize=15, markeredgewidth=2, label='pointing center')
    ax.set_xlim(fov_lm_deg / 2, -fov_lm_deg / 2)  # east to the left, as on-sky convention
    ax.set_ylim(-fov_lm_deg / 2, fov_lm_deg / 2)
    ax.set_xlabel('l [deg] (east-west)')
    ax.set_ylabel('m [deg] (north-south)')
    ax.set_title('Satellite tracks in the l,m plane')
    ax.set_aspect('equal')
    ax.legend()

    plt.show()
    if return_fig:
        return fig, ax
    return
    