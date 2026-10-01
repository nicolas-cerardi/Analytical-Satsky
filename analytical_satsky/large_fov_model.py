import numpy as np

import astropy.units as u
import astropy.constants as const
from astropy.coordinates import CIRS, SkyCoord

from analytical_satsky.model import compute_apparent_v

w_earth = 1/86164.098903691 * 2 * np.pi * u.rad / u.s

def local_basis_at_radec(ra_deg, dec_deg):
    """
    Compute the local orthonormal tangent basis at a given sky position.

    Parameters
    ----------
    ra_deg : float
        Right ascension of the reference point, in degrees.
    dec_deg : float
        Declination of the reference point, in degrees.

    Returns
    -------
    e_east : numpy.ndarray
        Unit vector pointing in the direction of increasing RA, shape ``(3,)``.
    e_north : numpy.ndarray
        Unit vector pointing in the direction of increasing DEC, shape ``(3,)``.
    e_los : numpy.ndarray
        Unit vector along the LOS, towards the reference point, shape ``(3,)``.
    """
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)

    e_los = np.array([
        np.cos(dec) * np.cos(ra),
        np.cos(dec) * np.sin(ra),
        np.sin(dec)
    ])

    # tangent basis
    e_east = np.array([
        -np.sin(ra),
         np.cos(ra),
         0.0
    ])

    e_north = np.array([
        -np.sin(dec) * np.cos(ra),
        -np.sin(dec) * np.sin(ra),
         np.cos(dec)
    ])

    # already normalized, but keep it clean
    e_east /= np.linalg.norm(e_east)
    e_north /= np.linalg.norm(e_north)
    e_los /= np.linalg.norm(e_los)

    return e_east, e_north, e_los

def radec_to_unitvec(ra_deg, dec_deg):
    """
    Convert a sky position to a unit vector in ICRS coordinates.

    Parameters
    ----------
    ra_deg : float
        Right ascension, in degrees.
    dec_deg : float
        Declination, in degrees.

    Returns
    -------
    numpy.ndarray
        Unit vector towards ``(ra_deg, dec_deg)``, shape ``(3,)``.
    """
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    return np.array([
        np.cos(dec) * np.cos(ra),
        np.cos(dec) * np.sin(ra),
        np.sin(dec)
    ])

def pointing_to_radec_circle(pointing_ra, pointing_dec, Lfov, Npoints=100):
    """
    Sample points on a circle of a given radius around a pointing direction.

    Parameters
    ----------
    pointing_ra : float
        Right ascension of the circle center, in degrees, ICRS.
    pointing_dec : float
        Declination of the circle center, in degrees, ICRS.
    Lfov : float
        Circle diameter, in degrees.
    Npoints : int, optional
        Number of points sampled on the circle. Default is 100.

    Returns
    -------
    circle_vecs : numpy.ndarray
        Unit vectors of the circle points, in ICRS coordinates, shape
        ``(Npoints, 3)``.
    circle_ra : numpy.ndarray
        Right ascension of the circle points, in radians, shape
        ``(Npoints,)``.
    circle_dec : numpy.ndarray
        Declination of the circle points, in radians, shape ``(Npoints,)``.
    """
    center = radec_to_unitvec(pointing_ra, pointing_dec)
    if np.abs(center[2]) < 0.9:
        ref = np.array([0.0, 0.0, 1.0])
    else:
        ref = np.array([1.0, 0.0, 0.0])

    v1 = np.cross(ref, center)
    v1 = v1 / np.linalg.norm(v1)

    v2 = np.cross(center, v1)
    v2 = v2 / np.linalg.norm(v2)

    # Circle points on the sphere
    phis = np.linspace(0, 2*np.pi, Npoints, endpoint=False)
    rho = np.deg2rad(Lfov / 2.0)
    cosphi = np.cos(phis)[:, None]
    sinphi = np.sin(phis)[:, None]

    b = (
        np.cos(rho) * center[None, :] +
        np.sin(rho) * (cosphi * v1[None, :] + sinphi * v2[None, :])
    )

    b = b / np.linalg.norm(b, axis=1)[:, None]

    x, y, z = b[:,0], b[:,1], b[:,2]
    circle_ra = np.arctan2(y, x)
    circle_dec = np.arcsin(z)

    return np.array(b), np.array(circle_ra), np.array(circle_dec)

def ICRS_to_CIRS(ra_icrs, dec_icrs, obstime, location):
    """
    Convert ICRS coordinates to CIRS, ERA-corrected for use in the analytical model.

    Parameters
    ----------
    ra_icrs : astropy.units.Quantity
        Right ascension, ICRS. Must be angular.
    dec_icrs : astropy.units.Quantity
        Declination, ICRS. Must be angular.
    obstime : astropy.time.Time
        Observation time.
    location : astropy.coordinates.EarthLocation
        Observer location.

    Returns
    -------
    lambda_cirs : astropy.units.Quantity
        CIRS right ascension minus the Earth rotation angle, in radians.
    dec_cirs : astropy.units.Quantity
        Declination in the CIRS frame, in radians.
    """
    coords_icrs = SkyCoord(ra=ra_icrs, dec=dec_icrs, frame='icrs')
    coords_cirs = coords_icrs.transform_to(CIRS(obstime=obstime))

    ra_cirs = coords_cirs.ra.to(u.rad)
    dec_cirs = coords_cirs.dec.to(u.rad)

    ERA = obstime.earth_rotation_angle(longitude=location.lon).to(u.rad)

    lambda_cirs = (ra_cirs - ERA).wrap_at(180 * u.deg).to(u.rad)

    return lambda_cirs, dec_cirs

def unitvec_to_lonlat(vecs):
    """
    Convert unit vectors to longitude/latitude angles.

    Parameters
    ----------
    vecs : numpy.ndarray
        Unit vectors, shape ``(N, 3)``.

    Returns
    -------
    lon : astropy.units.Quantity
        Longitude, in radians, shape ``(N,)``.
    lat : astropy.units.Quantity
        Latitude, in radians, shape ``(N,)``.
    """
    x = vecs[:, 0]
    y = vecs[:, 1]
    z = vecs[:, 2]
    lon = np.arctan2(y, x) * u.rad
    lat = np.arcsin(np.clip(z, -1.0, 1.0)) * u.rad
    return lon, lat

def wsat_to_theta01r01(wsat, circle_vecs_cirs, obs_time, location, altaz_frame, dt_plot):
    """
    Convert apparent satellite velocities to polar-plot coordinates.

    Displaces the circle points by ``wsat * dt_plot``, then converts both
    the original and displaced positions to Alt/Az and then to polar
    (theta, r) coordinates for plotting.

    Parameters
    ----------
    wsat : astropy.units.Quantity
        Apparent angular velocity vectors at the circle points, shape
        ``(N, 3)``.
    circle_vecs_cirs : array-like
        Unit vectors of the circle points, in the CIRS frame, shape
        ``(N, 3)``.
    obs_time : astropy.time.Time
        Observation time.
    location : astropy.coordinates.EarthLocation
        Observer location.
    altaz_frame : astropy.coordinates.AltAz
        Alt/Az frame to transform into.
    dt_plot : astropy.units.Quantity
        Small time step used to displace the circle points.

    Returns
    -------
    theta0 : numpy.ndarray
        Azimuth of the original position, in radians, shape ``(N,)``.
    r0 : numpy.ndarray
        Zenith distance of the original position, in degrees, shape
        ``(N,)``.
    theta1 : numpy.ndarray
        Azimuth of the displaced position, in radians, shape ``(N,)``.
    r1 : numpy.ndarray
        Zenith distance of the displaced position, in degrees, shape
        ``(N,)``.
    """
    # Original directions on the sphere
    b0 = np.asarray(circle_vecs_cirs, dtype=float)   # shape (N,3)

    # Small displacement using the apparent angular velocity
    db = (wsat * dt_plot).to_value(u.rad)          # shape (N,3)
    b1 = b0 + db
    b1 /= np.linalg.norm(b1, axis=1, keepdims=True)

    # Convert both sets of vectors back to spherical coordinates
    # IMPORTANT:
    # this lon/lat must match the convention used to build circle_vecs_cirs.
    lon0, lat0 = unitvec_to_lonlat(b0)
    lon1, lat1 = unitvec_to_lonlat(b1)

    # If your geometry vectors are built with lambda = -H (recommended),
    # then true CIRS RA is:
    ERA = obs_time.earth_rotation_angle(longitude=location.lon).to(u.rad)
    ra0 = (ERA + lon0).wrap_at(180 * u.deg)
    ra1 = (ERA + lon1).wrap_at(180 * u.deg)

    sky0_cirs = SkyCoord(ra=ra0, dec=lat0, frame=CIRS(obstime=obs_time))
    sky1_cirs = SkyCoord(ra=ra1, dec=lat1, frame=CIRS(obstime=obs_time))

    altaz0 = sky0_cirs.transform_to(altaz_frame)
    altaz1 = sky1_cirs.transform_to(altaz_frame)

    theta0 = altaz0.az.to_value(u.rad)
    theta1 = altaz1.az.to_value(u.rad)

    # polar radius = zenith distance, in radians
    r0 = (90 - altaz0.alt.deg)
    r1 = (90 - altaz1.alt.deg)
    return theta0, r0, theta1, r1

def compute_d_lat_lon_gcrs(location, dec, ra, hsat, obstime):
    '''
    Compute the distance and satellite GCRS position along a line of sight.

    Solves ``0 = d^2 + d 2 [z_CO sin(dec) + y_CO cos(dec) sin(ra) +
    x_CO cos(dec) cos(ra)] + ||CO||^2 - (R_earth+hsat)^2`` for ``d``, the
    distance from the observer to the shell along the line of sight.

    Parameters
    ----------
    location : astropy.coordinates.EarthLocation
        Observer location.
    dec : astropy.units.Quantity
        Declination of the line of sight, GCRS frame. Must be angular.
    ra : astropy.units.Quantity
        Right ascension of the line of sight, GCRS frame. Must be angular.
    hsat : astropy.units.Quantity
        Altitude of the shell.
    obstime : astropy.time.Time
        Observation time.

    Returns
    -------
    d : astropy.units.Quantity
        Distance from the observer to the shell along the line of sight.
    sat_lat : astropy.units.Quantity
        Satellite latitude, GCRS frame.
    sat_lon : astropy.units.Quantity
        Satellite longitude, GCRS frame.
    '''
    #1. compute obs location xyz in gcrs
    obs_xyz_gcrs = location.get_gcrs_posvel(obstime)[0].xyz.to(u.km)

    #2. compute the aeq, beq, ceq
    aeq = 1.0
    beq = 2 * (obs_xyz_gcrs[2]*np.sin(dec) + obs_xyz_gcrs[1]*np.cos(dec)*np.sin(ra) + obs_xyz_gcrs[0]*np.cos(dec)*np.cos(ra))
    ceq = np.sum(obs_xyz_gcrs**2) - (const.R_earth + hsat)**2
    discriminant = beq**2 - 4*aeq*ceq
    #there is one positive root, that we retain
    d = (-beq + np.sqrt(discriminant)) / 2

    #3. Then lat and lon of satellite in GCRS frame
    sat_lat = np.arcsin((obs_xyz_gcrs[2]+d*np.sin(dec)) / (const.R_earth + hsat))
    sat_lon = np.arctan2((obs_xyz_gcrs[1]+d*np.cos(dec)*np.sin(ra)), (obs_xyz_gcrs[0]+d*np.cos(dec)*np.cos(ra)))
    #print(sat_lat.unit, sat_lon.unit)
    return d, sat_lat, sat_lon

def compute_flux_vel(V_N, location, circle_dec, circle_ra, d, e_east0, e_north0, e_center, r_dot_E, r_dot_N, r_dot_C, tx, ty, obstime, dt):
    """
    Compute apparent velocity and satellite entry flux at FoV boundary points.

    Parameters
    ----------
    V_N : astropy.units.Quantity
        Geocentric velocity vectors of the satellites, shape ``(3, N)``.
    location : astropy.coordinates.EarthLocation
        Observer location.
    circle_dec : array-like
        Declination of the boundary points, in degrees, shape ``(N,)``.
    circle_ra : array-like
        Right ascension of the boundary points, in degrees, shape ``(N,)``.
    d : astropy.units.Quantity
        Distance from the observer to the shell at each boundary point,
        shape ``(N,)``.
    e_east0, e_north0, e_center : numpy.ndarray
        Local tangent basis at the pointing center, each shape ``(3,)``.
    r_dot_E, r_dot_N, r_dot_C : numpy.ndarray
        Projections of the boundary unit vectors onto ``e_east0``,
        ``e_north0``, ``e_center``, each shape ``(N,)``.
    tx, ty : numpy.ndarray
        Tangent vector components along the FoV boundary, each shape
        ``(N,)``.
    obstime : astropy.time.Time
        Observation time.
    dt : astropy.units.Quantity
        Time step used to convert apparent angular velocity into a
        displacement.

    Returns
    -------
    w_sat : astropy.units.Quantity
        Apparent angular velocity vectors at the boundary points, shape
        ``(N, 3)``.
    flux_vel : numpy.ndarray
        Rate of satellites entering the field of view across the boundary,
        shape ``(N,)``. Negative (exiting) contributions are clipped to
        zero.
    vx, vy : numpy.ndarray
        Apparent velocity components in the tangent-plane projection, each
        shape ``(N,)``.
    """
    v_obs_gcrs = location.get_gcrs_posvel(obstime)[1].xyz #w_earth*const.R_earth*np.cos(obslat.to(u.rad))*np.array([-np.sin(obslon.to(u.rad)), np.cos(obslon.to(u.rad)), 0])
    
    topo_V_N = V_N - v_obs_gcrs[:,np.newaxis]
    app_V_N = compute_apparent_v(topo_V_N, np.deg2rad(circle_dec)*u.rad, np.deg2rad(circle_ra)*u.rad)
    w_sat_N = app_V_N / d.to(u.km) *u.rad
    #print("w_sat_N in", w_sat_N.unit)
    dr = w_sat_N.T  # shape (Npoints, 3)
    dr_dot_E = dr @ e_east0
    dr_dot_N = dr @ e_north0
    dr_dot_C = dr @ e_center
    vx = (dr_dot_E * r_dot_C - r_dot_E * dr_dot_C) / (r_dot_C**2)
    vy = (dr_dot_N * r_dot_C - r_dot_N * dr_dot_C) / (r_dot_C**2)
    flux_vel = tx * (dt * vy) - ty * (dt * vx)
    flux_vel = np.where(flux_vel > 0.0, flux_vel, 0.0)
    
    return w_sat_N.T, flux_vel, vx, vy

def angular_distance(ra1, dec1, ra2, dec2):
    """
    Compute the angular distance between a reference point and a grid.

    Parameters
    ----------
    ra1 : float
        Right ascension of the reference point, in radians.
    dec1 : float
        Declination of the reference point, in radians.
    ra2 : numpy.ndarray
        Right ascension grid, in radians.
    dec2 : numpy.ndarray
        Declination grid, in radians.

    Returns
    -------
    numpy.ndarray
        Angular distance between ``(ra1, dec1)`` and each point of the
        grid, in radians, same shape as ``ra2``/``dec2``.
    """
    x = np.cos(dec1) * np.sin(dec2) - np.sin(dec1) * np.cos(dec2) * np.cos(ra2 - ra1)
    y = np.cos(dec2) * np.sin(ra2 - ra1)
    z = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra2 - ra1)
    return np.arctan2(np.sqrt(x**2 + y**2), z)

def satellite_positions_gcrs(satellite_catalogue, t_array_s):
    """
    Propagate a sampled satellite catalogue to GCRS positions over time.

    Builds each satellite's circular-orbit trajectory from its orbital
    elements, using the periapsis then inclination then longitude-of-
    ascending-node rotation sequence (true anomaly measured from periapsis
    as ``omega_orbit * t_array_s``, matching the phase convention used by
    ``SingleShellFoV``/``SingleShellFlux``'s ``compute_orbital_params_samples``).

    Parameters
    ----------
    satellite_catalogue : pandas.DataFrame
        Satellite catalogue as returned by ``IntegralObsModel.sample_satellites``
        (or any of the other ``sample_satellites`` methods in ``shell_obs.py``).
        Must have columns ``i``, ``h``, ``long_asc_node``, ``periapsis``, all
        plain floats (degrees for angles, km for altitude).
    t_array_s : numpy.ndarray
        Times at which to evaluate each satellite's position, in seconds
        since the catalogue's reference time (``t_init_mjd``), shape
        ``(n_t,)``.

    Returns
    -------
    numpy.ndarray
        Satellite positions in the GCRS frame, in km, shape
        ``(nsat, 3, n_t)``.
    """
    nsat = len(satellite_catalogue)
    n_t = len(t_array_s)

    R_earth_km = const.R_earth.to_value(u.km)
    GM_earth_km3_s2 = const.GM_earth.to_value(u.km**3 / u.s**2)

    h_km = satellite_catalogue['h'].values.astype(float)
    R_orbit_km = R_earth_km + h_km
    omega_orbit = np.sqrt(GM_earth_km3_s2 / R_orbit_km**3)

    i_rad = np.deg2rad(satellite_catalogue['i'].values.astype(float))
    Omega_rad = np.deg2rad(satellite_catalogue['long_asc_node'].values.astype(float))
    periapsis_rad = np.deg2rad(satellite_catalogue['periapsis'].values.astype(float))

    nu = omega_orbit[:, None] * t_array_s[None, :]  # (nsat, n_t), rad

    r_perifocal = np.stack([
        R_orbit_km[:, None] * np.cos(nu),
        R_orbit_km[:, None] * np.sin(nu),
        np.zeros((nsat, n_t)),
    ], axis=1)  # (nsat, 3, n_t), km

    def _rot_z(theta):
        c, s, z, o = np.cos(theta), np.sin(theta), np.zeros_like(theta), np.ones_like(theta)
        return np.stack([
            np.stack([c, -s, z], axis=-1),
            np.stack([s,  c, z], axis=-1),
            np.stack([z,  z, o], axis=-1),
        ], axis=-2)  # (..., 3, 3)

    def _rot_x(theta):
        c, s, z, o = np.cos(theta), np.sin(theta), np.zeros_like(theta), np.ones_like(theta)
        return np.stack([
            np.stack([o, z,  z], axis=-1),
            np.stack([z, c, -s], axis=-1),
            np.stack([z, s,  c], axis=-1),
        ], axis=-2)  # (..., 3, 3)

    R_total = _rot_z(Omega_rad) @ _rot_x(i_rad) @ _rot_z(periapsis_rad)  # (nsat, 3, 3)
    return np.einsum('sij,sjt->sit', R_total, r_perifocal)

def positions_to_lm(positions_gcrs, obsloc, times, e_east, e_north, e_center):
    """
    Project GCRS positions into l,m direction cosines from a moving observer.

    Parameters
    ----------
    positions_gcrs : numpy.ndarray
        Positions in the GCRS frame, in km, shape ``(nsat, 3, n_t)``.
    obsloc : astropy.coordinates.EarthLocation
        Observer location.
    times : astropy.time.Time
        Observation times, shape ``(n_t,)``, matching ``positions_gcrs``'s
        last axis.
    e_east, e_north, e_center : numpy.ndarray
        Local tangent basis at the pointing center (see
        ``local_basis_at_radec``), each shape ``(3,)``.

    Returns
    -------
    l_deg, m_deg : numpy.ndarray
        East-west / north-south direction cosines, in degrees (flat-sky
        approximation), shape ``(nsat, n_t)``.
    n : numpy.ndarray
        Third direction-cosine component, shape ``(nsat, n_t)``. Negative
        values mean the position is on the far/antipodal hemisphere from
        the pointing center, where ``l``/``m`` alone can't distinguish it
        from a genuine nearby crossing.
    """
    obs_pos_gcrs = obsloc.get_gcrs_posvel(times)[0].xyz.to_value(u.km)  # (3, n_t)
    rel = positions_gcrs - obs_pos_gcrs[None, :, :]
    rel_unit = rel / np.linalg.norm(rel, axis=1, keepdims=True)

    l_deg = np.rad2deg(np.einsum('sit,i->st', rel_unit, e_east))
    m_deg = np.rad2deg(np.einsum('sit,i->st', rel_unit, e_north))
    n = np.einsum('sit,i->st', rel_unit, e_center)
    return l_deg, m_deg, n
