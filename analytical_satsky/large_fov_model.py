import numpy as np

import astropy.units as u
import astropy.constants as const
from astropy.coordinates import CIRS, SkyCoord

from analytical_satsky.model import compute_cosalpha, compute_geocentric_vs, satellite_density, compute_apparent_v

w_earth = 1/86164.098903691 * 2 * np.pi * u.rad / u.s

def local_basis_at_radec(ra_deg, dec_deg):
    """
    Returns orthonormal basis at (ra, dec):
      - e_east  : increasing RA direction
      - e_north : increasing Dec direction
      - e_los   : line of sight
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
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    return np.array([
        np.cos(dec) * np.cos(ra),
        np.cos(dec) * np.sin(ra),
        np.sin(dec)
    ])

def pointing_to_radec_circle(pointing_ra, pointing_dec, Lfov, Npoints=100, convention='westeast'):
    """
    Returns circle_vecs, circle_ra, circle_dec
     - circle_vecs: shape (Npoints, 3) unit vectors of the circle points in ICRS coordinates
     - circle_ra: shape (Npoints,) RA of the circle points in degrees
     - circle_dec: shape (Npoints,) Dec of the circle points
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
    Converts ICRS coordinates to CIRS, including ERA correction, to be passed to the analytical model
    All inputs and outputs are astropy quantities in rad
    """
    coords_icrs = SkyCoord(ra=ra_icrs, dec=dec_icrs, frame='icrs')
    coords_cirs = coords_icrs.transform_to(CIRS(obstime=obstime))

    ra_cirs = coords_cirs.ra.to(u.rad)
    dec_cirs = coords_cirs.dec.to(u.rad)

    ERA = obstime.earth_rotation_angle(longitude=location.lon).to(u.rad)

    lambda_cirs = (ra_cirs - ERA).wrap_at(180 * u.deg).to(u.rad)

    return lambda_cirs, dec_cirs

def unitvec_to_lonlat(vecs):
    x = vecs[:, 0]
    y = vecs[:, 1]
    z = vecs[:, 2]
    lon = np.arctan2(y, x) * u.rad
    lat = np.arcsin(np.clip(z, -1.0, 1.0)) * u.rad
    return lon, lat

def wsat_to_theta01r01(wsat, circle_vecs_cirs, obs_time, location, altaz_frame, dt_plot):
    """
    Converts wsat at the circle points to the corresponding theta0, r0 (original position) and theta1, r1 (after small displacement) in the polar plot.
    This is done by applying the small displacement given by wsat to the circle_vecs_cirs, then converting both the original and displaced vectors to altaz coordinates, and then to polar coordinates.
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
    '''equation is 
    0 = d^2 + d 2 [z_CO sin(dec) + y_CO cos(dec) sin(ra) + x_CO cos(dec) cos(ra)] + ||CO||^2 - (R_earth+hsat)^2
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

def obs_and_shells_to_nsatflux(location, pointing_ra, pointing_dec, circle_vecs, circle_ra, circle_dec, i, Nsat, hsat, Lfov, obstime, Npoints=100, dt=1):
    d, sat_lat, sat_lon = compute_d_lat_lon_gcrs(location, circle_dec, circle_ra, hsat, obstime)
    cosalpha = compute_cosalpha(d, hsat)
    V_N, V_S = compute_geocentric_vs(i, sat_lat, sat_lon, hsat)
    rho_sat = satellite_density(Nsat, sat_lat, i, hsat, d.to(u.km), cosalpha).squeeze()
    
    e_east0, e_north0, e_center = local_basis_at_radec(pointing_ra, pointing_dec)
    den = circle_vecs @ e_center
    x = (circle_vecs @ e_east0) / den
    y = (circle_vecs @ e_north0) / den
    dl_norm = np.deg2rad(Lfov / 2.0) * u.rad * 2.0 * np.pi / Npoints # boundary angular element length [rad]
    r_dot_E = circle_vecs @ e_east0
    r_dot_N = circle_vecs @ e_north0
    r_dot_C = circle_vecs @ e_center
    tx = np.roll(x, -1) - np.roll(x, 1)
    ty = np.roll(y, -1) - np.roll(y, 1)
    tnorm = np.hypot(tx, ty)
    tx = dl_norm * tx / tnorm
    ty = dl_norm * ty / tnorm
    # Keep only entering contributions
    _, flux_N, _, _ = compute_flux_vel(V_N, location, circle_dec.to(u.deg).value, circle_ra.to(u.deg).value, d, 
                              e_east0, e_north0, e_center, r_dot_E, r_dot_N, r_dot_C, tx, ty, obstime, dt)
    _, flux_S, _, _ = compute_flux_vel(V_S, location, circle_dec.to(u.deg).value, circle_ra.to(u.deg).value, d,
                              e_east0, e_north0, e_center, r_dot_E, r_dot_N, r_dot_C, tx, ty, obstime, dt)
    return  rho_sat/2 *flux_N, rho_sat/2 *flux_S

def compute_orbital_parameters_gcrs(location, obstime, inclination, dec, ra, h_sat):
    _, satlat_gcrs, satlon_gcrs = compute_d_lat_lon_gcrs(location, dec, ra, h_sat, obstime)
    lonlambda = np.arcsin(np.tan(satlat_gcrs)/np.tan(inclination.to(u.rad))) #*u.rad
    
    Omega_N = (satlon_gcrs.to(u.rad) - lonlambda).to(u.rad) # in rad
    Omega_S = (satlon_gcrs.to(u.rad) + lonlambda).to(u.rad) + np.pi*u.rad # in rad
    periapsis_N = np.arccos(np.cos((lonlambda).to(u.rad))*np.cos(satlat_gcrs)).to(u.deg)*np.sign(satlat_gcrs)
    periapsis_S = np.pi*u.rad - np.arccos(np.cos((lonlambda).to(u.rad))*np.cos(satlat_gcrs)).to(u.deg)*np.sign(satlat_gcrs) #+ np.pi*u.rad
    return Omega_N, Omega_S, periapsis_N, periapsis_S

def angular_distance(ra1, dec1, ra2, dec2):
    """
    ra1 dec1 are just scalars, ra2 and dec2 are arrays of the same shape, and the output is an array of the same shape as ra2 and dec2
    """
    x = np.cos(dec1) * np.sin(dec2) - np.sin(dec1) * np.cos(dec2) * np.cos(ra2 - ra1)
    y = np.cos(dec2) * np.sin(ra2 - ra1)
    z = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra2 - ra1)
    return np.arctan2(np.sqrt(x**2 + y**2), z)
