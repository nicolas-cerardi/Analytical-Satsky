import numpy as np

import os
import scipy
import skyfield
from skyfield.positionlib import ICRF
import astropy.units as u
import astropy.constants as const
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time

w_earth = 1/86164.098903691 * 2 * np.pi * u.rad / u.s

def single_sat_density(lat, i, hsat):
    '''equation A.14
    P(lat, i hsat) = 1/(2 pi^2 (Re + hsat)^2 (sin^2 i - sin^2 lat)^1/2)
    '''
    res = np.zeros(lat.shape)/u.km/u.km
    res[np.abs(lat)<i] = 1/(2*np.pi**2 * (const.R_earth+hsat)**2 * (np.sin(i)**2-np.sin(lat[np.abs(lat)<i])**2)**.5)
    return res

def satellite_density(Nsat, lat, i, hsat, d, cosalpha):
    '''equation 3
    rho_sat = Nsat * P(lat, i, hsat) * (d^2 A)/(cos alpha), with A=1 here
    '''
    return Nsat*single_sat_density(lat, i, hsat)*d**2/cosalpha #*A

def compute_vdirection(Omega, i, lat, lon):
    ''' equation A.8 to A.11
    '''
    #unit vector from earth center to ascending node
    CA = np.stack([np.cos(Omega.to(u.rad)), np.sin(Omega.to(u.rad)), np.zeros(Omega.shape)])
    #unit vector from earth center to perigee
    CP = np.stack([np.cos(Omega.to(u.rad)+(np.pi/2)*u.rad)*np.cos(i.to(u.rad)), 
                   np.sin(Omega.to(u.rad)+(np.pi/2)*u.rad)*np.cos(i.to(u.rad)), 
                   np.sin(i.to(u.rad))*np.ones(Omega.shape)])
    #unit vector from earth center to satellite position
    CS = np.array([np.cos(lon.to(u.rad))*np.cos(lat.to(u.rad)),
                   np.sin(lon.to(u.rad))*np.cos(lat.to(u.rad)),
                   np.sin(lat.to(u.rad))*np.ones(Omega.shape)])
    # compute the velocity direction
    N = np.cross(CA, CP, axisa=0, axisb=0, axisc=0)
    NCS = np.cross(N, CS, axisa=0, axisb=0, axisc=0)
    return NCS

def compute_geocentric_vs(i, lat, lon, hsat):
    '''equation A.6, A.7 and A.11
    '''
    # compute the orbital velocity norm
    v_norm = np.sqrt((const.G*const.M_earth)/((const.R_earth+hsat).to(u.m))).to(u.km/u.s)

    # compute the velocity direction
    sin_lambda = np.tan(lat.to(u.rad))/np.tan(i)
    Omega_N = lon.to(u.rad) - np.arcsin(sin_lambda)
    Omega_S = lon.to(u.rad) + np.arcsin(sin_lambda) + np.pi*u.rad
    NCS_N = compute_vdirection(Omega_N, i, lat, lon)
    NCS_S = compute_vdirection(Omega_S, i, lat, lon)
    return NCS_N*v_norm, NCS_S*v_norm

def compute_topocentric_v(v, obs_lat, obs_lon):
    '''equation A.12
    '''
    v_obs = w_earth*const.R_earth*np.cos(obs_lat.to(u.rad))*np.array([-np.sin(obs_lon.to(u.rad)), np.cos(obs_lon.to(u.rad)), 0])
    #print(v.shape, v.unit, v_obs.unit)
    return v - v_obs[:,np.newaxis,np.newaxis]/u.rad

def compute_apparent_v(v, target_dec, target_ra):
    '''equation A.14
    '''
    
    OS = np.stack([np.cos(target_dec.to(u.rad))*np.cos(target_ra.to(u.rad)),
                   np.cos(target_dec.to(u.rad))*np.sin(target_ra.to(u.rad)),
                   np.sin(target_dec.to(u.rad))*np.ones((target_ra.shape))])
    v_dot_OS = np.einsum('ijk,ijk->jk', v, OS)
    
    OS_dot_OS = np.einsum('ijk,ijk->jk', OS, OS)
    
    
    v_los = (v_dot_OS / OS_dot_OS)[None, ...] * OS
    v_perp = v - v_los
    return v_perp

def compute_wsat(i, lat, lon, hsat, obs_lat, obs_lon, target_dec, target_ra):
    '''equation A.15
    '''
    nra, ndec = target_ra.size, target_dec.size
    target_ra = np.tile(target_ra[:,np.newaxis], (1, ndec))
    target_dec = np.tile(target_dec[np.newaxis,:], (nra, 1))
    
    V_N, V_S = compute_geocentric_vs(i, lat, lon, hsat)
    topo_V_N = compute_topocentric_v(V_N, obs_lat, obs_lon)
    topo_V_S = compute_topocentric_v(V_S, obs_lat, obs_lon)
    app_V_N = compute_apparent_v(topo_V_N, target_dec, target_ra)
    app_V_S = compute_apparent_v(topo_V_S, target_dec, target_ra)
    norm_V_N = np.sqrt(np.einsum('ijk,ijk->jk',app_V_N,app_V_N))
    norm_V_S = np.sqrt(np.einsum('ijk,ijk->jk',app_V_S,app_V_S))
    mean_app_V = (norm_V_N+norm_V_S)/2.
    return mean_app_V

def compute_d_phi(obs_phi, hsat, target_dec, target_ra):
    '''equation A.24
    d^2 + 2 Re (xos cos(obs_lat) + z sin(obs_lat)) - (hsat^2 + 2 Re hsat) = 0
    '''
    target_ra = target_ra[:,np.newaxis]
    target_dec = target_dec[np.newaxis,:]
    b = 2*const.R_earth*(np.cos(target_dec)*np.cos(target_ra)*np.cos(obs_phi) + np.sin(target_dec)*np.sin(obs_phi))
    c = -(hsat**2+2*hsat*const.R_earth)
    Delta = b**2-4*c #a=1
    r1, r2 = (-b-np.sqrt(Delta))/2, (-b+np.sqrt(Delta))/2.
    lat = np.arcsin((r2*np.sin(target_dec)+const.R_earth*np.sin(obs_phi))/(const.R_earth+hsat))
    sintheta = (r2*np.cos(target_dec)*np.sin(target_ra))/(np.cos(lat)*(const.R_earth+hsat))
    costheta = (r2*np.cos(target_dec)*np.cos(target_ra)+const.R_earth*np.cos(obs_phi))/(np.cos(lat)*(const.R_earth+hsat))
    lon = np.arctan2(sintheta, costheta)
    return r2, lat, lon

def compute_cosalpha(d, hsat):
    '''equation A.25
    cos alpha = (Rsat^2 + d^2 - Re^2)/(2 d Rsat)
    '''
    Rsat = const.R_earth+hsat
    return (Rsat**2+d**2-const.R_earth**2)/(2*d*Rsat)

def compute_nsats(rho_sat, Lfov, wsat, texp):
    '''equation 1
    Ntrail = rho_sat (pi R^2 + L wsat Tobs)
    '''
    return rho_sat*(np.pi*(Lfov/2)**2 + Lfov*wsat.to(u.deg/u.s)*texp).to(u.rad*u.rad)


def get_highest_time(start_day, obslat, obslon, obsheight, target_str):
    """ obslat and obslon in degrees. obsheight in meters. """
    end_day = start_day + 1.0
    times = Time(np.linspace(start_day, end_day, 360), format='jd', scale='utc')
    observer = EarthLocation(lat=obslat, lon=obslon, height=obsheight)
    target = SkyCoord(target_str, unit=(u.hourangle, u.deg))
    alt=[]
    for t in times:
        targetaltaz = target.transform_to(AltAz(obstime=t, location=observer))
        alt.append(targetaltaz.alt.to(u.deg).value)
    highest_idx = np.argmax(np.array(alt))
    return times[highest_idx], target.ra.deg, target.dec.deg #This is the JD of the center of the observing window

def compute_shell_satellite_density(obslat, obslon, i, Nsat, hsat, target_dec, target_ra, Lfov, obslength):
    '''
    Compute the satellite density for a given shell.

    Parameters
    ----------
    obslat : astropy Quantity
        Latitude of the observer (in rad).
    obslon : astropy Quantity
        Longitude of the observer (in rad).
    i : astropy Quantity
        Inclination of the shell (in rad).
    Nsat : int
        Number of satellites in the shell.
    hsat : astropy Quantity
        Altitude of the shell (in km).
    target_dec : astropy Quantity
        Declination of the target (in rad).
    target_ra : astropy Quantity
        Right ascension of the target (in rad).
    Lfov : astropy Quantity
        Field of view of the telescope (in rad).
    obslength : astropy Quantity
        Length of the observation (in s).

    Returns
    -------
    nsats : astropy Quantity
        Number of satellites expected at the target dec and ra during the observation.
    '''
    d, lat, lon = compute_d_phi(obslat, hsat, target_dec, target_ra)
    cosalpha = compute_cosalpha(d, hsat)
    P_sat = single_sat_density(lat, i, hsat) #
    rho_sat = satellite_density(Nsat, lat, i, hsat, d, cosalpha)
    wsat = compute_wsat(i, lat, lon, hsat, obslat, obslon, target_dec, target_ra)
    nsats = np.nan_to_num(compute_nsats(rho_sat, Lfov, wsat/d.to(u.m)*u.rad, obslength))
    return nsats.to(u.rad*u.rad)

def compute_total_satellite_density(obslat, obslon, shells, target_dec, target_ra, Lfov, obslength):
    '''
    Compute the total satellite density for a list of shells.

    Parameters
    ----------
    obslat : astropy Quantity
        Latitude of the observer (in rad).
    obslon : astropy Quantity
        Longitude of the observer (in rad).
    shells: pandas DataFrame
        DataFrame containing the shell parameters (inclination in deg, altitude in km, number of satellites).
    target_dec : astropy Quantity
        Declination of the target (in rad).
    target_ra : astropy Quantity
        Right ascension of the target (in rad).
    Lfov : astropy Quantity
        Field of view of the telescope (in rad).
    obslength : astropy Quantity
        Length of the observation (in s).

    Returns
    -------
    total_nsats : astropy Quantity
        Total number of satellites expected at the target dec and ra during the observation.
    '''
    total_nsats = 0
    nshells = shells.shape[0]
    for idx in range(nshells):
        i = (shells['i'][idx]*u.deg).to(u.rad)
        Nsat = shells['n'][idx]
        hsat = shells['h'][idx] * u.km
        nsats = compute_shell_satellite_density(obslat, obslon, i, Nsat, hsat, target_dec, target_ra, Lfov, obslength)
        total_nsats += nsats
    return total_nsats


def run_analytical_ska(obslat, obslon, obstime, obslength, target_dec, target_ra, shells_i, shells_h, shells_n, Lfov, nstat=100):
    location = EarthLocation(lat=obslat, lon=obslon)
    altaz_frame = AltAz(obstime=obstime, location=location)
    LST = obstime.sidereal_time('apparent', longitude=location.lon)
    target_dec = target_dec.to(u.rad)
    target_ra = target_ra.to(u.rad) - LST.to(u.rad)
    all_inits, all_ts = [],[]
    result=0
    nshells = len(shells_i)
    for idx in range(nshells):
        i = (shells_i[idx]*u.deg).to(u.rad)
        Nsat = shells_n[idx]
        hsat = shells_h[idx] * u.km
        d, lat, lon = compute_d_phi(obslat, hsat, target_dec, target_ra)
        cosalpha = compute_cosalpha(d, hsat)
        P_sat = single_sat_density(lat, i, hsat) #
        rho_sat = satellite_density(Nsat, lat, i, hsat, d, cosalpha)
        wsat = compute_wsat(i, lat, lon, hsat, obslat, obslon, target_dec, target_ra)
        
        nsats = np.nan_to_num(compute_nsats(rho_sat, Lfov, wsat/d.to(u.m)*u.rad, obslength))
        n_sample = np.random.poisson(nsats.value.squeeze(), size=(nstat,)).squeeze()
        if n_sample.any()>0:
            tmp_inits, tmp_ts = draw_passes(n_sample, Lfov, wsat, d, obslength)
            all_inits.append(tmp_inits)
            all_ts.append(tmp_ts.value)
    all_inits = np.concatenate(all_inits, axis=1)
    all_ts = np.concatenate(all_ts, axis=1)
    return all_inits, all_ts

def draw_passes(n_samp, Lfov, wsat, d, obslength):
    nstat = n_samp.size
    nsampmax = n_samp.max()
    hs = np.random.uniform(low=-Lfov.value/2, high=Lfov.value/2, size=(nstat, nsampmax))
    inits = np.random.uniform(low=0., high=obslength.value, size=(nstat, nsampmax))
    ts = 2*np.cos(np.arcsin(hs/(Lfov.value/2)))*(Lfov.to(u.rad)/2)/(wsat/d.to(u.m)*u.rad)
    for j in range(nstat):
        ts[j, n_samp[j]:] = -1 *u.s
    return inits, ts

def compute_starlink_fraction(obslength, ntimestep, all_inits, all_ts):
    timeframe = np.linspace(0,obslength.to(u.s).value,ntimestep)
    nsat_at_t = np.zeros((all_inits.shape[0], ntimestep))
    #print(nsat_at_t.shape, all_inits.shape, timeframe.shape)
    for j in range(all_inits.shape[0]):
        for init, tpass in zip(all_inits[j], all_ts[j]):
            nsat_at_t[j, (timeframe>init)*(timeframe<init+tpass)] += 1
    return np.mean(nsat_at_t>=1, axis=1)