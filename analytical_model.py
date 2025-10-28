import numpy as np

import os
import scipy
import skyfield
from skyfield.positionlib import ICRF
import astropy.units as u
import astropy.constants as const

#R_earth = 6378 * u.km
#M_earth = 5.9722*1e24 *u.kg
#G = 6.6743 *1e-11*(u.m*u.m*u.m/u.kg/u.s/u.s)
w_earth = 7.2921159*1e-5 * u.rad / u.s
v_equator = 465.1 * u.m / u.s

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
    #print(Omega.shape, lat.shape, phi.shape)
    CA = np.stack([np.cos(Omega.to(u.rad)), np.sin(Omega.to(u.rad)), np.zeros(Omega.shape)]) #np.array([np.cos(Omega.to(u.rad)), np.sin(Omega.to(u.rad)), 0])
    #print(CA.shape)
    CP = np.stack([np.cos(Omega.to(u.rad)+(np.pi/2)*u.rad)*np.cos(i.to(u.rad)), 
                   np.sin(Omega.to(u.rad)+(np.pi/2)*u.rad)*np.cos(i.to(u.rad)), 
                   np.sin(i.to(u.rad))*np.ones(Omega.shape)])
    CS = np.array([np.cos(lon.to(u.rad))*np.cos(lat.to(u.rad)),
                   np.sin(lon.to(u.rad))*np.cos(lat.to(u.rad)),
                   np.sin(lat.to(u.rad))*np.ones(Omega.shape)])
    N = np.cross(CA, CP, axisa=0, axisb=0, axisc=0)
    NCS = np.cross(N, CS, axisa=0, axisb=0, axisc=0)
    #print(NCS.shape)
    return NCS

def compute_geocentric_vs(i, lat, lon, hsat):
    v_norm = np.sqrt((const.G*const.M_earth)/((const.R_earth+hsat).to(u.m)))
    
    sin_lambda = np.tan(lat.to(u.rad))/np.tan(i)
    Omega_N = lon.to(u.rad) - np.arcsin(sin_lambda)
    Omega_S = lon.to(u.rad) + np.arcsin(sin_lambda) + np.pi*u.rad
    NCS_N = compute_vdirection(Omega_N, i, lat, lon)
    NCS_S = compute_vdirection(Omega_S, i, lat, lon)
    return NCS_N*v_norm, NCS_S*v_norm

def compute_topocentric_v(v, obs_lat, obs_lon):
    v_obs = w_earth*const.R_earth*np.cos(obs_lat.to(u.rad))*np.array([-np.sin(obs_lon.to(u.rad)), np.cos(obs_lon.to(u.rad)), 0])
    #print(v.shape, v.unit, v_obs.unit)
    return v - v_obs[:,np.newaxis,np.newaxis]/u.rad

def compute_apparent_v(v, target_dec, target_ra):
    #print(target_dec.shape, target_ra.shape)
    OS = np.stack([np.cos(target_dec.to(u.rad))*np.cos(target_ra.to(u.rad)),
                   np.cos(target_dec.to(u.rad))*np.sin(target_ra.to(u.rad)),
                   np.sin(target_dec.to(u.rad))*np.ones((target_ra.shape))])
    v_dot_OS = np.einsum('ijk,ijk->jk', v, OS)
    #print('v_dot_OS', v_dot_OS.shape)
    OS_dot_OS = np.einsum('ijk,ijk->jk', OS, OS)
    
    #print('OS_dot_OS', OS_dot_OS.shape)
    v_los = (v_dot_OS / OS_dot_OS)[None, ...] * OS
    v_perp = v - v_los
    return v_perp

def compute_wsat(i, lat, lon, hsat, obs_lat, obs_lon, target_dec, target_ra):
    #print(target_ra.shape, target_dec.shape)
    nra, ndec = target_ra.size, target_dec.size
    target_ra = np.tile(target_ra[:,np.newaxis], (1, ndec))
    target_dec = np.tile(target_dec[np.newaxis,:], (nra, 1))
    #print(target_ra.shape, target_dec.shape)
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