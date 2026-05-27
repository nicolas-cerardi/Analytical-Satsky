import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.coordinates import EarthLocation
import pandas as pd
from astropy.units import Quantity
from functools import cached_property

from .model import satellite_density, compute_d_phi, compute_cosalpha, compute_wsat, compute_nsats, draw_passes, compute_geocentric_vs
from .large_fov_model import angular_distance, compute_d_lat_lon_gcrs, pointing_to_radec_circle, compute_flux_vel, local_basis_at_radec

class SingleShellObs:
    """
    Cached model for the contribution of a single satellite shell.

    This class wraps the core analytical functions used to describe the
    geometry, apparent motion, and projected density of satellites belonging
    to one orbital shell. Derived quantities are evaluated lazily using
    ``cached_property`` in order to avoid recomputing intermediate quantities
    shared by multiple calculations.

    Parameters
    ----------
    obsloc : astropy.coordinates.EarthLocation
        Observer location on Earth.
    shell : pandas.Series
        Row of the shell catalogue describing a single orbital shell.
        The following fields are expected:

        - ``"i"`` : orbital inclination in degrees
        - ``"h"`` : shell altitude in km
        - ``"n"`` : number of satellites in the shell

    target_dec : astropy.units.Quantity
        Declination of the target line of sight. Must be angular.
    target_lha : astropy.units.Quantity
        Local hour angle of the target line of sight. Must be angular.
    Lfov : astropy.units.Quantity
        Effective telescope field of view. Must be angular.
    tobs : astropy.units.Quantity
        Observation duration. Must be convertible to seconds.

    Attributes
    ----------
    obsloc : astropy.coordinates.EarthLocation
        Observer location.
    shell : pandas.Series
        Shell description.
    target_dec : astropy.units.Quantity
        Target declination in radians.
    target_lha : astropy.units.Quantity
        Target local hour angle in radians.
    Lfov : astropy.units.Quantity
        Field of view in radians.
    tobs : astropy.units.Quantity
        Observation duration.
    i : astropy.units.Quantity
        Orbital inclination in radians.
    Nsat : int
        Number of satellites in the shell.
    hsat : astropy.units.Quantity
        Shell altitude.
    obslat : astropy.units.Quantity
        Observer latitude in radians.

    Notes
    -----
    Most derived quantities are implemented as ``cached_property`` objects.
    They are therefore computed only once upon first access and then stored
    on the instance for reuse.
    """

    def __init__(
        self,
        obsloc: EarthLocation,
        shell: pd.Series,
        target_dec: Quantity,
        target_lha: Quantity,
        Lfov: Quantity,
        tobs: Quantity,
    ) -> None:

        self.obsloc = obsloc
        self.shell = shell
        self.target_dec = target_dec.to(u.rad)
        self.target_lha = target_lha.to(u.rad)
        self.Lfov = Lfov.to(u.rad)
        self.tobs = tobs

        self.i = (shell["i"] * u.deg).to(u.rad)
        self.Nsat = shell["n"]
        self.hsat = shell["h"] * u.km

        self.obslat = np.arcsin(
            obsloc.z / np.sqrt(obsloc.x**2 + obsloc.y**2 + obsloc.z**2)
        )
    @cached_property
    def geometry(self):
        return compute_d_phi(
            self.obslat,
            self.hsat,
            self.target_dec,
            self.target_lha,
        )

    @cached_property
    def d(self):
        return self.geometry[0]

    @cached_property
    def lat(self):
        return self.geometry[1]

    @cached_property
    def lon(self):
        return self.geometry[2]

    @cached_property
    def cosalpha(self):
        return compute_cosalpha(self.d, self.hsat)

    @cached_property
    def rho_sat(self):
        return satellite_density(
            self.Nsat,
            self.lat,
            self.i,
            self.hsat,
            self.d,
            self.cosalpha,
        )

    @cached_property
    def wsat(self):
        return compute_wsat(
            self.i,
            self.lat,
            self.lon,
            self.hsat,
            self.obslat,
            self.target_dec,
            self.target_lha,
        )

    @cached_property
    def nsats(self):
        return np.nan_to_num(
            compute_nsats(
                self.rho_sat,
                self.Lfov,
                self.wsat / self.d * u.rad,
                self.tobs,
            )
        )
    
class MultiShellObs:
    """
    Cached analytical model for multiple orbital shells.

    One ``SingleShellObs`` instance is created internally for each row
    of the input shell catalogue.

    Parameters
    ----------
    obsloc : astropy.coordinates.EarthLocation
        Observer location on Earth.

    shells_df : pandas.DataFrame
        Table describing the orbital shells of the constellation.

    target_dec : astropy.units.Quantity
        Declination of the target line of sight. Must be angular.

    target_lha : astropy.units.Quantity
        Local hour angle of the target line of sight. Must be angular.

    Lfov : astropy.units.Quantity
        Effective telescope field of view. Must be angular.

    tobs : astropy.units.Quantity
        Observation duration. Must be convertible to seconds.

    Attributes
    ----------
    obsloc : astropy.coordinates.EarthLocation
        Observer location.

    shells_df : pandas.DataFrame
        Input shell catalogue.

    target_dec : astropy.units.Quantity
        Target declination in radians.

    target_lha : astropy.units.Quantity
        Target local hour angle in radians.

    Lfov : astropy.units.Quantity
        Field of view in radians.

    tobs : astropy.units.Quantity
        Observation duration.

    shell_models : list[SingleShellObs]
        List of ``SingleShellObs`` models, one for each orbital shell.

    nsats_per_shell : list
        Expected number of satellite crossings for each shell.

    total_satellite_density : astropy.units.Quantity
        Total projected satellite density summed over all shells.

    total_nsats : astropy.units.Quantity
        Total expected number of satellites crossing the field of view
        during the observation.

    Methods
    -------
    sample_passes(nstat)
        Draw stochastic satellite crossing events from the analytical
        model.
    """
    def __init__(self, obsloc, shells_df, target_dec, target_lha, Lfov, tobs):
        self.obsloc = obsloc
        self.shells_df = shells_df
        self.target_dec = target_dec.to(u.rad)
        self.target_lha = target_lha.to(u.rad)
        self.Lfov = Lfov.to(u.rad)
        self.tobs = tobs

        self.shell_models = [
            SingleShellObs(
                obsloc=self.obsloc,
                shell=shell,
                target_dec=self.target_dec,
                target_lha=self.target_lha,
                Lfov=self.Lfov,
                tobs=self.tobs,
            )
            for _, shell in self.shells_df.iterrows()
        ]

    @cached_property
    def nsats_per_shell(self):
        return [shell_model.nsats for shell_model in self.shell_models]

    @cached_property
    def total_satellite_density(self):
        return sum(shell_model.rho_sat for shell_model in self.shell_models)

    @cached_property
    def total_nsats(self):
        return sum(self.nsats_per_shell)

    def sample_passes(
        self,
        nstat: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Draw stochastic satellite crossing events from the analytical model.

        Parameters
        ----------
        nstat : int
            Number of independent statistical realisations to generate.

        Returns
        -------
        all_inits : numpy.ndarray
            Ingress times into the effective beam for all sampled events.

            The returned array has shape ``(nstat, nevents)``, where
            ``nevents`` is the total number of sampled events aggregated
            over all orbital shells.

            Values are expressed in seconds.

        all_ts : numpy.ndarray
            Fly-through durations across the effective beam for all sampled
            events.

            The returned array has the same shape as ``all_inits``.

            Values are expressed in seconds.

        Notes
        -----
        The number of sampled events varies between statistical realisations
        because crossings are drawn from Poisson statistics independently
        for each shell.

        This method combines the contributions from all orbital shells
        included in the ``MultiShellObs`` instance.

        The returned arrays are intended to be used with
        ``compute_occupancy_fraction()``.
        """
        all_inits = []
        all_ts = []

        for shell_model in self.shell_models:
            n_sample = np.random.poisson(
                shell_model.nsats.value.squeeze(),
                size=nstat,
            ).squeeze()

            if np.any(n_sample > 0):
                tmp_inits, tmp_ts = draw_passes(
                    n_sample,
                    shell_model.Lfov,
                    shell_model.wsat,
                    shell_model.d,
                    shell_model.tobs,
                )

                all_inits.append(tmp_inits)
                all_ts.append(tmp_ts.to(u.s).value)
        all_inits = np.concatenate(all_inits, axis=1)
        all_ts = np.concatenate(all_ts, axis=1)
        return all_inits, all_ts
    
class SingleShellFoV:
    def __init__(
        self,
        obs,
        shell,
        ndec=None,
        nra=None,
        Lfov=None,
        t_mjd=None,
        ra_grid=None,
        dec_grid=None,
        FoV_mask=None,
        obsloc=None,
        pointing_ra=None,
        pointing_dec=None,
        dOmega=None
    ):
        self.obs = obs
        if obsloc is not None:
            self.obsloc = obsloc
        else:
            self.obsloc = EarthLocation.from_geocentric(*obs.ITRF.compute()[0]*u.m)

        if pointing_ra is not None and pointing_dec is not None:
            self.pointing_ra = pointing_ra
            self.pointing_dec = pointing_dec
        else:
            self.pointing_ra = (obs.ra.compute()*u.deg).to(u.rad)
            self.pointing_dec = (obs.dec.compute()*u.deg).to(u.rad)

        self.shell = shell
        self.ndec = ndec
        self.nra = nra
        self.Lfov = Lfov #be sure that this is in radians
        self.t_mjd = t_mjd

        self.i = (shell["i"] * u.deg).to(u.rad)
        self.Nsat = shell["n"]
        self.hsat = shell["h"] * u.km
        self.obslat = np.arcsin(
            self.obsloc.z / np.sqrt(self.obsloc.x**2 + self.obsloc.y**2 + self.obsloc.z**2)
        )

        if ra_grid is None and dec_grid is None:
            if ndec is None or nra is None:
                raise ValueError("Needs ndec/nra if ra_grid/dec_grid are not provided.")
            range_dec = np.linspace(-np.pi/2, np.pi/2, ndec, endpoint=False) * u.rad
            range_ra = np.linspace(0, 2*np.pi, nra, endpoint=False) * u.rad
            ra_grid, dec_grid = np.meshgrid(range_ra, range_dec, indexing="ij")
            angdist = angular_distance(self.pointing_ra, self.pointing_dec, ra_grid, dec_grid)
            FoV_mask = angdist < self.Lfov / 2
            ddec = np.pi/ndec * u.rad
            dra = 2*np.pi/nra * np.cos(dec_grid[FoV_mask]) * u.rad
            dOmega = dra * ddec
        else:
            assert FoV_mask is not None and dOmega is not None, "If ra_grid and dec_grid are provided, FoV_mask and dOmega must also be provided."


        self.ra_grid = ra_grid
        self.dec_grid = dec_grid
        self.FoV_mask = FoV_mask
        self.dOmega = dOmega

    @cached_property
    def geometry(self):
        return compute_d_lat_lon_gcrs(
            self.obsloc,
            self.dec_grid[self.FoV_mask],
            self.ra_grid[self.FoV_mask],
            self.hsat,
            self.t_mjd
        )

    @cached_property
    def d(self):
        return self.geometry[0]

    @cached_property
    def lat(self):
        return self.geometry[1]

    @cached_property
    def lon(self):
        return self.geometry[2]

    @cached_property
    def cosalpha(self):
        return compute_cosalpha(self.d, self.hsat)

    @cached_property
    def rho_sat(self):
        return satellite_density(
            self.Nsat,
            self.lat,
            self.i,
            self.hsat,
            self.d,
            self.cosalpha,
        )
    
    @cached_property
    def nsat_shell_obs(self):
        return np.nan_to_num(self.rho_sat.decompose() * self.dOmega.to(u.rad**2))
    
    def sample_satellites(self):
        satellite_catalogue = pd.DataFrame(columns=['N', 'i', 'h', 'long_asc_node', 'periapsis', "decstart", "rastart", "tstart", "towards"])

        Xsat_shell_obs_N = np.random.poisson(self.nsat_shell_obs/2)
        for icoord, xsat in enumerate(Xsat_shell_obs_N):
            if xsat > 0:
                satellite_catalogue = self.compute_orbital_params_samples(satellite_catalogue, icoord, xsat, towards_north=True)

        Xsat_shell_obs_S = np.random.poisson(self.nsat_shell_obs/2)
        for icoord, xsat in enumerate(Xsat_shell_obs_S):
            if xsat > 0:
                satellite_catalogue = self.compute_orbital_params_samples(satellite_catalogue, icoord, xsat, towards_north=False)
        return satellite_catalogue
    
    def compute_orbital_params_samples(self, satellite_catalogue, icoord, xsat, towards_north=True):
        # compute Omega and periapsis for each sample in the catalogue
        lonlambda = np.arcsin(np.tan(self.lat[icoord])/np.tan(self.i))

        if towards_north:
            Omega = (self.lon[icoord] - lonlambda) #should be naturally in rad
            periapsis = np.arccos(np.cos((lonlambda).to(u.rad))*np.cos(self.lat[icoord])).to(u.deg)*np.sign(self.lat[icoord])     
        else:
            Omega = (self.lon[icoord] + lonlambda) + np.pi*u.rad # in rad
            periapsis = np.pi*u.rad - np.arccos(np.cos((lonlambda).to(u.rad))*np.cos(self.lat[icoord])).to(u.deg)*np.sign(self.lat[icoord])
        #omega = np.sqrt(const.GM_earth/(const.R_earth + self.h_sat)**3).to(1/u.s)*u.rad
        #periapsis -= (self.t_offset * omega)
        satellite_catalogue.loc[-1] = [xsat, self.i.to(u.deg).value, self.hsat.to(u.km).value, Omega.to(u.deg).value, periapsis.to(u.deg).value, self.dec_grid[self.FoV_mask][icoord], self.ra_grid[self.FoV_mask][icoord], 0., ("N" if towards_north else "S")]
        satellite_catalogue.index = satellite_catalogue.index + 1
        return satellite_catalogue.sort_index()
    
    
class MultiShellFoV:
    def __init__(self, obs, shells_df, ndec, nra, Lfov, t_mjd):
        self.obs = obs
        self.ndec = ndec
        self.nra = nra
        self.Lfov = Lfov.to(u.rad)
        self.t_mjd = t_mjd

        self.obsloc = EarthLocation.from_geocentric(*obs.ITRF.compute()[0] * u.m)
        self.pointing_ra = (obs.ra.compute() * u.deg).to(u.rad)
        self.pointing_dec = (obs.dec.compute() * u.deg).to(u.rad)

        range_dec = np.linspace(-np.pi/2, np.pi/2, ndec, endpoint=False) * u.rad
        range_ra = np.linspace(0, 2*np.pi, nra, endpoint=False) * u.rad
        self.ra_grid, self.dec_grid = np.meshgrid(range_ra, range_dec, indexing="ij")

        angdist = angular_distance(
            self.pointing_ra,
            self.pointing_dec,
            self.ra_grid,
            self.dec_grid,
        )
        self.FoV_mask = angdist < self.Lfov / 2

        ddec = np.pi/ndec * u.rad
        dra = 2*np.pi/nra * np.cos(self.dec_grid[self.FoV_mask]) * u.rad
        self.dOmega = dra * ddec    

        self.shell_models = [
            SingleShellFoV(
                obs=obs,
                shell=shell,
                ndec=ndec,
                nra=nra,
                Lfov=self.Lfov,
                t_mjd=t_mjd,
                ra_grid=self.ra_grid,
                dec_grid=self.dec_grid,
                FoV_mask=self.FoV_mask,
                obsloc=self.obsloc,
                pointing_ra=self.pointing_ra,
                pointing_dec=self.pointing_dec,
                dOmega=self.dOmega
            )
            for _, shell in shells_df.iterrows()
        ]

    @cached_property
    def rho_sat(self):
        return sum(np.nan_to_num(shell_model.rho_sat) for shell_model in self.shell_models)
    
    def sample_satellites(self):
        satellite_catalogue = pd.DataFrame(columns=['N', 'i', 'h', 'long_asc_node', 'periapsis', "decstart", "rastart", "tstart", "towards"])
        for shell_model in self.shell_models:
            satellite_catalogue = pd.concat([shell_model.sample_satellites(), satellite_catalogue], ignore_index=True)
        return satellite_catalogue.sort_values(by="tstart").reset_index(drop=True)
    
    

class SingleShellFlux:
    def __init__(
        self,
        obs,
        shell,
        Npoints=None,
        Lfov=None,
        t_mjd=None,
        t_init_mjd=None,
        dt=None,
        obsloc=None,
        pointing_ra=None,
        pointing_dec=None,
        circle_vecs=None,
        circle_ra=None,
        circle_dec=None,
    ):
        self.obs = obs
        self.obsloc = (
            obsloc
            if obsloc is not None
            else EarthLocation.from_geocentric(*obs.ITRF.compute()[0] * u.m)
        )
        if pointing_ra is not None and pointing_dec is not None:
            self.pointing_ra = pointing_ra
            self.pointing_dec = pointing_dec
        else:
            self.pointing_ra = (obs.pointing_ra * u.deg).value
            self.pointing_dec = (obs.pointing_dec * u.deg).value

        
        self.Npoints = Npoints
        self.Lfov = Lfov #be sure that this is in radians

        if circle_vecs is None or circle_ra is None or circle_dec is None:
            if Npoints is None or Lfov is None:
                raise ValueError(
                    "Needs Npoints/Lfov if circle_vecs/circle_ra/circle_dec are not provided."
                )

            circle_vecs, circle_ra, circle_dec = pointing_to_radec_circle(
                self.pointing_ra,
                self.pointing_dec,
                self.Lfov,
                Npoints=self.Npoints,
            )
        
        self.circle_vecs = circle_vecs
        self.circle_ra = circle_ra * u.rad
        self.circle_dec = circle_dec * u.rad

        self.shell = shell
        self.i = (shell["i"] * u.deg).to(u.rad)
        self.Nsat = shell["n"]
        self.hsat = shell["h"] * u.km

        self.t_mjd = t_mjd #should be an astropy Time object in MJD format
        self.t_offset = (self.t_mjd - t_init_mjd).to(u.s)
        self.dt = dt #should be in seconds

    @cached_property
    def geometry(self):
        return compute_d_lat_lon_gcrs(
            self.obsloc,
            self.circle_dec,
            self.circle_ra,
            self.hsat,
            self.t_mjd
        )

    @cached_property
    def d(self):
        return self.geometry[0]

    @cached_property
    def lat(self):
        return self.geometry[1]

    @cached_property
    def lon(self):
        return self.geometry[2]

    @cached_property
    def cosalpha(self):
        return compute_cosalpha(self.d, self.hsat)

    @cached_property
    def rho_sat(self):
        return satellite_density(
            self.Nsat,
            self.lat,
            self.i,
            self.hsat,
            self.d,
            self.cosalpha,
        )

    @cached_property
    def shell_nsat_fluxes(self):
        V_N, V_S = compute_geocentric_vs(self.i, self.lat, self.lon, self.hsat)

        e_east0, e_north0, e_center = local_basis_at_radec(self.pointing_ra, self.pointing_dec)
        den = self.circle_vecs @ e_center
        x = (self.circle_vecs @ e_east0) / den
        y = (self.circle_vecs @ e_north0) / den
        dl_norm = self.Lfov / 2.0 * 2.0 * np.pi / self.Npoints # boundary angular element length [rad]
        r_dot_E = self.circle_vecs @ e_east0
        r_dot_N = self.circle_vecs @ e_north0
        r_dot_C = self.circle_vecs @ e_center
        tx = np.roll(x, -1) - np.roll(x, 1)
        ty = np.roll(y, -1) - np.roll(y, 1)
        tnorm = np.hypot(tx, ty)
        tx = dl_norm * tx / tnorm
        ty = dl_norm * ty / tnorm
        # Keep only entering contributions
        _, flux_N, _, _ = compute_flux_vel(V_N, self.obsloc, self.circle_dec.to(u.deg).value, self.circle_ra.to(u.deg).value, self.d, 
                                e_east0, e_north0, e_center, r_dot_E, r_dot_N, r_dot_C, tx, ty, self.t_mjd, self.dt)
        _, flux_S, _, _ = compute_flux_vel(V_S, self.obsloc, self.circle_dec.to(u.deg).value, self.circle_ra.to(u.deg).value, self.d,
                                e_east0, e_north0, e_center, r_dot_E, r_dot_N, r_dot_C, tx, ty, self.t_mjd, self.dt)
        return self.rho_sat.squeeze()/2 * flux_N, self.rho_sat.squeeze()/2 * flux_S
    
    def sample_satellites(self):
        satellite_catalogue = pd.DataFrame(columns=['N', 'i', 'h', 'long_asc_node', 'periapsis', "decstart", "rastart", "tstart", "towards"])

        satflux_N, satflux_S = self.shell_nsat_fluxes
        Xsat_shell_obs_N = np.random.poisson(np.nan_to_num(satflux_N).decompose())
        for icoord, xsat in enumerate(Xsat_shell_obs_N):
            if xsat > 0:
                satellite_catalogue = self.compute_orbital_params_samples(satellite_catalogue, icoord, xsat, towards_north=True)

        Xsat_shell_obs_S = np.random.poisson(np.nan_to_num(satflux_S).decompose())
        for icoord, xsat in enumerate(Xsat_shell_obs_S):
            if xsat > 0:
                satellite_catalogue = self.compute_orbital_params_samples(satellite_catalogue, icoord, xsat, towards_north=False)
        return satellite_catalogue
    
    def compute_orbital_params_samples(self, satellite_catalogue, icoord, xsat, towards_north=True):
        # compute Omega and periapsis for each sample in the catalogue
        lonlambda = np.arcsin(np.tan(self.lat[icoord])/np.tan(self.i))

        if towards_north:
            Omega = (self.lon[icoord] - lonlambda) #should be naturally in rad
            periapsis = np.arccos(np.cos((lonlambda).to(u.rad))*np.cos(self.lat[icoord])).to(u.deg)*np.sign(self.lat[icoord])     
        else:
            Omega = (self.lon[icoord] + lonlambda) + np.pi*u.rad # in rad
            periapsis = np.pi*u.rad - np.arccos(np.cos((lonlambda).to(u.rad))*np.cos(self.lat[icoord])).to(u.deg)*np.sign(self.lat[icoord])
        omega = np.sqrt(const.GM_earth/(const.R_earth + self.hsat)**3).to(1/u.s)*u.rad
        periapsis -= (self.t_offset * omega)
        satellite_catalogue.loc[-1] = [xsat, self.i.to(u.deg).value, self.hsat.to(u.km).value, Omega.to(u.deg).value, periapsis.to(u.deg).value, self.circle_dec[icoord], self.circle_ra[icoord], self.t_offset.to(u.s).value, ("N" if towards_north else "S")]
        satellite_catalogue.index = satellite_catalogue.index + 1
        return satellite_catalogue.sort_index()

class MultiShellFlux:
    def __init__(self, obs, shells_df, Npoints, Lfov, t_mjd, t_init_mjd, dt):
        self.obs = obs
        self.Npoints = Npoints
        self.Lfov = Lfov
        self.t_mjd = t_mjd
        self.t_init_mjd = t_init_mjd
        self.dt = dt

        self.obsloc = EarthLocation.from_geocentric(*obs.ITRF.compute()[0] * u.m)

        self.pointing_ra = (obs.pointing_ra * u.deg).value
        self.pointing_dec = (obs.pointing_dec * u.deg).value

        self.circle_vecs, self.circle_ra, self.circle_dec = pointing_to_radec_circle(
            self.pointing_ra,
            self.pointing_dec,
            self.Lfov,
            Npoints=self.Npoints,
        )

        self.shell_models = [
            SingleShellFlux(
                obs=obs,
                shell=shell,
                Npoints=Npoints,
                Lfov=Lfov,
                t_mjd=t_mjd,
                t_init_mjd=t_init_mjd,
                dt=dt,
                obsloc=self.obsloc,
                pointing_ra=self.pointing_ra,
                pointing_dec=self.pointing_dec,
                circle_vecs=self.circle_vecs,
                circle_ra=self.circle_ra,
                circle_dec=self.circle_dec,
            )
            for _, shell in shells_df.iterrows()
        ]

    @cached_property
    def shell_nsat_fluxes(self):
        flux_N = []
        flux_S = []

        for shell_model in self.shell_models:
            shell_flux_N, shell_flux_S = shell_model.shell_nsat_fluxes
            flux_N.append(shell_flux_N)
            flux_S.append(shell_flux_S)

        return sum(flux_N), sum(flux_S)

    def sample_satellites(self):
        satellite_catalogue = pd.DataFrame(
            columns=[
                "N",
                "i",
                "h",
                "long_asc_node",
                "periapsis",
                "decstart",
                "rastart",
                "tstart",
                "towards",
            ]
        )

        for shell_model in self.shell_models:
            satellite_catalogue = pd.concat(
                [shell_model.sample_satellites(), satellite_catalogue],
                ignore_index=True,
            )

        return satellite_catalogue.sort_values(by="tstart").reset_index(drop=True)

class IntegralObsModel:
    pass