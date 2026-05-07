import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
import pandas as pd
from astropy.units import Quantity
from functools import cached_property

from .model import satellite_density, compute_d_phi, compute_cosalpha, compute_wsat, compute_nsats, draw_passes

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

    def sample_passes(self, nstat):
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