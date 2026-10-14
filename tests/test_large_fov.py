"""
This testing file provides tests for the large-FoV / satellite-catalogue sampling path:
 - the shell_obs.py classes (``MultiShellFoV``, ``MultiShellFlux``, ``IntegralObsModel``)
 - the pure geometry functions in large_fov_model.py they are built on
"""

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time

from analytical_satsky import MultiShellFoV, MultiShellFlux, IntegralObsModel
from analytical_satsky.large_fov_model import (
    radec_to_unitvec,
    local_basis_at_radec,
    angular_distance,
    pointing_to_radec_circle,
    satellite_positions_gcrs,
    positions_to_lm,
)

CATALOGUE_COLUMNS = [
    "N", "i", "h", "long_asc_node", "periapsis", "decstart", "rastart", "tstart", "towards",
]


class _DummyComputable:
    """Mimics a dask-delayed-like object exposing ``.compute()``, as used by
    tabsim's ``Observation.ITRF``/``.ra``/``.dec``."""

    def __init__(self, value):
        self._value = value

    def compute(self):
        return self._value


class _DummyObs:
    """Minimal stand-in for a tabsim ``Observation``, exposing only the
    ``.ITRF``/``.ra``/``.dec`` attributes that ``MultiShellFoV``/``MultiShellFlux``/
    ``IntegralObsModel`` read from it."""

    def __init__(self, obsloc: EarthLocation, ra_deg: float, dec_deg: float):
        itrf_m = np.array([
            obsloc.x.to_value(u.m),
            obsloc.y.to_value(u.m),
            obsloc.z.to_value(u.m),
        ])
        self.ITRF = _DummyComputable(itrf_m[None, :])  # shape (1, 3), like obs.ITRF.compute()
        self.ra = _DummyComputable(ra_deg)
        self.dec = _DummyComputable(dec_deg)


def _make_dummy_obs():
    obsloc = EarthLocation(lat=30 * u.deg, lon=0 * u.deg, height=0 * u.m)
    return _DummyObs(obsloc, ra_deg=83.0, dec_deg=-30.0)


def _assert_empty_catalogue(catalogue):
    assert len(catalogue) == 0
    assert list(catalogue.columns) == CATALOGUE_COLUMNS


def test_multi_shell_fov_sample_satellites_empty():
    # A shell with n=0 satellites always samples zero crossings, so
    # nsat_shell_obs is exactly zero everywhere and rng.poisson(0) is
    # deterministically 0 regardless of seed.
    obs = _make_dummy_obs()
    shells_df = pd.DataFrame({"n": [0], "h": [550.0], "i": [53.0]})

    multi_shell_fov = MultiShellFoV(
        obs=obs,
        shells_df=shells_df,
        ndec=20,
        nra=40,
        Lfov=10.0 * u.deg,
        t_mjd=Time("2026-01-01T00:00:00"),
    )

    catalogue = multi_shell_fov.sample_satellites(rng=np.random.default_rng(42))
    _assert_empty_catalogue(catalogue)


def test_multi_shell_flux_sample_satellites_empty():
    obs = _make_dummy_obs()
    shells_df = pd.DataFrame({"n": [0], "h": [550.0], "i": [53.0]})
    t_init_mjd = Time("2026-01-01T00:00:00")

    multi_shell_flux = MultiShellFlux(
        obs=obs,
        shells_df=shells_df,
        Npoints=50,
        Lfov=10.0 * u.deg,
        t_mjd=t_init_mjd,
        t_init_mjd=t_init_mjd,
        dt=1.0 * u.s,
    )

    catalogue = multi_shell_flux.sample_satellites(rng=np.random.default_rng(42))
    _assert_empty_catalogue(catalogue)


def test_integral_obs_model_sample_satellites_empty():
    obs = _make_dummy_obs()
    shells_df = pd.DataFrame({"n": [0], "h": [550.0], "i": [53.0]})
    t_init_mjd = Time("2026-01-01T00:00:00")

    integral_obs_model = IntegralObsModel(
        obs=obs,
        shells_df=shells_df,
        Lfov=10.0 * u.deg,
        t_exp_mjd=0.01,  # short exposure, keeps the per-timestep loop fast
        t_init_mjd=t_init_mjd,
        ndec=20,
        nra=40,
        Npoints=50,
        dt=60.0 * u.s,
    )

    catalogue = integral_obs_model.sample_satellites(seed=42)
    _assert_empty_catalogue(catalogue)


def test_integral_obs_model_sample_satellites_regression_snapshot():
    # Row count and per-row values for a fixed seed/config, frozen as a
    # stability guard against unintended changes to the sampling pipeline
    # (initial FoV snapshot + time-stepped flux sampling, orbital element
    # computation, rng draws). These values were computed with the current
    # implementation and are not independently derived by hand.
    obs = _make_dummy_obs()
    shells_df = pd.DataFrame({"n": [1000.0], "h": [550.0], "i": [53.0]})
    t_init_mjd = Time("2026-01-01T00:00:00")

    integral_obs_model = IntegralObsModel(
        obs=obs,
        shells_df=shells_df,
        Lfov=10.0 * u.deg,
        t_exp_mjd=0.01,  # short exposure, keeps the per-timestep loop fast
        t_init_mjd=t_init_mjd,
        ndec=20,
        nra=40,
        Npoints=50,
        dt=60.0 * u.s,
    )

    catalogue = integral_obs_model.sample_satellites(seed=42)

    assert len(catalogue) == 3

    expected_i = [53.0, 53.0, 53.0]
    expected_h = [550.0, 550.0, 550.0]
    expected_long_asc_node = [81.003180, 83.917466, 83.361406]
    expected_periapsis = [-5.338444, -27.675629, -27.532380]
    expected_decstart_rad = [-0.5323475012959175, -0.6106834267025876, -0.5965919929786228]
    expected_rastart_rad = [1.3480981020846348, 1.4419433244499613, 1.392141039112944]
    expected_tstart = [531.6923076923097, 864.0000000000008, 864.0000000000008]
    expected_towards = ["N", "N", "N"]

    assert np.allclose(catalogue["i"], expected_i)
    assert np.allclose(catalogue["h"], expected_h)
    assert np.allclose(catalogue["long_asc_node"], expected_long_asc_node, rtol=1e-6)
    assert np.allclose(catalogue["periapsis"], expected_periapsis, rtol=1e-6)
    assert np.allclose(
        [q.to_value(u.rad) for q in catalogue["decstart"]], expected_decstart_rad, rtol=1e-8
    )
    assert np.allclose(
        [q.to_value(u.rad) for q in catalogue["rastart"]], expected_rastart_rad, rtol=1e-8
    )
    assert np.allclose(catalogue["tstart"], expected_tstart, rtol=1e-8)
    assert catalogue["towards"].tolist() == expected_towards


def test_radec_to_unitvec_axis_points():
    assert np.allclose(radec_to_unitvec(0.0, 0.0), [1.0, 0.0, 0.0], atol=1e-12)
    assert np.allclose(radec_to_unitvec(90.0, 0.0), [0.0, 1.0, 0.0], atol=1e-12)
    assert np.allclose(radec_to_unitvec(0.0, 90.0), [0.0, 0.0, 1.0], atol=1e-12)


def test_local_basis_at_radec_orthonormal():
    for ra_deg, dec_deg in [(0.0, 0.0), (83.0, -30.0), (200.0, 60.0), (10.0, -89.0)]:
        e_east, e_north, e_los = local_basis_at_radec(ra_deg, dec_deg)

        for v in (e_east, e_north, e_los):
            assert np.isclose(np.linalg.norm(v), 1.0)

        assert np.isclose(np.dot(e_east, e_north), 0.0, atol=1e-12)
        assert np.isclose(np.dot(e_east, e_los), 0.0, atol=1e-12)
        assert np.isclose(np.dot(e_north, e_los), 0.0, atol=1e-12)

        # e_los is the same line-of-sight unit vector radec_to_unitvec computes.
        assert np.allclose(e_los, radec_to_unitvec(ra_deg, dec_deg))


def test_angular_distance_known_cases():
    # angular_distance's ra/dec arguments are in radians, unlike
    # radec_to_unitvec/local_basis_at_radec above (degrees).
    assert np.isclose(angular_distance(0.0, 0.0, 0.0, 0.0), 0.0, atol=1e-12)
    assert np.isclose(angular_distance(0.0, 0.0, np.pi, 0.0), np.pi)  # antipodal
    assert np.isclose(angular_distance(0.0, 0.0, np.pi / 2, 0.0), np.pi / 2)  # quarter turn, equator
    assert np.isclose(angular_distance(0.0, np.pi / 2, 1.234, 0.0), np.pi / 2)  # pole to equator


def test_pointing_to_radec_circle_radius_invariant():
    pointing_ra_deg, pointing_dec_deg = 83.0, -30.0
    Lfov_deg = 10.0
    Npoints = 200

    circle_vecs, circle_ra, circle_dec = pointing_to_radec_circle(
        pointing_ra_deg, pointing_dec_deg, Lfov_deg, Npoints=Npoints
    )

    assert np.allclose(np.linalg.norm(circle_vecs, axis=1), 1.0)

    # Every returned boundary point sits at exactly Lfov/2 from the pointing center.
    separations = angular_distance(
        np.deg2rad(pointing_ra_deg), np.deg2rad(pointing_dec_deg), circle_ra, circle_dec
    )
    assert np.allclose(separations, np.deg2rad(Lfov_deg / 2.0))


def test_sampled_satellites_enter_the_fov():
    # Independent check that every sampled satellite's orbit (i, h,
    # long_asc_node, periapsis) genuinely brings it within Lfov/2 of the
    # pointing center at some point during the exposure -- not just that the
    # catalogue's own rastart/decstart/tstart bookkeeping is self-consistent
    # (that's covered by the snapshot test above). A 3-shell config is used
    # (rather than a single shell) so this is exercised across more than a
    # couple of satellites and more than one inclination.
    obs = _make_dummy_obs()
    shells_df = pd.DataFrame({
        "n": [1000.0, 800.0, 600.0],
        "h": [550.0, 650.0, 450.0],
        "i": [53.0, 70.0, 97.6],
    })
    t_init_mjd = Time("2026-01-01T00:00:00")
    Lfov_deg = 10.0

    integral_obs_model = IntegralObsModel(
        obs=obs,
        shells_df=shells_df,
        Lfov=Lfov_deg * u.deg,
        t_exp_mjd=0.01,
        t_init_mjd=t_init_mjd,
        ndec=20,
        nra=40,
        Npoints=50,
        dt=5.0 * u.s,  # fine enough to not clump every crossing at the last timestep
    )
    catalogue = integral_obs_model.sample_satellites(seed=97)
    assert len(catalogue) > 0

    t_exp_s = integral_obs_model.t_exp_mjd * 86400
    # ~1 sample/second: measured dwell times for this config are ~20-35s, so
    # this leaves a comfortable margin against missing a brief crossing
    # between grid points.
    n_t = int(t_exp_s) + 1
    t_array_s = np.linspace(0, t_exp_s, n_t)

    positions_gcrs = satellite_positions_gcrs(catalogue, t_array_s)

    pointing_ra_deg = integral_obs_model.initial_multi_shell_fov.pointing_ra.to_value(u.deg)
    pointing_dec_deg = integral_obs_model.initial_multi_shell_fov.pointing_dec.to_value(u.deg)
    e_east, e_north, e_center = local_basis_at_radec(pointing_ra_deg, pointing_dec_deg)

    # obsloc.get_gcrs_posvel is evaluated at every sample in obs_times, so the
    # observer's own motion from Earth's rotation is accounted for here, not
    # just the satellites' motion.
    obs_times = integral_obs_model.t_init_mjd + t_array_s * u.s
    _, _, n = positions_to_lm(positions_gcrs, integral_obs_model.obsloc, obs_times, e_east, e_north, e_center)

    # True angular separation from the pointing center (not the flat-sky l,m
    # approximation), so a satellite on the far/antipodal hemisphere (n<0)
    # correctly reads as far away rather than falsely close.
    separation_deg = np.rad2deg(np.arccos(np.clip(n, -1.0, 1.0)))
    min_separation_deg = separation_deg.min(axis=1)

    # A satellite sampled right as it crosses the boundary can sit exactly at
    # Lfov/2 (not strictly inside it), hence <= with a small tolerance.
    assert np.all(min_separation_deg <= Lfov_deg / 2.0 + 1e-6)
