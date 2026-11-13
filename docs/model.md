# Analytical Model

This package implements two analytical models of satellite orbital shells:

- The small FoV model, similar to [Bassa et al. 2022](https://www.aanda.org/articles/aa/abs/2022/01/aa42101-21/aa42101-21.html), originally developed for optical astronomy.
- The large FoV model, dropping important assumptions that do not hold for radio telescopes that have a substantially larger FoV.

Both models take the following assumption:

- Satellite orbits are assumed to be purely circular, which results in perfectly spherical shells.

## Small FoV Model

For a complete description of this, see [Bassa et al. 2022](https://www.aanda.org/articles/aa/abs/2022/01/aa42101-21/aa42101-21.html). Note that this version of the model takes the following assumptions:

- Satellite shells are defined w.r.t. the Earth's instantaneous rotation axis. This allows to drop the dependency on the time of the observation, and to substitute the target right ascension by its local hour angle.
- Constant Earth radius: we assume a latitude-independent Earth radius, which value is taken from Astropy. This shapes the equation implemented in `compute_d_phi`.
- Spatial and time integration (see function `compute_nsats`): in the instrument FoV we assume straight trajectories, constant velocity and constant satellite density per shell.

## Large FoV Model

- Observer-satellite distance $d$: the model solves for the equation

    $$d^2 + 2 (z_{\rm OB} \sin \delta + y_{\rm OB} \cos \delta \cos \alpha + z_{\rm OB} \cos \delta \cos \alpha) d + ||OB|| - (R_{\rm earth} + h_{\rm sat})^2$$

- The number of satellites entering the FoV during the observation is separated in $N_{\rm obs} = N_{\rm init} + N_{\rm flux}$, with $N_{\rm init}$ the initial instantaneous number of satellites, computed as the integral of the density over the FoV:

    $$N_{\rm init} = \int_{\rm FoV} \rho_{\rm sat} \, {\rm d}\Omega ,$$

    and $N_{\rm flux}$ the flux of ingressing satellites:

    $$N_{\rm flux} = \oint \rho_{\rm sat} \, \vec{v}_{\rm sat} \times R_{\rm FoV} \, \vec{u}_r \, {\rm d}\theta .$$
