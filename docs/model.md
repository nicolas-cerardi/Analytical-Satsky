# Model

## Main assumpions

The current version of the model take the following assumptions:

* Satellite orbits are assumed to be purely circular, which results in perfectly spherical shells.
* Satellites shells are defined w.r.t the earth instantaneaous rotation axis. This allows to drop the dependency on the time of the observation, and to substitute the target right ascension by its local hour angle.
* Constant earth radius: we assume a latitude independant earth radius, which value is taken from Astropy. This shapes the equation implemented in `compute_d_phi`.
* Spatial and time integration (see function `compute_nsats`): in the instrument FoV we assume straight trajectories, constant velocity and constant satellite density per shell.
