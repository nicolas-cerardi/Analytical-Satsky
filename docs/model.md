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

For a shell of $N_{\rm sat}$ satellites at altitude $h_{\rm sat}$ and inclination $i$, seen by an observer at latitude $\phi_{\rm obs}$ pointing at declination $\delta$ and local hour angle $H$, the model computes:

- The observer-satellite distance $d$ along the line of sight, as the positive root of

    $$d^2 + 2 R_{\rm earth} \left(\cos\delta \cos H \cos\phi_{\rm obs} + \sin\delta \sin\phi_{\rm obs}\right) d - \left(h_{\rm sat}^2 + 2 R_{\rm earth} h_{\rm sat}\right) = 0 ,$$

    which also gives the latitude $\phi_{\rm sat}$ of the point where the line of sight crosses the shell.

- The probability density for a single satellite to lie at latitude $\phi_{\rm sat}$, per unit area of the shell:

    $$P(\phi_{\rm sat}) = \frac{1}{2 \pi^2 (R_{\rm earth} + h_{\rm sat})^2 \sqrt{\sin^2 i - \sin^2 \phi_{\rm sat}}} \quad {\rm for}\ |\phi_{\rm sat}| < i ,$$

    and $P = 0$ otherwise.

- The satellite density on the sky, per unit solid angle:

    $$\rho_{\rm sat} = N_{\rm sat} \, P(\phi_{\rm sat}) \, \frac{d^2}{\cos\alpha} , \qquad \cos\alpha = \frac{(R_{\rm earth} + h_{\rm sat})^2 + d^2 - R_{\rm earth}^2}{2 \, d \, (R_{\rm earth} + h_{\rm sat})} ,$$

    with $\alpha$ the impact angle between the line of sight and the normal to the shell.

- The number of satellites seen in a circular FoV of diameter $L_{\rm FoV}$ during an observation of duration $T_{\rm obs}$:

    $$N_{\rm obs} = \rho_{\rm sat} \left( \pi \left(\frac{L_{\rm FoV}}{2}\right)^2 + L_{\rm FoV} \, w_{\rm sat} \, T_{\rm obs} \right) ,$$

    where $w_{\rm sat}$ is the apparent angular velocity of the satellites, i.e. their velocity relative to the observer, projected perpendicular to the line of sight and divided by $d$, averaged over the north- and south-moving satellites. The first term counts the satellites already in the FoV at the start of the observation, the second those entering it afterwards.

The total over a constellation is the sum of $N_{\rm obs}$ over its shells.

## Large FoV Model

This repository also implements an extension to the previous model, dropping two major assumptions. First, the density $\rho_{sat}$ is generally not constant within the FoV, and has strong singularities at shell edges. Following the same logic, satellite ingress probability is not constant along the FoV perimeter. The larger the FoV, the higher the error of the constant density approximation from the original model. We hence rewrite the above equations with integrals. In addition, we drop the fixed earth radius assumption, and adopt the ellipsoidal model, and we drop the straight trajectory assumption.

- Observer-satellite distance $d$: the model solves for the equation

    $$d^2 + 2\,(x_{\rm OB} \cos\delta \cos\alpha + y_{\rm OB} \cos\delta \sin\alpha + z_{\rm OB} \sin\delta)\, d + ||OB||^2 - (R_{\rm earth} + h_{\rm sat})^2 = 0 ,$$

    where $OB = (x_{\rm OB}, y_{\rm OB}, z_{\rm OB})$ is the vector from the Earth's centre to the observer, in the GCRS frame at the time of the observation, and $(\alpha, \delta)$ are the right ascension and declination of the line of sight. The positive root is kept.

- The number of satellites entering the FoV during the observation is separated in $N_{\rm obs} = N_{\rm init} + N_{\rm flux}$, with $N_{\rm init}$ the initial instantaneous number of satellites, computed as the integral of the density over the FoV:

    $$N_{\rm init} = \int_{\rm FoV} \rho_{\rm sat} \, {\rm d}\Omega ,$$

    and $N_{\rm flux}$ the flux of ingressing satellites:

    $$N_{\rm flux} = \int_0^{T_{\rm obs}} \oint \sum_{k \in \{N,S\}} \frac{\rho_{\rm sat}}{2}\, \max\left(0, -\vec{w}_k \cdot \hat{n}\right) R_{\rm FoV}\, {\rm d}\theta\, {\rm d}t$$

## Coordinate conventions

The two models do not describe the target in the same way:

| | Small FoV model | Large FoV model |
|---|---|---|
| Target position | Declination and local hour angle | True right ascension and declination |
| Time | None: the result does not depend on the time of the observation | Explicit observation time (MJD) |
| Frame | Rotating with the Earth, observer on the reference meridian | GCRS |
| Classes | `SingleShellObs`, `MultiShellObs` | `SingleShellFoV`, `MultiShellFoV`, `SingleShellFlux`, `MultiShellFlux`, `IntegralObsModel` |

A target given in right ascension and declination must therefore be converted to a local hour angle at the time of the observation before being passed to the small FoV model.
