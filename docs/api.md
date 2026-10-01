# API Reference

This page documents the public functions exposed by Analytical-Satsky.

```python
from analytical_satsky import *
```

## Overview

The public API is organized into three categories:

### Constellation management

* [`list_constellations`](#analytical_satsky.list_constellations)
* [`load_constellation`](#analytical_satsky.load_constellation)

### Satellite density and occupancy modelling

* [`SingleShellObs`](#analytical_satsky.SingleShellObs)
* [`MultiShellObs`](#analytical_satsky.MultiShellObs)
* [`compute_occupancy_fraction`](#analytical_satsky.compute_occupancy_fraction)

### Large FoV model and satellite sampling

* [`SingleShellFoV`](#analytical_satsky.SingleShellFoV)
* [`MultiShellFoV`](#analytical_satsky.MultiShellFoV)
* [`SingleShellFlux`](#analytical_satsky.SingleShellFlux)
* [`MultiShellFlux`](#analytical_satsky.MultiShellFlux)
* [`IntegralObsModel`](#analytical_satsky.IntegralObsModel)

### Visualisation

* [`plot_sky_map`](#analytical_satsky.plot_sky_map)

---

# Constellation management

::: analytical_satsky.list_constellations

::: analytical_satsky.load_constellation

For more on constellations see the [constellation page](constellations.md).

---

# Satellite density and occupancy modelling

::: analytical_satsky.SingleShellObs

::: analytical_satsky.MultiShellObs

::: analytical_satsky.compute_occupancy_fraction

For information on the model assumptions see the [model page](model.md).

---

# Large FoV model and satellite sampling

This path parametrizes targets by true RA/Dec and an explicit observation time, instead of the
local-hour-angle convention used above.

::: analytical_satsky.SingleShellFoV

::: analytical_satsky.MultiShellFoV

::: analytical_satsky.SingleShellFlux

::: analytical_satsky.MultiShellFlux

::: analytical_satsky.IntegralObsModel

---

# Visualisation

::: analytical_satsky.plot_sky_map
