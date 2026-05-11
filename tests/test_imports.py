def test_package_imports():
    import analytical_satsky
    from analytical_satsky import (
        list_constellations,
        load_constellation,
        compute_occupancy_fraction,
        plot_sky_map,
        MultiShellObs,
        SingleShellObs
    )
    