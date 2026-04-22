def test_package_imports():
    import analytical_satsky
    from analytical_satsky import (
        list_constellations,
        load_constellation,
        compute_total_satellite_density,
        compute_shell_satellite_density,
        simulate_exposed_time,
        compute_exposure_fraction,
        ra_to_lha,
        plot_sky_map,
    )
    