import importlib

import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_allclose

import skyborn.calc.srh as srh_package
import skyborn.calc.srh.core as srh_core


def _sample_profile():
    """metpy 1.7.1 docstring sounding for storm_relative_helicity."""
    height = np.array([250.0, 700.0, 1500.0, 3100.0, 5720.0, 7120.0])
    wdir_deg = np.array([165.0, 180.0, 190.0, 210.0, 220.0, 250.0])
    speed_knots = np.array([5.0, 15.0, 20.0, 30.0, 50.0, 60.0])
    rad = np.deg2rad(wdir_deg)
    speed_ms = speed_knots * 0.5144444444444445
    # meteorological "from" convention, as metpy's wind_components
    u = -speed_ms * np.sin(rad)
    v = -speed_ms * np.cos(rad)
    return height, u, v


def _sample_grid():
    height, u, v = _sample_profile()
    return (
        np.broadcast_to(height[:, None, None], (6, 2, 3)).copy(),
        np.broadcast_to(u[:, None, None], (6, 2, 3)).copy(),
        np.broadcast_to(v[:, None, None], (6, 2, 3)).copy(),
    )


def test_calc_package_exports_srh_submodule_and_function():
    calc_module = importlib.reload(importlib.import_module("skyborn.calc"))
    assert calc_module.srh is srh_package
    assert (
        calc_module.calculate_storm_relative_helicity
        is srh_package.calculate_storm_relative_helicity
    )


def test_srh_matches_metpy_reference():
    """metpy docstring case: depth=1km, storm=(7, 7) m/s -> 49.6086162."""
    height, u, v = _sample_profile()
    result = srh_package.calculate_storm_relative_helicity(
        height, u, v, depth=1000.0, storm_u=7.0, storm_v=7.0
    )
    assert_allclose(result[0], 49.6086162, rtol=1e-9)
    assert_allclose(result[1], 0.0, atol=1e-12)
    assert_allclose(result[2], 49.6086162, rtol=1e-9)


def test_srh_default_storm_is_zero():
    height, u, v = _sample_profile()
    result = srh_package.calculate_storm_relative_helicity(height, u, v, depth=1000.0)
    # metpy with storm=0
    assert_allclose(result[0], 14.6158, rtol=1e-3)
    assert_allclose(result[2], 14.6158, rtol=1e-3)


def test_srh_public_api_matches_backend():
    height, u, v = _sample_profile()
    public = srh_package.calculate_storm_relative_helicity(
        height, u, v, depth=1000.0, storm_u=7.0, storm_v=7.0
    )
    backend = srh_core.srh_profile(height, u, v, 1000.0, 0.0, 7.0, 7.0)
    assert_allclose(public, backend, rtol=1e-12, atol=1e-12)


def test_srh_grid_matches_backend_and_xarray():
    h3, u3, v3 = _sample_grid()

    grid = srh_core.srh_grid(h3, u3, v3, 1000.0, 0.0, 7.0, 7.0)
    profile = srh_core.srh_profile(h3[:, 0, 0], u3[:, 0, 0], v3[:, 0, 0],
                                   1000.0, 0.0, 7.0, 7.0)
    for i, arr in enumerate(grid):
        assert_allclose(arr[0, 0], profile[i], rtol=1e-12, atol=1e-12)

    h_xr = xr.DataArray(
        h3,
        dims=["level", "lat", "lon"],
        coords={"level": h3[:, 0, 0], "lat": [0.0, 10.0], "lon": [100.0, 110.0, 120.0]},
    )
    u_xr = xr.DataArray(u3, dims=["level", "lat", "lon"])
    v_xr = xr.DataArray(v3, dims=["level", "lat", "lon"])

    public = srh_package.calculate_storm_relative_helicity(
        h_xr, u_xr, v_xr, depth=1000.0, storm_u=7.0, storm_v=7.0
    )
    assert all(isinstance(r, xr.DataArray) for r in public)
    assert public[0].dims == ("lat", "lon")
    for i, arr in enumerate(grid):
        assert_allclose(public[i].values, arr, rtol=1e-12, atol=1e-12)


def test_srh_bottom_offset():
    """A non-zero bottom changes the integration layer."""
    height, u, v = _sample_profile()
    surface = srh_package.calculate_storm_relative_helicity(
        height, u, v, depth=1000.0, storm_u=7.0, storm_v=7.0
    )
    elevated = srh_package.calculate_storm_relative_helicity(
        height, u, v, depth=1000.0, bottom=200.0, storm_u=7.0, storm_v=7.0
    )
    assert elevated[0] != surface[0]


def test_srh_nan_handling():
    height, u, v = _sample_profile()
    h2 = height.copy()
    u2 = u.copy()
    v2 = v.copy()
    h2[2] = np.nan
    u2[2] = np.nan
    v2[2] = np.nan
    result = srh_package.calculate_storm_relative_helicity(
        h2, u2, v2, depth=1000.0, storm_u=7.0, storm_v=7.0
    )
    assert all(np.isfinite(r) for r in result)


def test_srh_input_validation():
    with pytest.raises(ValueError, match="height must be 1D or 3D"):
        srh_package.calculate_storm_relative_helicity(
            np.ones((2, 2), dtype=float),
            np.ones((2, 2), dtype=float),
            np.ones((2, 2), dtype=float),
            depth=1000.0,
        )

    with pytest.raises(ValueError, match="height must be 1D"):
        srh_core.srh_profile(
            np.ones((2, 2), dtype=float),
            np.ones((2, 2), dtype=float),
            np.ones((2, 2), dtype=float),
            1000.0,
        )
