import importlib

import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_allclose

import skyborn.calc.dcape as dcape_package
import skyborn.calc.dcape.core as dcape_core


def _sample_profile():
    pressure = np.array([1000.0, 850.0, 700.0, 600.0, 500.0], dtype=float)
    temperature = np.array([30.0, 22.0, 12.0, 5.0, -3.0], dtype=float)
    dewpoint = np.array([22.0, 14.0, 3.0, -4.0, -12.0], dtype=float)
    return pressure, temperature, dewpoint


def test_calc_package_exports_dcape_submodule_and_function():
    calc_module = importlib.reload(importlib.import_module("skyborn.calc"))

    assert calc_module.dcape is dcape_package
    assert calc_module.calculate_dcape is dcape_package.calculate_dcape


def test_dcape_public_api_matches_backend_profile():
    pressure, temperature, dewpoint = _sample_profile()

    public_result = dcape_package.calculate_dcape(pressure, temperature, dewpoint)
    backend_result = dcape_core.dcape_profile(pressure, temperature, dewpoint)

    assert_allclose(public_result, backend_result, rtol=1e-12, atol=1e-12)


def test_dcape_public_api_matches_backend_grid_and_xarray():
    pressure, temperature, dewpoint = _sample_profile()
    pressure_3d = np.broadcast_to(pressure[:, None, None], (5, 2, 3)).copy()
    temperature_3d = np.broadcast_to(temperature[:, None, None], (5, 2, 3)).copy()
    dewpoint_3d = np.broadcast_to(dewpoint[:, None, None], (5, 2, 3)).copy()

    grid_result = dcape_package.calculate_dcape(
        pressure_3d, temperature_3d, dewpoint_3d
    )
    backend_grid = dcape_core.dcape_grid(pressure_3d, temperature_3d, dewpoint_3d)

    assert_allclose(grid_result, backend_grid, rtol=1e-12, atol=1e-12)

    pressure_xr = xr.DataArray(
        pressure_3d,
        dims=["level", "lat", "lon"],
        coords={"level": pressure, "lat": [0.0, 10.0], "lon": [100.0, 110.0, 120.0]},
        attrs={"units": "hPa"},
    )
    temperature_xr = xr.DataArray(temperature_3d, dims=["level", "lat", "lon"])
    dewpoint_xr = xr.DataArray(dewpoint_3d, dims=["level", "lat", "lon"])

    xarray_result = dcape_package.calculate_dcape(
        pressure_xr, temperature_xr, dewpoint_xr
    )

    assert isinstance(xarray_result, xr.DataArray)
    assert xarray_result.dims == ("lat", "lon")
    assert_allclose(xarray_result.values, backend_grid, rtol=1e-12, atol=1e-12)


def test_dcape_backend_is_consistent_between_profile_and_grid():
    pressure, temperature, dewpoint = _sample_profile()
    pressure_3d = np.broadcast_to(pressure[:, None, None], (5, 2, 3)).copy()
    temperature_3d = np.broadcast_to(temperature[:, None, None], (5, 2, 3)).copy()
    dewpoint_3d = np.broadcast_to(dewpoint[:, None, None], (5, 2, 3)).copy()

    profile_result = dcape_core.dcape_profile(pressure, temperature, dewpoint)
    grid_result = dcape_core.dcape_grid(pressure_3d, temperature_3d, dewpoint_3d)

    assert_allclose(
        grid_result,
        np.full((2, 3), profile_result, dtype=np.float64),
        rtol=1e-12,
        atol=1e-12,
    )


def test_dcape_input_validation():
    with pytest.raises(ValueError, match="pressure must be 1D"):
        dcape_package.calculate_dcape(
            np.ones((2, 2), dtype=float),
            np.ones((2, 2), dtype=float),
            np.ones((2, 2), dtype=float),
        )
