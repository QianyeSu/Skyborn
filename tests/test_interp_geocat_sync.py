"""Focused regression tests for the GeoCAT interpolation sync."""

import numpy as np
import xarray as xr
from numpy.testing import assert_allclose, assert_array_equal

from skyborn.interp.interpolation import (
    _pressure_from_hybrid,
    _temp_extrapolate,
    delta_pressure_hybrid,
    interp_hybrid_to_pressure,
    pressure_at_hybrid_levels,
)


def test_pressure_helpers_broadcast_with_default_xarray_dims():
    """Public and legacy pressure helpers should agree under dim-name collisions."""

    ps = xr.DataArray([100000.0, 95000.0])
    hya = xr.DataArray([0.0, 10000.0, 50000.0])
    hyb = xr.DataArray([1.0, 0.8, 0.0])

    pressure_public = pressure_at_hybrid_levels(ps, hya, hyb)
    pressure_private = _pressure_from_hybrid(ps, hya, hyb)
    dph = delta_pressure_hybrid(ps, hya, hyb)

    assert pressure_public.shape == (3, 2)
    assert_array_equal(pressure_public.values, pressure_private.values)
    assert dph.shape == (2, 2)
    assert np.all(dph.values >= 0)


def test_temp_extrapolate_supports_legacy_and_explicit_t_bot_calls():
    """Temperature extrapolation should accept both Skyborn and GeoCAT call styles."""

    full = xr.DataArray(
        np.array([[[280.0]], [[260.0]]]),
        dims=["lev", "lat", "lon"],
        coords={"lev": [85000.0, 70000.0], "lat": [0.0], "lon": [0.0]},
    )
    t_bot = full.isel(lev=-1, drop=True)
    p_sfc = xr.DataArray([[85000.0]], dims=["lat", "lon"], coords=t_bot.coords)
    ps = xr.DataArray([[100000.0]], dims=["lat", "lon"], coords=t_bot.coords)
    phi_sfc = xr.DataArray([[0.0]], dims=["lat", "lon"], coords=t_bot.coords)

    explicit = _temp_extrapolate(t_bot, 100000.0, p_sfc, ps, phi_sfc)
    legacy = _temp_extrapolate(full, "lev", 100000.0, p_sfc, ps, phi_sfc)

    assert_allclose(explicit.values, legacy.values)


def test_interp_hybrid_to_pressure_accepts_different_coeff_dimension_name():
    """Hybrid coefficients should be renamed onto the data lev dim when needed."""

    data = xr.DataArray(
        np.arange(2 * 4 * 3, dtype=float).reshape(2, 4, 3),
        dims=["time", "lev", "x"],
        coords={"time": [0, 1], "lev": [0, 1, 2, 3], "x": [0, 1, 2]},
        attrs={"units": "K", "long_name": "temperature"},
    )
    ps = xr.DataArray(
        np.full((2, 3), 100000.0),
        dims=["time", "x"],
        coords={"time": data.time, "x": data.x},
    )
    hyam = xr.DataArray([0.0, 0.1, 0.3, 0.5], dims=["hlev"])
    hybm = xr.DataArray([1.0, 0.8, 0.4, 0.0], dims=["hlev"])

    out = interp_hybrid_to_pressure(
        data=data,
        ps=ps,
        hyam=hyam,
        hybm=hybm,
        new_levels=np.array([95000.0, 70000.0]),
        lev_dim="lev",
    )

    assert out.dims == ("time", "plev", "x")
    assert_array_equal(out.plev.values, np.array([95000.0, 70000.0]))
    assert out.attrs["units"] == "K"
    assert out.attrs["long_name"] == "temperature"


def test_interp_hybrid_to_pressure_temperature_extrapolation_uses_explicit_t_bot():
    """Temperature extrapolation should work with the explicit bottom-temperature path."""

    data = xr.DataArray(
        np.array([[[300.0]], [[280.0]], [[260.0]]]),
        dims=["lev", "lat", "lon"],
        coords={"lev": [0, 1, 2], "lat": [0.0], "lon": [0.0]},
    )
    ps = xr.DataArray(
        [[101325.0]],
        dims=["lat", "lon"],
        coords={"lat": [0.0], "lon": [0.0]},
    )
    hyam = xr.DataArray([0.0, 0.3, 0.6], dims=["lev"])
    hybm = xr.DataArray([1.0, 0.5, 0.0], dims=["lev"])
    t_bot = data.isel(lev=-1, drop=True)
    phi_sfc = xr.DataArray(
        [[0.0]], dims=["lat", "lon"], coords={"lat": [0.0], "lon": [0.0]}
    )

    out = interp_hybrid_to_pressure(
        data=data,
        ps=ps,
        hyam=hyam,
        hybm=hybm,
        new_levels=np.array([110000.0, 70000.0]),
        lev_dim="lev",
        extrapolate=True,
        variable="temperature",
        t_bot=t_bot,
        phi_sfc=phi_sfc,
    )

    assert out.shape == (2, 1, 1)
    assert np.all(np.isfinite(out.values))
