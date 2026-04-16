"""
Enhanced tests for skyborn.interp.interpolation module.

This module provides additional test coverage for edge cases, performance scenarios,
and specific code paths that may not be covered by the existing test suite.
Focuses on improving code coverage for the interpolation functionality.
"""

import warnings

import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_allclose, assert_array_almost_equal, assert_array_equal

import skyborn.interp.interpolation as interpolation_mod
from skyborn.interp.interpolation import (
    __pres_lev_mandatory__,
    _align_hybrid_level_dimension,
    _func_interpolate,
    _interpolate_mb,
    _is_pint_backed,
    _pressure_from_hybrid,
    _rename_colliding_coeff_dim,
    _sigma_from_hybrid,
    _strip_unexpected_pint_units,
    delta_pressure_hybrid,
    interp_hybrid_to_pressure,
    interp_multidim,
    interp_sigma_to_hybrid,
    pressure_at_hybrid_levels,
)

# Try to import private functions - they may not be available
try:
    from skyborn.interp.interpolation import (
        _geo_height_extrapolate,
        _post_interp_multidim,
        _pre_interp_multidim,
        _temp_extrapolate,
        _vertical_remap,
        _vertical_remap_extrap,
    )

    PRIVATE_FUNCTIONS_AVAILABLE = True
except ImportError:
    PRIVATE_FUNCTIONS_AVAILABLE = False


class TestInterpolationEdgeCases:
    """Test edge cases and boundary conditions for interpolation functions."""

    def test_pressure_from_hybrid_scalar_inputs(self):
        """Test pressure calculation with scalar inputs."""
        ps = xr.DataArray(100000.0)  # Scalar surface pressure
        hya = xr.DataArray([0.0, 0.5])  # 2 levels
        hyb = xr.DataArray([1.0, 0.0])
        p0 = 100000.0

        pressure = _pressure_from_hybrid(ps, hya, hyb, p0)

        assert pressure.shape == (2,)
        assert pressure[0] == 100000.0  # Surface level
        assert pressure[1] == 50000.0  # Top level

    def test_pressure_from_hybrid_zero_surface_pressure(self):
        """Test pressure calculation with zero surface pressure."""
        ps = xr.DataArray([0.0, 50000.0], dims=["x"])
        hya = xr.DataArray([0.0, 25000.0], dims=["lev"])
        hyb = xr.DataArray([1.0, 0.5], dims=["lev"])
        p0 = 100000.0

        pressure = _pressure_from_hybrid(ps, hya, hyb, p0)

        # First point should have zero pressure at surface
        assert pressure.shape == (2, 2)
        assert pressure[0, 0] == 0.0
        assert pressure[1, 0] > 0.0  # Second point should be positive

    def test_sigma_from_hybrid_edge_values(self):
        """Test sigma calculation with edge values."""
        ps = xr.DataArray([101325.0], dims=["x"])
        # Test with hya=0 (pure sigma) and hyb=0 (pure pressure)
        hya = xr.DataArray([0.0, 50000.0 / 101325.0, 1.0], dims=["lev"])
        hyb = xr.DataArray([1.0, 0.0, 0.0], dims=["lev"])
        p0 = 101325.0

        sigma = _sigma_from_hybrid(ps, hya, hyb, p0)

        assert sigma.shape == (3, 1)
        assert_array_almost_equal(sigma[0, 0], 1.0, decimal=10)  # Pure sigma level
        assert_array_almost_equal(
            sigma[1, 0], 50000.0 / 101325.0, decimal=6
        )  # Pure pressure
        # At reference pressure
        assert_array_almost_equal(sigma[2, 0], 1.0, decimal=10)

    def test_func_interpolate_function_properties(self):
        """Test properties of interpolation functions."""
        # Test linear interpolation
        func_linear = _func_interpolate("linear")

        # Create test data
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([10.0, 20.0, 30.0])
        xi = np.array([1.5, 2.5])

        result = func_linear(xi, x, y)
        expected = np.array([15.0, 25.0])  # Linear interpolation
        assert_array_almost_equal(result, expected)

        # Test log interpolation
        func_log = _func_interpolate("log")

        # Positive values for log interpolation
        x_log = np.array([1.0, 10.0, 100.0])
        y_log = np.array([1.0, 10.0, 100.0])
        xi_log = np.array([3.16227766])  # sqrt(10)

        result_log = func_log(xi_log, x_log, y_log)
        assert len(result_log) == 1
        assert result_log[0] > 1.0 and result_log[0] < 100.0

    def test_hybrid_to_pressure_exact_surface_match(self):
        """Test hybrid to pressure interpolation for an exact surface-pressure match."""
        # Minimal supported two-level data
        data = xr.DataArray(
            [
                [[280.0, 285.0], [275.0, 290.0]],
                [[240.0, 245.0], [235.0, 250.0]],
            ],
            dims=["lev", "lat", "lon"],
            coords={"lat": [0, 30], "lon": [0, 90]},
        )

        ps = xr.DataArray(
            [[101325.0, 101325.0], [100000.0, 100000.0]],
            dims=["lat", "lon"],
            coords={"lat": [0, 30], "lon": [0, 90]},
        )

        hya = xr.DataArray([0.0, 0.2], dims=["lev"])
        hyb = xr.DataArray([1.0, 0.0], dims=["lev"])

        new_levels = np.array([101325.0, 95000.0])

        result = interp_hybrid_to_pressure(
            data=data,
            ps=ps,
            hyam=hya,
            hybm=hyb,
            new_levels=new_levels,
            lev_dim="lev",
        )

        assert result.shape == (2, 2, 2)  # plev, lat, lon
        # Only the grid row whose surface pressure matches the target level is exact.
        assert_array_almost_equal(result[0, 0, :], [280.0, 285.0])
        assert np.all(np.isnan(result[0, 1, :]))

    def test_hybrid_to_pressure_temperature_extrapolation(self):
        """Test temperature extrapolation in hybrid to pressure interpolation."""
        # Create test data with temperature
        nlev = 5
        data = xr.DataArray(
            np.array([[[300.0], [280.0], [260.0], [240.0], [220.0]]]).T,  # 1x5x1
            dims=["lat", "lev", "lon"],
            coords={"lat": [0], "lev": range(nlev), "lon": [0]},
        )

        ps = xr.DataArray(
            [[101325.0]], dims=["lat", "lon"], coords={"lat": [0], "lon": [0]}
        )
        hya = xr.DataArray(np.linspace(0.0, 0.2, nlev), dims=["lev"])
        hyb = xr.DataArray(np.linspace(1.0, 0.0, nlev), dims=["lev"])
        t_bot = data.isel(lev=-1, drop=True)
        phi_sfc = xr.DataArray(
            [[0.0]], dims=["lat", "lon"], coords={"lat": [0], "lon": [0]}
        )

        # Request levels both within and outside the data range
        # surface, mid, top
        new_levels = np.array([110000.0, 70000.0, 20000.0])

        result = interp_hybrid_to_pressure(
            data=data,
            ps=ps,
            hyam=hya,
            hybm=hyb,
            new_levels=new_levels,
            lev_dim="lev",
            extrapolate=True,
            variable="temperature",
            t_bot=t_bot,
            phi_sfc=phi_sfc,
        )

        assert result.shape == (1, 3, 1)  # lat, plev, lon
        assert np.all(np.isfinite(result.values))

        # Temperature should increase towards surface (extrapolation)
        # Surface warmer than mid-level
        assert result[0, 0, 0] >= result[0, 1, 0]

    def test_hybrid_to_pressure_geopotential_extrapolation(self):
        """Test geopotential height extrapolation."""
        # Create geopotential height data (increases with altitude)
        nlev = 4
        data = xr.DataArray(
            np.array([[[1000.0], [5000.0], [10000.0], [15000.0]]]).T,  # m
            dims=["lat", "lev", "lon"],
            coords={"lat": [0], "lev": range(nlev), "lon": [0]},
        )

        ps = xr.DataArray(
            [[100000.0]], dims=["lat", "lon"], coords={"lat": [0], "lon": [0]}
        )
        hya = xr.DataArray(np.linspace(0.0, 0.1, nlev), dims=["lev"])
        hyb = xr.DataArray(np.linspace(1.0, 0.0, nlev), dims=["lev"])
        t_bot = xr.DataArray(
            [[280.0]], dims=["lat", "lon"], coords={"lat": [0], "lon": [0]}
        )
        phi_sfc = xr.DataArray(
            [[0.0]], dims=["lat", "lon"], coords={"lat": [0], "lon": [0]}
        )

        new_levels = np.array([105000.0, 50000.0, 10000.0])

        result = interp_hybrid_to_pressure(
            data=data,
            ps=ps,
            hyam=hya,
            hybm=hyb,
            new_levels=new_levels,
            lev_dim="lev",
            extrapolate=True,
            variable="geopotential",
            t_bot=t_bot,
            phi_sfc=phi_sfc,
        )

        assert result.shape == (1, 3, 1)
        # Geopotential should increase with altitude (decrease with pressure)
        assert result[0, 0, 0] <= result[0, 1, 0] <= result[0, 2, 0]

    def test_fortran_path_prefers_preallocated_output_buffer(self, monkeypatch):
        """Test that the Fortran fast path reuses a caller-allocated output buffer."""
        data = xr.DataArray(
            [
                [[280.0, 285.0]],
                [[240.0, 245.0]],
            ],
            dims=["lev", "lat", "lon"],
            coords={"lat": [0], "lon": [0, 90]},
        )
        ps = xr.DataArray(
            [[100000.0, 100000.0]],
            dims=["lat", "lon"],
            coords={"lat": [0], "lon": [0, 90]},
        )
        hyam = xr.DataArray([0.0, 0.2], dims=["lev"])
        hybm = xr.DataArray([1.0, 0.0], dims=["lev"])
        new_levels = np.array([90000.0, 80000.0], dtype=np.float64)

        captured = {}

        def fake_nodes_into(
            dati,
            dato,
            hbcofa,
            hbcofb,
            p0,
            plevo,
            intyp,
            psfc,
            spvl,
            kxtrp,
        ):
            captured["dato_shape"] = dato.shape
            captured["dato_f"] = dato.flags.f_contiguous
            captured["psfc_shape"] = psfc.shape
            dato[:, :] = 123.0

        def fail_return_allocating_nodes(*args, **kwargs):
            raise AssertionError("return-allocating vinth2p wrapper should not be used")

        monkeypatch.setattr(
            interpolation_mod, "_dvinth2p_nodes_pa_into", fake_nodes_into
        )
        monkeypatch.setattr(
            interpolation_mod,
            "_dvinth2p_nodes_pa",
            fail_return_allocating_nodes,
        )

        result = interpolation_mod._interp_hybrid_to_pressure_fortran(
            data=data,
            ps=ps,
            hyam=hyam,
            hybm=hybm,
            p0=100000.0,
            new_levels=new_levels,
            lev_dim="lev",
            method="linear",
            extrapolate=False,
            variable=None,
            t_bot=None,
            phi_sfc=None,
        )

        assert captured["dato_shape"] == (2, 2)
        assert captured["dato_f"] is True
        assert captured["psfc_shape"] == (2,)
        assert result.shape == (2, 1, 2)
        assert_array_equal(result.values, np.full((2, 1, 2), 123.0))

    @pytest.mark.skipif(
        not PRIVATE_FUNCTIONS_AVAILABLE, reason="Private functions not available"
    )
    def test_sigma_to_hybrid_edge_coordinates(self):
        """Test sigma to hybrid with edge coordinate values."""
        # Single time step, simple spatial grid
        data = xr.DataArray(
            [[[250.0, 260.0], [270.0, 280.0], [290.0, 300.0]]],  # 1x3x2
            dims=["time", "sigma", "spatial"],
            coords={"time": [0], "sigma": [0.2, 0.6, 1.0], "spatial": [0, 1]},
        )

        sig_coords = xr.DataArray([0.2, 0.6, 1.0], dims=["sigma"])
        ps = xr.DataArray([[101325.0, 100000.0]], dims=["time", "spatial"])

        # Target hybrid levels
        hya = xr.DataArray([0.0, 0.5], dims=["hlev"])
        hyb = xr.DataArray([1.0, 0.0], dims=["hlev"])

        result = interp_sigma_to_hybrid(
            data=data, sig_coords=sig_coords, ps=ps, hyam=hya, hybm=hyb, lev_dim="sigma"
        )

        assert result.shape == (1, 2, 2)  # time, hlev, spatial
        assert np.all(np.isfinite(result.values))

    def test_multidim_interpolation_single_point(self):
        """Test multidimensional interpolation to single output point."""
        # 3x3 input grid
        lat_in = np.array([0, 30, 60])
        lon_in = np.array([0, 90, 180])
        data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

        data_in = xr.DataArray(
            data, dims=["lat", "lon"], coords={"lat": lat_in, "lon": lon_in}
        )

        # Single output point
        lat_out = np.array([15])
        lon_out = np.array([45])

        result = interp_multidim(data_in=data_in, lat_out=lat_out, lon_out=lon_out)

        assert result.shape == (1, 1)
        assert np.isfinite(result.values[0, 0])

    def test_multidim_interpolation_boundary_points(self):
        """Test interpolation at grid boundary points."""
        lat_in = np.array([-90, 0, 90])
        lon_in = np.array([0, 180, 360])  # Note: 360 = 0
        data = np.array([[1, 2, 1], [3, 4, 3], [5, 6, 5]])

        data_in = xr.DataArray(
            data, dims=["lat", "lon"], coords={"lat": lat_in, "lon": lon_in}
        )

        # Request points exactly on grid boundaries
        lat_out = np.array([-90, 90])
        lon_out = np.array([0, 360])

        result = interp_multidim(data_in=data_in, lat_out=lat_out, lon_out=lon_out)

        assert result.shape == (2, 2)
        # Should match grid values at boundaries
        assert result[0, 0] == 1  # lat=-90, lon=0
        assert result[1, 0] == 5  # lat=90, lon=0


class TestInterpolationNumericalStability:
    """Test numerical stability and precision of interpolation functions."""

    def test_pressure_calculation_precision(self):
        """Test numerical precision in pressure calculations."""
        # Use high precision values
        ps = xr.DataArray([101325.0])
        hya = xr.DataArray([0.0, 0.001, 0.01])  # Small dimensionless increments
        hyb = xr.DataArray([1.0, 0.998, 0.95])
        p0 = 101325.0

        pressure = _pressure_from_hybrid(ps, hya, hyb, p0)

        # Check for reasonable precision
        assert pressure.shape == (3, 1)
        assert np.abs(pressure[0, 0] - 101325.0) < 1e-10
        assert pressure[1, 0] < pressure[0, 0]  # Should be decreasing
        assert pressure[2, 0] < pressure[1, 0]

    def test_public_hybrid_helpers_match_private_wrapper(self):
        """Test public hybrid helpers against the legacy private wrapper."""
        ps = xr.DataArray([100000.0, 95000.0], dims=["x"])
        hya = xr.DataArray([0.0, 10000.0, 50000.0], dims=["lev"])
        hyb = xr.DataArray([1.0, 0.8, 0.0], dims=["lev"])

        pressure_public = pressure_at_hybrid_levels(ps, hya, hyb)
        pressure_private = _pressure_from_hybrid(ps, hya, hyb)
        dph = delta_pressure_hybrid(ps, hya, hyb)

        assert_array_equal(pressure_public.values, pressure_private.values)
        assert dph.shape == (2, 2)
        assert np.all(dph.values >= 0)

    def test_interpolation_with_nan_values(self):
        """Test interpolation behavior with NaN values."""
        # Data with NaN values
        data = xr.DataArray(
            [[[250.0, np.nan], [np.nan, 280.0], [290.0, 300.0]]],
            dims=["time", "lev", "spatial"],
            coords={"time": [0], "lev": [0, 1, 2], "spatial": [0, 1]},
        )

        ps = xr.DataArray([[101325.0, 100000.0]], dims=["time", "spatial"])
        hya = xr.DataArray([0.0, 50000.0, 80000.0], dims=["lev"])
        hyb = xr.DataArray([1.0, 0.5, 0.2], dims=["lev"])

        new_levels = np.array([75000.0])

        # Should handle NaN values gracefully
        result = interp_hybrid_to_pressure(
            data=data, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels, lev_dim="lev"
        )

        assert result.shape == (1, 1, 2)
        # Where input had NaN, output should be NaN or interpolated from valid values
        # Exact behavior depends on implementation

    def test_extreme_pressure_ratios(self):
        """Test interpolation with extreme pressure ratios."""
        # Very high and very low pressures
        data = xr.DataArray(
            [[[200.0], [250.0], [300.0]]],
            dims=["x", "lev", "y"],
            coords={"x": [0], "lev": [0, 1, 2], "y": [0]},
        )

        ps = xr.DataArray([[100000.0]], dims=["x", "y"])
        # Extreme hybrid coordinates
        hya = xr.DataArray([0.0, 0.0, 0.0], dims=["lev"])
        hyb = xr.DataArray([1.0, 0.5, 0.1], dims=["lev"])

        new_levels = np.array([50000.0])

        result = interp_hybrid_to_pressure(
            data=data,
            ps=ps,
            hyam=hya,
            hybm=hyb,
            new_levels=new_levels,
            lev_dim="lev",
            method="log",
        )

        assert result.shape == (1, 1, 1)
        assert np.isfinite(result.values[0, 0, 0])


class TestInterpolationSpecialCases:
    """Test special cases and less common code paths."""

    def test_interp_multidim_with_datetime_coordinates(self):
        """Test multidimensional interpolation preserving datetime coordinates."""
        import pandas as pd

        # Create data with time coordinate
        time_coord = pd.date_range("2000-01-01", periods=3, freq="D")
        lat_in = np.array([0, 30, 60])
        lon_in = np.array([0, 90])

        data = np.random.randn(3, 3, 2)  # time, lat, lon

        data_in = xr.DataArray(
            data,
            dims=["time", "lat", "lon"],
            coords={"time": time_coord, "lat": lat_in, "lon": lon_in},
        )

        # Interpolate spatially but keep time
        lat_out = np.array([15, 45])
        lon_out = np.array([45])

        result = interp_multidim(data_in=data_in, lat_out=lat_out, lon_out=lon_out)

        # Should preserve time coordinate
        assert "time" in result.coords
        assert len(result.time) == 3
        assert result.shape == (3, 2, 1)  # time, lat_out, lon_out

    def test_sigma_to_hybrid_monotonic_check(self):
        """Test sigma to hybrid conversion with non-monotonic coordinates."""
        # Non-monotonic sigma coordinates (should still work)
        data = xr.DataArray(
            [[[250.0], [280.0], [260.0]]],  # 1x3x1
            dims=["time", "sigma", "spatial"],
            coords={"time": [0], "sigma": [0.2, 0.8, 0.6], "spatial": [0]},
        )

        sig_coords = xr.DataArray([0.2, 0.8, 0.6], dims=["sigma"])
        ps = xr.DataArray([[100000.0]], dims=["time", "spatial"])

        hya = xr.DataArray([0.0, 0.2], dims=["hlev"])
        hyb = xr.DataArray([0.8, 0.0], dims=["hlev"])

        # Should handle non-monotonic coordinates
        result = interp_sigma_to_hybrid(
            data=data, sig_coords=sig_coords, ps=ps, hyam=hya, hybm=hyb, lev_dim="sigma"
        )

        assert result.shape == (1, 2, 1)
        assert np.all(np.isfinite(result.values))

    def test_interpolation_with_unlimited_dimensions(self):
        """Test interpolation with datasets having unlimited dimensions."""
        # Simulate data that might come from NetCDF with unlimited time
        time_vals = np.arange(5)
        lev_vals = np.arange(3)

        # Create dataset similar to climate model output
        data = xr.DataArray(
            250 + np.random.randn(5, 3, 1, 1) * 20,
            dims=["time", "lev", "lat", "lon"],
            coords={"time": time_vals, "lev": lev_vals, "lat": [0], "lon": [0]},
        )

        ps = xr.DataArray(
            101325 + np.random.randn(5, 1, 1) * 1000,
            dims=["time", "lat", "lon"],
            coords={"time": time_vals, "lat": [0], "lon": [0]},
        )

        hya = xr.DataArray([0.0, 25000.0, 50000.0], dims=["lev"])
        hyb = xr.DataArray([1.0, 0.5, 0.0], dims=["lev"])

        new_levels = np.array([75000.0, 25000.0])

        result = interp_hybrid_to_pressure(
            data=data, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels, lev_dim="lev"
        )

        assert result.shape == (5, 2, 1, 1)  # time, plev, lat, lon
        assert "time" in result.dims
        assert_array_equal(result.time.values, time_vals)


class TestInterpolationGeoCatSyncCoverage:
    """Focused regression coverage for the GeoCAT sync paths."""

    def test_pressure_helpers_broadcast_with_default_xarray_dims(self):
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

    @pytest.mark.skipif(
        not PRIVATE_FUNCTIONS_AVAILABLE, reason="Private functions not available"
    )
    def test_temp_extrapolate_supports_legacy_and_explicit_t_bot_calls(self):
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

    def test_interp_hybrid_to_pressure_accepts_different_coeff_dimension_name(self):
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

    def test_interp_hybrid_to_pressure_temperature_extrapolation_uses_explicit_t_bot(
        self,
    ):
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
            new_levels=np.array([110000.0, 100000.0]),
            lev_dim="lev",
            extrapolate=True,
            variable="temperature",
            t_bot=t_bot,
            phi_sfc=phi_sfc,
        )

        assert out.shape == (2, 1, 1)
        assert np.all(np.isfinite(out.values))


class TestInterpolationHelperCoverage:
    """Extra helper and branch coverage for interpolation internals."""

    def test_interpolate_mb_wrapper_uses_requested_method(self):
        """The map-block helper should dispatch to the selected interpolation func."""
        data = np.array([280.0, 250.0, 220.0])
        curr_levels = np.array([1000.0, 700.0, 400.0])
        new_levels = np.array([850.0, 550.0])

        result = _interpolate_mb(data, curr_levels, new_levels, axis=0, method="linear")

        assert_array_almost_equal(result, np.array([265.0, 235.0]))

    def test_rename_colliding_coeff_dim_passthrough_for_non_dataarrays(self):
        """Non-xarray inputs should be returned unchanged by the rename guard."""
        hya = np.array([0.0, 0.2])
        hyb = np.array([1.0, 0.0])

        out_hya, out_hyb = _rename_colliding_coeff_dim(np.array([100000.0]), hya, hyb)

        assert out_hya is hya
        assert out_hyb is hyb

    def test_rename_colliding_coeff_dim_uses_internal_name_when_lev_exists(self):
        """A colliding coeff dim should fall back to the private helper dim name."""
        target = xr.DataArray(
            np.ones((3, 4)),
            dims=["dim_0", "lev"],
            coords={"dim_0": range(3), "lev": range(4)},
        )
        hya = xr.DataArray([0.0, 0.1], dims=["dim_0"])
        hyb = xr.DataArray([1.0, 0.0], dims=["dim_0"])

        out_hya, out_hyb = _rename_colliding_coeff_dim(target, hya, hyb)

        assert out_hya.dims == ("__hybrid_lev__",)
        assert out_hyb.dims == ("__hybrid_lev__",)

    def test_pressure_from_hybrid_rejects_invalid_input_types(self):
        """Hybrid-pressure helper should reject non-array and mixed-type inputs."""
        with pytest.raises(
            TypeError, match="must be xarray DataArrays or numpy arrays"
        ):
            pressure_at_hybrid_levels([100000.0], np.array([0.0]), np.array([1.0]))

        with pytest.raises(TypeError, match="must all be the same type"):
            pressure_at_hybrid_levels(
                xr.DataArray([100000.0], dims=["x"]),
                np.array([0.0]),
                np.array([1.0]),
            )

    def test_pressure_from_hybrid_numpy_validation_and_execution(self):
        """Numpy inputs should validate shape and broadcast into level-first output."""
        with pytest.raises(ValueError, match="dimension mismatch"):
            pressure_at_hybrid_levels(
                np.array([100000.0]), np.array([0.0, 0.2]), np.array([1.0])
            )

        with pytest.raises(ValueError, match="1-dimensional"):
            pressure_at_hybrid_levels(
                np.array([100000.0]),
                np.array([[0.0, 0.2]]),
                np.array([[1.0, 0.0]]),
            )

        pressure = pressure_at_hybrid_levels(
            np.array([100000.0, 90000.0]),
            np.array([0.0, 0.2]),
            np.array([1.0, 0.0]),
        )

        assert pressure.shape == (2, 2)
        assert_array_equal(
            pressure, np.array([[100000.0, 90000.0], [20000.0, 20000.0]])
        )

    def test_pressure_from_hybrid_warns_on_coeff_dim_name_mismatch(self):
        """Different xarray coeff dim names should warn and be renamed."""
        ps = xr.DataArray([100000.0, 95000.0], dims=["x"])
        hya = xr.DataArray([0.0, 0.2], dims=["a"])
        hyb = xr.DataArray([1.0, 0.0], dims=["b"])

        with pytest.warns(UserWarning, match="different dimension names"):
            pressure = pressure_at_hybrid_levels(ps, hya, hyb)

        assert pressure.dims == ("a", "x")
        assert pressure.shape == (2, 2)

    def test_delta_pressure_hybrid_rejects_invalid_inputs(self):
        """Layer-thickness helper should validate input types, shapes, and dimensions."""
        with pytest.raises(
            TypeError, match="Inputs must be xarray DataArrays or numpy arrays"
        ):
            delta_pressure_hybrid([100000.0], np.array([0.0]), np.array([1.0]))

        with pytest.raises(TypeError, match="p0 must be a scalar numeric value"):
            delta_pressure_hybrid(
                np.array([100000.0]), np.array([0.0]), np.array([1.0]), p0="bad"
            )

        with pytest.raises(ValueError, match="dimension mismatch"):
            delta_pressure_hybrid(
                np.array([100000.0]), np.array([0.0, 0.2]), np.array([1.0])
            )

        with pytest.raises(ValueError, match="1-dimensional"):
            delta_pressure_hybrid(
                np.array([100000.0]),
                np.array([[0.0, 0.2]]),
                np.array([[1.0, 0.0]]),
            )

    def test_delta_pressure_hybrid_numpy_and_dataarray_coeff_paths(self):
        """Delta-pressure helper should support both numpy and xarray pressure inputs."""
        hya = np.array([0.0, 0.2, 0.5])
        hyb = np.array([1.0, 0.5, 0.0])

        dph_np = delta_pressure_hybrid(np.array([100000.0, 90000.0]), hya, hyb)
        assert dph_np.shape == (2, 2)
        assert_array_equal(dph_np, np.array([[30000.0, 25000.0], [20000.0, 15000.0]]))

        dph_da = delta_pressure_hybrid(
            xr.DataArray([100000.0, 90000.0], dims=["x"]), hya, hyb
        )
        assert dph_da.dims == ("lev", "x")
        assert_array_equal(dph_da.values, dph_np)

    def test_delta_pressure_hybrid_warns_on_coeff_dim_name_mismatch(self):
        """Xarray coeff dim mismatches should be normalized with a warning."""
        ps = xr.DataArray([100000.0], dims=["x"])
        hya = xr.DataArray([0.0, 0.2], dims=["a"])
        hyb = xr.DataArray([1.0, 0.0], dims=["b"])

        with pytest.warns(UserWarning, match="different dimension names"):
            dph = delta_pressure_hybrid(ps, hya, hyb)

        assert dph.dims == ("a", "x")
        assert dph.shape == (1, 1)

    def test_sigma_from_hybrid_warns_on_coeff_dim_name_mismatch(self):
        """Sigma helper should also warn and rename mismatched coeff dims."""
        ps = xr.DataArray([101325.0], dims=["x"])
        hya = xr.DataArray([0.0, 0.2], dims=["a"])
        hyb = xr.DataArray([1.0, 0.0], dims=["b"])

        with pytest.warns(UserWarning, match="different dimension names"):
            sigma = _sigma_from_hybrid(ps, hya, hyb)

        assert sigma.dims == ("a", "x")
        assert sigma.shape == (2, 1)

    @pytest.mark.skipif(
        not PRIVATE_FUNCTIONS_AVAILABLE, reason="Private functions not available"
    )
    def test_temp_extrapolate_rejects_invalid_call_signature(self):
        """Temperature extrapolation should reject unsupported legacy signatures."""
        with pytest.raises(TypeError, match="accepts either"):
            _temp_extrapolate(xr.DataArray([280.0], dims=["lev"]), 100000.0)

    def test_pint_detection_and_strip_helpers(self):
        """Pint helper utilities should preserve expected quantity semantics."""
        pint = pytest.importorskip("pint")
        ureg = pint.UnitRegistry()
        quantity = xr.DataArray(np.array([1.0, 2.0]) * ureg.kelvin, dims=["x"])

        assert _is_pint_backed(np.array([1.0, 2.0])) is False
        assert _is_pint_backed(quantity) is True

        kept = _strip_unexpected_pint_units(quantity.copy(deep=False), in_pint=True)
        assert hasattr(kept.data, "magnitude")

        stripped = _strip_unexpected_pint_units(
            quantity.copy(deep=False), in_pint=False
        )
        assert isinstance(stripped.data, np.ndarray)
        assert_array_equal(stripped.values, np.array([1.0, 2.0]))

    def test_align_hybrid_level_dimension_validation_and_warning_paths(self):
        """Coefficient alignment helper should validate and rename as needed."""
        with pytest.raises(ValueError, match="dimension mismatch"):
            _align_hybrid_level_dimension(
                xr.DataArray([0.0, 0.2], dims=["a"]),
                xr.DataArray([1.0], dims=["a"]),
                "lev",
            )

        with pytest.raises(ValueError, match="one-dimensional arrays"):
            _align_hybrid_level_dimension(
                xr.DataArray([[0.0, 0.2]], dims=["a", "b"]),
                xr.DataArray([[1.0, 0.0]], dims=["a", "b"]),
                "lev",
            )

        with pytest.warns(UserWarning, match="different dimension names"):
            hyam, hybm = _align_hybrid_level_dimension(
                xr.DataArray([0.0, 0.2], dims=["a"]),
                xr.DataArray([1.0, 0.0], dims=["b"]),
                "lev",
            )

        assert hyam.dims == ("lev",)
        assert hybm.dims == ("lev",)

    def test_interp_hybrid_to_pressure_requires_dataarray_inputs(self):
        """Public hybrid-pressure interpolation should reject non-xarray inputs."""
        data = xr.DataArray(
            np.arange(6, dtype=float).reshape(3, 2),
            dims=["lev", "x"],
            coords={"lev": [0, 1, 2], "x": [0, 1]},
        )
        hyam = xr.DataArray([0.0, 0.2, 0.4], dims=["lev"])
        hybm = xr.DataArray([1.0, 0.5, 0.0], dims=["lev"])

        with pytest.raises(TypeError, match="must be xarray DataArray objects"):
            interp_hybrid_to_pressure(
                data=data,
                ps=np.full((2,), 100000.0),
                hyam=hyam,
                hybm=hybm,
                new_levels=np.array([95000.0]),
                lev_dim="lev",
            )

    def test_interp_hybrid_to_pressure_autodetects_vertical_with_cf_accessor(self):
        """Vertical axis should be auto-detected when cf_xarray metadata is present."""
        pytest.importorskip("cf_xarray")
        import cf_xarray  # noqa: F401

        data = xr.DataArray(
            np.arange(2 * 3 * 2, dtype=float).reshape(2, 3, 2),
            dims=["time", "model_level", "x"],
            coords={
                "time": [0, 1],
                "model_level": xr.DataArray(
                    [1, 2, 3], dims=["model_level"], attrs={"positive": "down"}
                ),
                "x": [0, 1],
            },
        )
        ps = xr.DataArray(
            np.full((2, 2), 100000.0),
            dims=["time", "x"],
            coords={"time": [0, 1], "x": [0, 1]},
        )
        hyam = xr.DataArray([0.0, 0.2, 0.4], dims=["hlev"])
        hybm = xr.DataArray([1.0, 0.5, 0.0], dims=["hlev"])

        result = interp_hybrid_to_pressure(
            data=data,
            ps=ps,
            hyam=hyam,
            hybm=hybm,
            new_levels=np.array([95000.0, 70000.0]),
        )

        assert result.dims == ("time", "plev", "x")
        assert result.shape == (2, 2, 2)

    def test_interp_hybrid_to_pressure_dask_single_vertical_chunk_preserves_metadata(
        self,
    ):
        """Single-chunk vertical dask inputs should use the xarray map-block path."""
        pytest.importorskip("dask")

        data = xr.DataArray(
            np.arange(2 * 3 * 2, dtype=float).reshape(2, 3, 2),
            dims=["time", "lev", "x"],
            coords={"time": [0, 1], "lev": [0, 1, 2], "x": [0, 1]},
            name="foo",
            attrs={"units": "K"},
        ).chunk({"time": 1, "lev": 3, "x": 2})
        ps = xr.DataArray(
            np.full((2, 2), 100000.0),
            dims=["time", "x"],
            coords={"time": [0, 1], "x": [0, 1]},
        ).chunk({"time": 1, "x": 2})
        hyam = xr.DataArray([0.0, 0.2, 0.4], dims=["lev"])
        hybm = xr.DataArray([1.0, 0.5, 0.0], dims=["lev"])

        result = interp_hybrid_to_pressure(
            data=data,
            ps=ps,
            hyam=hyam,
            hybm=hybm,
            new_levels=np.array([95000.0, 70000.0]),
            lev_dim="lev",
        )

        assert hasattr(result.data, "chunks")
        assert result.name == "foo"
        assert result.attrs["units"] == "K"
        assert result.compute().shape == (2, 2, 2)

    def test_interp_hybrid_to_pressure_dask_warns_for_chunked_vertical_dimension(self):
        """Multi-chunk vertical dask inputs should warn and use the fallback path."""
        pytest.importorskip("dask")

        data = xr.DataArray(
            np.arange(2 * 3 * 2, dtype=float).reshape(2, 3, 2),
            dims=["time", "lev", "x"],
            coords={"time": [0, 1], "lev": [0, 1, 2], "x": [0, 1]},
        ).chunk({"time": 1, "lev": 1, "x": 2})
        ps = xr.DataArray(
            np.full((2, 2), 100000.0),
            dims=["time", "x"],
            coords={"time": [0, 1], "x": [0, 1]},
        ).chunk({"time": 1, "x": 2})
        hyam = xr.DataArray([0.0, 0.2, 0.4], dims=["lev"])
        hybm = xr.DataArray([1.0, 0.5, 0.0], dims=["lev"])

        with pytest.warns(UserWarning, match="Chunking along lev"):
            result = interp_hybrid_to_pressure(
                data=data,
                ps=ps,
                hyam=hyam,
                hybm=hybm,
                new_levels=np.array([95000.0, 70000.0]),
                lev_dim="lev",
            )

        assert hasattr(result.data, "chunks")
        assert result.compute().shape == (2, 2, 2)

    def test_interp_hybrid_to_pressure_chunks_eager_data_when_other_inputs_are_dask(
        self,
    ):
        """Fallback dask path should chunk eager data when companion inputs are dask-backed."""
        pytest.importorskip("dask")

        data = xr.DataArray(
            np.arange(2 * 3 * 2, dtype=float).reshape(2, 3, 2),
            dims=["time", "lev", "x"],
            coords={"time": [0, 1], "lev": [0, 1, 2], "x": [0, 1]},
        )
        ps = xr.DataArray(
            np.full((2, 2), 100000.0),
            dims=["time", "x"],
            coords={"time": [0, 1], "x": [0, 1]},
        ).chunk({"time": 1, "x": 2})
        hyam = xr.DataArray([0.0, 0.2, 0.4], dims=["lev"])
        hybm = xr.DataArray([1.0, 0.5, 0.0], dims=["lev"])

        result = interp_hybrid_to_pressure(
            data=data,
            ps=ps,
            hyam=hyam,
            hybm=hybm,
            new_levels=np.array([95000.0, 70000.0]),
            lev_dim="lev",
        )

        assert hasattr(result.data, "chunks")
        assert result.compute().shape == (2, 2, 2)

    def test_interp_hybrid_to_pressure_preserves_metadata_when_xarray_map_blocks_succeeds(
        self, monkeypatch
    ):
        """When xr.map_blocks returns a DataArray directly, metadata should be copied."""
        pytest.importorskip("dask")

        data = xr.DataArray(
            np.arange(2 * 3 * 2, dtype=float).reshape(2, 3, 2),
            dims=["time", "lev", "x"],
            coords={"time": [0, 1], "lev": [0, 1, 2], "x": [0, 1]},
            name="foo",
            attrs={"units": "K", "long_name": "temperature"},
        ).chunk({"time": 1, "lev": 3, "x": 2})
        ps = xr.DataArray(
            np.full((2, 2), 100000.0),
            dims=["time", "x"],
            coords={"time": [0, 1], "x": [0, 1]},
        )
        hyam = xr.DataArray([0.0, 0.2, 0.4], dims=["lev"])
        hybm = xr.DataArray([1.0, 0.5, 0.0], dims=["lev"])

        def fake_map_blocks(func, arr, args=(), kwargs=None, template=None):
            del func, arr, args, kwargs, template
            return xr.DataArray(
                np.zeros((2, 2, 2)),
                dims=["time", "plev", "x"],
                coords={"time": [0, 1], "plev": [95000.0, 70000.0], "x": [0, 1]},
            )

        monkeypatch.setattr(xr, "map_blocks", fake_map_blocks)

        result = interp_hybrid_to_pressure(
            data=data,
            ps=ps,
            hyam=hyam,
            hybm=hybm,
            new_levels=np.array([95000.0, 70000.0]),
            lev_dim="lev",
        )

        assert result.name == "foo"
        assert result.attrs["units"] == "K"
        assert result.attrs["long_name"] == "temperature"
        assert result.dims == ("time", "plev", "x")
        assert result.shape == (2, 2, 2)

    def test_interp_sigma_to_hybrid_autodetects_vertical_with_cf_accessor(self):
        """Sigma interpolation should also honor cf_xarray vertical metadata."""
        pytest.importorskip("cf_xarray")
        import cf_xarray  # noqa: F401

        data = xr.DataArray(
            np.arange(2 * 3 * 2, dtype=float).reshape(2, 3, 2),
            dims=["time", "model_level", "x"],
            coords={
                "time": [0, 1],
                "model_level": xr.DataArray(
                    [0.2, 0.6, 1.0], dims=["model_level"], attrs={"positive": "down"}
                ),
                "x": [0, 1],
            },
        )
        sig_coords = xr.DataArray([0.2, 0.6, 1.0], dims=["model_level"])
        ps = xr.DataArray(
            np.full((2, 2), 100000.0),
            dims=["time", "x"],
            coords={"time": [0, 1], "x": [0, 1]},
        )
        hyam = xr.DataArray([0.0, 0.2, 0.4], dims=["hlev"])
        hybm = xr.DataArray([1.0, 0.5, 0.0], dims=["hlev"])

        result = interp_sigma_to_hybrid(
            data=data,
            sig_coords=sig_coords,
            ps=ps,
            hyam=hyam,
            hybm=hybm,
        )

        assert result.dims == ("time", "hlev", "x")
        assert result.shape == (2, 3, 2)

    def test_interp_sigma_to_hybrid_requires_explicit_vertical_when_cf_fails(self):
        """Without lev_dim or CF metadata, sigma interpolation should fail clearly."""
        data = xr.DataArray(
            np.arange(2 * 3 * 2, dtype=float).reshape(2, 3, 2),
            dims=["time", "model_level", "x"],
            coords={"time": [0, 1], "model_level": [0.2, 0.6, 1.0], "x": [0, 1]},
        )
        sig_coords = xr.DataArray([0.2, 0.6, 1.0], dims=["model_level"])
        ps = xr.DataArray(
            np.full((2, 2), 100000.0),
            dims=["time", "x"],
            coords={"time": [0, 1], "x": [0, 1]},
        )
        hyam = xr.DataArray([0.0, 0.2, 0.4], dims=["hlev"])
        hybm = xr.DataArray([1.0, 0.5, 0.0], dims=["hlev"])

        with pytest.raises(
            ValueError, match="Unable to determine vertical dimension name"
        ):
            interp_sigma_to_hybrid(
                data=data,
                sig_coords=sig_coords,
                ps=ps,
                hyam=hyam,
                hybm=hybm,
            )


# Performance and stress tests
@pytest.mark.slow
class TestInterpolationStress:
    """Stress tests for interpolation functions with large or complex data."""

    def test_large_dataset_interpolation(self):
        """Test interpolation with large datasets to check memory usage."""
        # Large dataset dimensions
        time, lev, lat, lon = 50, 20, 90, 180

        # Create data in chunks to avoid memory issues
        data = xr.DataArray(
            250 + 20 * np.random.randn(time, lev, lat, lon),
            dims=["time", "lev", "lat", "lon"],
            coords={
                "time": np.arange(time),
                "lev": np.arange(lev),
                "lat": np.linspace(-89, 89, lat),
                "lon": np.linspace(0, 358, lon),
            },
        )

        ps = xr.DataArray(
            101325 + 2000 * np.random.randn(time, lat, lon), dims=["time", "lat", "lon"]
        )

        hya = xr.DataArray(np.linspace(0, 50000, lev), dims=["lev"])
        hyb = xr.DataArray(np.linspace(1.0, 0.0, lev), dims=["lev"])

        new_levels = np.array([100000, 70000, 50000, 30000, 10000])

        # Should complete without memory errors
        result = interp_hybrid_to_pressure(
            data=data, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels, lev_dim="lev"
        )

        assert result.shape == (50, 5, 90, 180)

        # Check some basic properties
        assert np.all(np.isfinite(result.values[~np.isnan(result.values)]))

    def test_high_resolution_spatial_interpolation(self):
        """Test multidimensional interpolation with high resolution output."""
        # Coarse input grid
        lat_in = np.linspace(-60, 60, 25)
        lon_in = np.linspace(0, 360, 50, endpoint=False)

        data = np.random.randn(25, 50)
        data_in = xr.DataArray(
            data, dims=["lat", "lon"], coords={"lat": lat_in, "lon": lon_in}
        )

        # Fine output grid
        lat_out = np.linspace(-59, 59, 100)
        lon_out = np.linspace(0, 359, 200)

        result = interp_multidim(data_in=data_in, lat_out=lat_out, lon_out=lon_out)

        assert result.shape == (100, 200)
        assert np.all(np.isfinite(result.values[~np.isnan(result.values)]))
        assert np.isfinite(result.values).mean() > 0.95


if __name__ == "__main__":
    # Run enhanced tests
    pytest.main([__file__, "-v", "--tb=short"])
