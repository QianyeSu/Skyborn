"""
Enhanced tests for skyborn.interp.interpolation module.

This module provides additional test coverage for edge cases, performance scenarios,
and specific code paths that may not be covered by the existing test suite.
Focuses on improving code coverage for the interpolation functionality.
"""

import pytest
import numpy as np
import xarray as xr
from numpy.testing import assert_array_almost_equal, assert_array_equal
import warnings

from skyborn.interp.interpolation import (
    interp_hybrid_to_pressure,
    interp_sigma_to_hybrid,
    interp_multidim,
    _pressure_from_hybrid,
    _sigma_from_hybrid,
    _func_interpolate,
    __pres_lev_mandatory__,
)

# Try to import private functions - they may not be available
try:
    from skyborn.interp.interpolation import (
        _pre_interp_multidim,
        _post_interp_multidim,
        _vertical_remap,
        _temp_extrapolate,
        _geo_height_extrapolate,
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
        hya = xr.DataArray([0.0, 50000.0])  # 2 levels
        hyb = xr.DataArray([1.0, 0.0])
        p0 = 100000.0

        pressure = _pressure_from_hybrid(ps, hya, hyb, p0)

        assert pressure.shape == (2,)
        assert pressure[0] == 100000.0  # Surface level
        assert pressure[1] == 50000.0  # Top level

    def test_pressure_from_hybrid_zero_surface_pressure(self):
        """Test pressure calculation with zero surface pressure."""
        ps = xr.DataArray([0.0, 50000.0])
        hya = xr.DataArray([0.0, 25000.0])
        hyb = xr.DataArray([1.0, 0.5])
        p0 = 100000.0

        pressure = _pressure_from_hybrid(ps, hya, hyb, p0)

        # First point should have zero pressure at surface
        assert pressure[0, 0] == 0.0
        assert pressure[1, 0] > 0.0  # Second point should be positive

    def test_sigma_from_hybrid_edge_values(self):
        """Test sigma calculation with edge values."""
        ps = xr.DataArray([101325.0])
        # Test with hya=0 (pure sigma) and hyb=0 (pure pressure)
        hya = xr.DataArray([0.0, 50000.0, 101325.0])
        hyb = xr.DataArray([1.0, 0.0, 0.0])
        p0 = 101325.0

        sigma = _sigma_from_hybrid(ps, hya, hyb, p0)

        assert sigma.shape == (1, 3)
        assert_array_almost_equal(sigma[0, 0], 1.0, decimal=10)  # Pure sigma level
        assert_array_almost_equal(
            sigma[0, 1], 50000.0 / 101325.0, decimal=6
        )  # Pure pressure
        # At reference pressure
        assert_array_almost_equal(sigma[0, 2], 1.0, decimal=10)

    def test_func_interpolate_function_properties(self):
        """Test properties of interpolation functions."""
        # Test linear interpolation
        func_linear = _func_interpolate("linear")

        # Create test data
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([10.0, 20.0, 30.0])
        xi = np.array([1.5, 2.5])

        result = func_linear(x, y, xi)
        expected = np.array([15.0, 25.0])  # Linear interpolation
        assert_array_almost_equal(result, expected)

        # Test log interpolation
        func_log = _func_interpolate("log")

        # Positive values for log interpolation
        x_log = np.array([1.0, 10.0, 100.0])
        y_log = np.array([1.0, 10.0, 100.0])
        xi_log = np.array([3.16227766])  # sqrt(10)

        result_log = func_log(x_log, y_log, xi_log)
        assert len(result_log) == 1
        assert result_log[0] > 1.0 and result_log[0] < 100.0

    def test_hybrid_to_pressure_single_level(self):
        """Test hybrid to pressure interpolation with single vertical level."""
        # Single level data
        data = xr.DataArray(
            [[280.0, 285.0], [275.0, 290.0]],  # 2x2 spatial grid
            dims=["lat", "lon"],
            coords={"lat": [0, 30], "lon": [0, 90]},
        )

        ps = xr.DataArray([101325.0, 100000.0], dims=["lat"], coords={"lat": [0, 30]})

        # Single hybrid level (surface)
        hya = xr.DataArray([0.0], dims=["lev"])
        hyb = xr.DataArray([1.0], dims=["lev"])

        # Add level dimension to data
        data_3d = data.expand_dims({"lev": [0]})

        new_levels = np.array([101325.0, 95000.0])

        result = interp_hybrid_to_pressure(
            data=data_3d, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels
        )

        assert result.shape == (2, 2, 2)  # plev, lat, lon
        # At surface pressure, should match original data
        assert_array_almost_equal(result[0, :, 0], [280.0, 275.0])

    def test_hybrid_to_pressure_temperature_extrapolation(self):
        """Test temperature extrapolation in hybrid to pressure interpolation."""
        # Create test data with temperature
        nlev = 5
        data = xr.DataArray(
            np.array([[[300.0], [280.0], [260.0], [240.0], [220.0]]]).T,  # 1x5x1
            dims=["lat", "lev", "lon"],
            coords={"lat": [0], "lev": range(nlev), "lon": [0]},
        )

        ps = xr.DataArray([101325.0], dims=["lat"], coords={"lat": [0]})
        hya = xr.DataArray(np.linspace(0, 50000, nlev), dims=["lev"])
        hyb = xr.DataArray(np.linspace(1.0, 0.0, nlev), dims=["lev"])

        # Request levels both within and outside the data range
        # surface, mid, top
        new_levels = np.array([110000.0, 70000.0, 20000.0])

        result = interp_hybrid_to_pressure(
            data=data,
            ps=ps,
            hyam=hya,
            hybm=hyb,
            new_levels=new_levels,
            extrapolate=True,
            variable="temperature",
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

        ps = xr.DataArray([100000.0], dims=["lat"])
        hya = xr.DataArray(np.linspace(0, 30000, nlev), dims=["lev"])
        hyb = xr.DataArray(np.linspace(1.0, 0.0, nlev), dims=["lev"])

        new_levels = np.array([105000.0, 50000.0, 10000.0])

        result = interp_hybrid_to_pressure(
            data=data,
            ps=ps,
            hyam=hya,
            hybm=hyb,
            new_levels=new_levels,
            extrapolate=True,
            variable="geopotential",
        )

        assert result.shape == (1, 3, 1)
        # Geopotential should increase with altitude (decrease with pressure)
        assert result[0, 0, 0] <= result[0, 1, 0] <= result[0, 2, 0]

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
        hya = xr.DataArray([0.0, 50000.0], dims=["hlev"])
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
        hya = xr.DataArray([0.0, 1.0, 100.0])  # Small increments
        hyb = xr.DataArray([1.0, 0.99999, 0.999])
        p0 = 101325.0

        pressure = _pressure_from_hybrid(ps, hya, hyb, p0)

        # Check for reasonable precision
        assert np.abs(pressure[0, 0] - 101325.0) < 1e-10
        assert pressure[0, 1] < pressure[0, 0]  # Should be decreasing
        assert pressure[0, 2] < pressure[0, 1]

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
            data=data, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels
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
        hya = xr.DataArray([0.0, 99000.0, 99999.0], dims=["lev"])
        hyb = xr.DataArray([1.0, 0.01, 0.00001], dims=["lev"])

        new_levels = np.array([50000.0])

        result = interp_hybrid_to_pressure(
            data=data, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels, method="log"
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

        hya = xr.DataArray([0.0, 50000.0], dims=["hlev"])
        hyb = xr.DataArray([1.0, 0.0], dims=["hlev"])

        # Should handle non-monotonic coordinates
        result = interp_sigma_to_hybrid(
            data=data, sig_coords=sig_coords, ps=ps, hyam=hya, hybm=hyb
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
            data=data, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels
        )

        assert result.shape == (5, 2, 1, 1)  # time, plev, lat, lon
        assert "time" in result.dims
        assert_array_equal(result.time.values, time_vals)


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
            data=data, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels
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
        assert np.all(np.isfinite(result.values))


if __name__ == "__main__":
    # Run enhanced tests
    pytest.main([__file__, "-v", "--tb=short"])
