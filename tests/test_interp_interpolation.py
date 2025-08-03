"""
Tests for skyborn.interp.interpolation module.

This module tests the interpolation functionality including hybrid-sigma
to pressure level interpolation and multidimensional spatial interpolation.
"""

import pytest
import numpy as np
import xarray as xr
from numpy.testing import assert_array_almost_equal, assert_array_equal

from skyborn.interp.interpolation import (
    interp_hybrid_to_pressure,
    interp_sigma_to_hybrid,
    interp_multidim,
    _pressure_from_hybrid,
    _sigma_from_hybrid,
    _func_interpolate,
)


class TestInterpolationHelpers:
    """Test helper functions for interpolation."""

    def test_pressure_from_hybrid(self):
        """Test pressure calculation from hybrid coordinates."""
        # Simple test data
        ps = np.array([100000, 95000])  # Pa
        hya = np.array([0.0, 0.1, 0.5])  # hybrid A coefficients
        hyb = np.array([1.0, 0.9, 0.5])  # hybrid B coefficients
        p0 = 100000.0

        # Convert to xarray
        ps_da = xr.DataArray(ps, dims=["x"])
        hya_da = xr.DataArray(hya, dims=["lev"])
        hyb_da = xr.DataArray(hyb, dims=["lev"])

        pressure = _pressure_from_hybrid(ps_da, hya_da, hyb_da, p0)

        # Check shape and basic properties
        assert pressure.shape == (2, 3)  # x, lev
        assert np.all(pressure > 0)  # All pressures should be positive
        assert np.all(pressure <= 100000)  # Should not exceed reference pressure

    def test_sigma_from_hybrid(self):
        """Test sigma calculation from hybrid coordinates."""
        ps = np.array([100000, 95000])  # Pa
        hya = np.array([0.0, 0.1, 0.5])
        hyb = np.array([1.0, 0.9, 0.5])
        p0 = 100000.0

        # Convert to xarray
        ps_da = xr.DataArray(ps, dims=["x"])
        hya_da = xr.DataArray(hya, dims=["lev"])
        hyb_da = xr.DataArray(hyb, dims=["lev"])

        sigma = _sigma_from_hybrid(ps_da, hya_da, hyb_da, p0)

        # Check shape and basic properties
        assert sigma.shape == (2, 3)  # x, lev
        assert np.all(sigma >= 0)  # Sigma should be non-negative
        assert np.all(sigma <= 1.2)  # Should be close to [0, 1] range

    def test_func_interpolate(self):
        """Test interpolation function selection."""
        # Test linear interpolation function
        func_linear = _func_interpolate("linear")
        assert func_linear is not None

        # Test log interpolation function
        func_log = _func_interpolate("log")
        assert func_log is not None

        # Test invalid method
        with pytest.raises(ValueError, match="Unknown interpolation method"):
            _func_interpolate("invalid_method")


class TestHybridToPressureInterpolation:
    """Test hybrid-sigma to pressure level interpolation."""

    @pytest.fixture
    def sample_hybrid_data(self):
        """Create sample hybrid-sigma data for testing."""
        # Dimensions
        time = 5
        lev = 10
        lat = 20
        lon = 30

        # Coordinates
        time_coord = np.arange(time)
        lev_coord = np.arange(lev)
        lat_coord = np.linspace(-90, 90, lat)
        lon_coord = np.linspace(0, 357.5, lon)

        # Create realistic hybrid coefficients
        hya = np.linspace(0, 50000, lev)  # Pa
        hyb = np.linspace(1.0, 0.0, lev)  # dimensionless

        # Surface pressure (varying in space and time)
        ps_base = 101325.0  # Standard atmospheric pressure
        ps = ps_base + np.random.randn(time, lat, lon) * 1000

        # Sample temperature data
        temp_data = 250 + 50 * np.random.randn(time, lev, lat, lon)

        # Create xarray objects
        data = xr.DataArray(
            temp_data,
            dims=["time", "lev", "lat", "lon"],
            coords={
                "time": time_coord,
                "lev": lev_coord,
                "lat": lat_coord,
                "lon": lon_coord,
            },
            attrs={"units": "K", "long_name": "Temperature"},
        )

        ps_da = xr.DataArray(
            ps,
            dims=["time", "lat", "lon"],
            coords={"time": time_coord, "lat": lat_coord, "lon": lon_coord},
            attrs={"units": "Pa", "long_name": "Surface Pressure"},
        )

        hya_da = xr.DataArray(hya, dims=["lev"], coords={"lev": lev_coord})
        hyb_da = xr.DataArray(hyb, dims=["lev"], coords={"lev": lev_coord})

        return data, ps_da, hya_da, hyb_da

    def test_interp_hybrid_to_pressure_basic(self, sample_hybrid_data):
        """Test basic hybrid to pressure interpolation."""
        data, ps, hya, hyb = sample_hybrid_data

        # Use a subset of standard pressure levels
        new_levels = np.array([100000, 85000, 70000, 50000, 30000])  # Pa

        result = interp_hybrid_to_pressure(
            data=data, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels, lev_dim="lev"
        )

        # Check output structure
        assert "plev" in result.dims
        assert "lev" not in result.dims
        assert len(result.plev) == len(new_levels)
        assert result.shape == (5, 5, 20, 30)  # time, plev, lat, lon

        # Check that pressure coordinates are correct
        assert_array_equal(result.plev.values, new_levels)

        # Check that metadata is preserved
        assert result.attrs["units"] == "K"
        assert result.attrs["long_name"] == "Temperature"

    def test_interp_hybrid_to_pressure_methods(self, sample_hybrid_data):
        """Test different interpolation methods."""
        data, ps, hya, hyb = sample_hybrid_data
        new_levels = np.array([100000, 50000, 30000])

        # Test linear interpolation
        result_linear = interp_hybrid_to_pressure(
            data=data, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels, method="linear"
        )

        # Test log interpolation
        result_log = interp_hybrid_to_pressure(
            data=data, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels, method="log"
        )

        # Both should have same shape
        assert result_linear.shape == result_log.shape

        # Results should be different (unless by coincidence)
        assert not np.allclose(result_linear.values, result_log.values)

    def test_interp_hybrid_to_pressure_extrapolation(self, sample_hybrid_data):
        """Test extrapolation functionality."""
        data, ps, hya, hyb = sample_hybrid_data
        new_levels = np.array([100000, 85000, 70000])

        # Test with extrapolation enabled
        result = interp_hybrid_to_pressure(
            data=data,
            ps=ps,
            hyam=hya,
            hybm=hyb,
            new_levels=new_levels,
            extrapolate=True,
            variable="other",  # Use simple extrapolation
        )

        assert result.shape == (5, 3, 20, 30)
        assert np.all(np.isfinite(result.values))

    def test_interp_hybrid_to_pressure_validation(self, sample_hybrid_data):
        """Test input validation."""
        data, ps, hya, hyb = sample_hybrid_data
        new_levels = np.array([100000, 50000])

        # Test invalid method
        with pytest.raises(ValueError, match="Unknown interpolation method"):
            interp_hybrid_to_pressure(
                data=data,
                ps=ps,
                hyam=hya,
                hybm=hyb,
                new_levels=new_levels,
                method="invalid",
            )

        # Test extrapolation without variable
        with pytest.raises(ValueError, match="If `extrapolate` is True"):
            interp_hybrid_to_pressure(
                data=data,
                ps=ps,
                hyam=hya,
                hybm=hyb,
                new_levels=new_levels,
                extrapolate=True,
            )

        # Test invalid variable
        with pytest.raises(ValueError, match="accepted values are"):
            interp_hybrid_to_pressure(
                data=data,
                ps=ps,
                hyam=hya,
                hybm=hyb,
                new_levels=new_levels,
                extrapolate=True,
                variable="invalid_variable",
            )


class TestSigmaToHybridInterpolation:
    """Test sigma to hybrid coordinate interpolation."""

    @pytest.fixture
    def sample_sigma_data(self):
        """Create sample sigma coordinate data."""
        # Dimensions
        time = 3
        sigma_lev = 8
        lat = 15
        lon = 20

        # Sigma coordinates (0 at top, 1 at surface)
        sigma_coords = np.linspace(0.1, 1.0, sigma_lev)

        # Sample data
        data = 250 + 50 * np.random.randn(time, sigma_lev, lat, lon)

        # Surface pressure
        ps = 101325 + np.random.randn(time, lat, lon) * 1000

        # Target hybrid coefficients
        nlev_hybrid = 6
        hya = np.linspace(0, 30000, nlev_hybrid)
        hyb = np.linspace(1.0, 0.0, nlev_hybrid)

        # Create xarray objects
        data_da = xr.DataArray(
            data,
            dims=["time", "sigma", "lat", "lon"],
            coords={
                "time": np.arange(time),
                "sigma": sigma_coords,
                "lat": np.linspace(-45, 45, lat),
                "lon": np.linspace(0, 357.5, lon),
            },
        )

        ps_da = xr.DataArray(ps, dims=["time", "lat", "lon"], coords=data_da.coords)

        hya_da = xr.DataArray(hya, dims=["hlev"])
        hyb_da = xr.DataArray(hyb, dims=["hlev"])
        sig_da = xr.DataArray(sigma_coords, dims=["sigma"])

        return data_da, sig_da, ps_da, hya_da, hyb_da

    def test_interp_sigma_to_hybrid_basic(self, sample_sigma_data):
        """Test basic sigma to hybrid interpolation."""
        data, sig_coords, ps, hya, hyb = sample_sigma_data

        result = interp_sigma_to_hybrid(
            data=data, sig_coords=sig_coords, ps=ps, hyam=hya, hybm=hyb, lev_dim="sigma"
        )

        # Check output structure
        assert "hlev" in result.dims
        assert "sigma" not in result.dims
        assert len(result.hlev) == len(hya)
        assert result.shape == (3, 6, 15, 20)  # time, hlev, lat, lon


class TestMultidimensionalInterpolation:
    """Test multidimensional spatial interpolation."""

    def test_interp_multidim_basic(self):
        """Test basic multidimensional interpolation."""
        # Create test data
        lat_in = np.array([0, 30, 60, 90])
        lon_in = np.array([0, 90, 180, 270])
        data = np.random.randn(4, 4)

        data_in = xr.DataArray(
            data, dims=["lat", "lon"], coords={"lat": lat_in, "lon": lon_in}
        )

        # Output coordinates
        lat_out = np.array([15, 45, 75])
        lon_out = np.array([45, 135, 225, 315])

        result = interp_multidim(data_in=data_in, lat_out=lat_out, lon_out=lon_out)

        # Check output shape
        assert result.shape == (3, 4)  # lat_out, lon_out
        assert_array_equal(result.lat.values, lat_out)
        assert_array_equal(result.lon.values, lon_out)

    def test_interp_multidim_numpy_input(self):
        """Test multidimensional interpolation with numpy arrays."""
        lat_in = np.array([0, 30, 60])
        lon_in = np.array([0, 120, 240])
        data = np.random.randn(3, 3)

        lat_out = np.array([15, 45])
        lon_out = np.array([60, 180])

        result = interp_multidim(
            data_in=data, lat_in=lat_in, lon_in=lon_in, lat_out=lat_out, lon_out=lon_out
        )

        assert result.shape == (2, 2)
        assert isinstance(result, xr.DataArray)

    def test_interp_multidim_cyclic(self):
        """Test multidimensional interpolation with cyclic boundary."""
        lat_in = np.array([-90, 0, 90])
        lon_in = np.array([0, 180])  # Only half the globe
        data = np.array([[1, 2], [3, 4], [5, 6]])

        data_in = xr.DataArray(
            data, dims=["lat", "lon"], coords={"lat": lat_in, "lon": lon_in}
        )

        # Request data at 360 degrees (should wrap to 0)
        lat_out = np.array([0])
        lon_out = np.array([360])

        result = interp_multidim(
            data_in=data_in, lat_out=lat_out, lon_out=lon_out, cyclic=True
        )

        assert result.shape == (1, 1)
        # Should be close to the value at lon=0
        assert not np.isnan(result.values[0, 0])

    def test_interp_multidim_missing_values(self):
        """Test handling of missing values."""
        lat_in = np.array([0, 30, 60])
        lon_in = np.array([0, 90, 180])
        data = np.array([[1, 2, 3], [4, 99, 6], [7, 8, 9]])  # 99 is missing

        data_in = xr.DataArray(
            data, dims=["lat", "lon"], coords={"lat": lat_in, "lon": lon_in}
        )

        lat_out = np.array([15, 45])
        lon_out = np.array([45, 135])

        result = interp_multidim(
            data_in=data_in, lat_out=lat_out, lon_out=lon_out, missing_val=99
        )

        assert result.shape == (2, 2)

    def test_interp_multidim_validation(self):
        """Test input validation for multidimensional interpolation."""
        data = np.random.randn(3, 3)
        lat_out = np.array([15, 45])
        lon_out = np.array([60, 180])

        # Test missing coordinates for numpy input
        with pytest.raises(ValueError, match="lat_in and lon_in must be provided"):
            interp_multidim(data_in=data, lat_out=lat_out, lon_out=lon_out)


class TestInterpolationIntegration:
    """Integration tests for interpolation module."""

    def test_interpolation_with_climate_data(self, sample_climate_data):
        """Test interpolation with realistic climate data."""
        temp = sample_climate_data["temperature"]

        # Create fake hybrid coordinate data
        nlev = 10
        ps = xr.DataArray(
            101325 + np.random.randn(12, 73, 144) * 1000,
            dims=["time", "lat", "lon"],
            coords=temp.coords,
        )

        # Add hybrid level dimension to temperature
        temp_hybrid = xr.concat([temp] * nlev, dim="lev")
        temp_hybrid = temp_hybrid.assign_coords(lev=np.arange(nlev))

        hya = xr.DataArray(np.linspace(0, 50000, nlev), dims=["lev"])
        hyb = xr.DataArray(np.linspace(1.0, 0.0, nlev), dims=["lev"])

        # Test hybrid to pressure interpolation
        new_levels = np.array([100000, 85000, 70000, 50000])
        result = interp_hybrid_to_pressure(
            data=temp_hybrid, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels
        )

        assert result.shape == (12, 4, 73, 144)  # time, plev, lat, lon
        assert np.all(np.isfinite(result.values))

    def test_interpolation_error_handling(self):
        """Test comprehensive error handling."""
        # Create minimal test data
        data = xr.DataArray(
            np.random.randn(3, 5, 5),
            dims=["lev", "lat", "lon"],
            coords={
                "lev": np.arange(3),
                "lat": np.linspace(-60, 60, 5),
                "lon": np.linspace(0, 288, 5),
            },
        )

        ps = xr.DataArray(
            101325 + np.random.randn(5, 5),
            dims=["lat", "lon"],
            coords={"lat": data.lat, "lon": data.lon},
        )

        hya = xr.DataArray([0, 25000, 50000], dims=["lev"])
        hyb = xr.DataArray([1.0, 0.5, 0.0], dims=["lev"])

        # Test with invalid pressure levels (negative)
        with pytest.raises((ValueError, RuntimeError)):
            interp_hybrid_to_pressure(
                data=data,
                ps=ps,
                hyam=hya,
                hybm=hyb,
                new_levels=np.array([-1000, 50000]),
            )


# Performance tests (marked as slow)
@pytest.mark.slow
class TestInterpolationPerformance:
    """Performance tests for interpolation module."""

    def test_hybrid_to_pressure_large_data(self):
        """Test hybrid to pressure interpolation with large datasets."""
        # Large dataset
        time, lev, lat, lon = 100, 50, 180, 360

        data = xr.DataArray(
            250 + 50 * np.random.randn(time, lev, lat, lon),
            dims=["time", "lev", "lat", "lon"],
            coords={
                "time": np.arange(time),
                "lev": np.arange(lev),
                "lat": np.linspace(-90, 90, lat),
                "lon": np.linspace(0, 357.5, lon),
            },
        )

        ps = xr.DataArray(
            101325 + np.random.randn(time, lat, lon) * 1000,
            dims=["time", "lat", "lon"],
            coords=data.coords,
        )

        hya = xr.DataArray(np.linspace(0, 50000, lev), dims=["lev"])
        hyb = xr.DataArray(np.linspace(1.0, 0.0, lev), dims=["lev"])

        new_levels = np.array([100000, 85000, 70000, 50000, 30000])

        # Should complete without memory issues
        result = interp_hybrid_to_pressure(
            data=data, ps=ps, hyam=hya, hybm=hyb, new_levels=new_levels
        )

        assert result.shape == (100, 5, 180, 360)
        assert np.all(np.isfinite(result.values))


if __name__ == "__main__":
    # Quick test runner
    pytest.main([__file__, "-v"])
