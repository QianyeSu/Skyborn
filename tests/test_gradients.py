"""
Tests for skyborn.gradients module.

This module tests the gradient calculation functionality for atmospheric
and oceanic data on spherical coordinates.
"""

import pytest
import numpy as np
import xarray as xr
from numpy.testing import assert_array_almost_equal, assert_array_equal

from skyborn.gradients import (
    calculate_gradient,
    calculate_meridional_gradient,
    calculate_zonal_gradient,
    calculate_vertical_gradient,
    EARTH_RADIUS,
)


class TestCalculateGradient:
    """Test the base gradient calculation function."""

    def test_calculate_gradient_1d(self):
        """Test gradient calculation for 1D array."""
        field = np.array([1, 4, 9, 16, 25])  # x^2 where x = [1,2,3,4,5]
        coordinates = np.array([1, 2, 3, 4, 5])

        gradient = calculate_gradient(field, coordinates, axis=0)

        # For x^2, gradient should be approximately 2*x
        # At x=1: gradient ≈ 2*1 = 2, at x=2: gradient ≈ 2*2 = 4, etc.
        expected = np.array([3, 5, 7, 9, 7])  # Numerical approximation

        assert gradient.shape == field.shape
        assert np.allclose(gradient, expected, atol=1.0)

    def test_calculate_gradient_2d(self):
        """Test gradient calculation for 2D array."""
        # Create a simple 2D field
        field = np.array([[1, 2, 3], [2, 4, 6], [3, 6, 9]])
        coordinates = np.array([0, 1, 2])

        # Calculate gradient along axis 1 (columns)
        gradient = calculate_gradient(field, coordinates, axis=1)

        assert gradient.shape == field.shape
        assert np.all(np.isfinite(gradient))

        # For a linear field, gradient should be approximately constant
        # Check that interior points have similar gradients
        assert np.allclose(gradient[:, 1], gradient[:, 1][0], rtol=0.1)

    def test_calculate_gradient_multidimensional(self):
        """Test gradient calculation for multidimensional arrays."""
        # Create a 4D field (time, level, lat, lon)
        time, level, lat, lon = 3, 5, 10, 15
        field = np.random.randn(time, level, lat, lon)
        coordinates = np.linspace(0, 360, lon)

        gradient = calculate_gradient(field, coordinates, axis=-1)

        assert gradient.shape == field.shape
        assert np.all(np.isfinite(gradient))

    def test_calculate_gradient_latitude_coordinates(self):
        """Test gradient calculation with latitude coordinates."""
        field = np.array([[1, 2, 3], [2, 4, 6]])
        latitudes = np.array([-30, 0, 30])  # Valid latitude range

        gradient = calculate_gradient(field, latitudes, axis=1)

        assert gradient.shape == field.shape
        assert np.all(np.isfinite(gradient))

    def test_calculate_gradient_validation(self):
        """Test input validation for gradient calculation."""
        field = np.random.randn(5, 10)
        coordinates = np.linspace(0, 100, 8)  # Wrong size

        with pytest.raises(ValueError, match="Coordinate array size"):
            calculate_gradient(field, coordinates, axis=1)

    def test_calculate_gradient_edge_cases(self):
        """Test edge cases for gradient calculation."""
        # Test with single point
        field = np.array([5.0])
        coordinates = np.array([0.0])

        gradient = calculate_gradient(field, coordinates, axis=0)
        assert gradient.shape == (1,)
        assert gradient[0] == 0.0  # Gradient of single point should be 0

        # Test with two points
        field = np.array([1.0, 3.0])
        coordinates = np.array([0.0, 2.0])

        gradient = calculate_gradient(field, coordinates, axis=0)
        assert gradient.shape == (2,)
        assert np.allclose(gradient, [1.0, 1.0])  # Should be 1.0 for both points


class TestMeridionalGradient:
    """Test meridional (latitudinal) gradient calculation."""

    def test_calculate_meridional_gradient_basic(self):
        """Test basic meridional gradient calculation."""
        # Create a field that varies linearly with latitude
        lats = np.linspace(-90, 90, 19)  # 10-degree intervals
        field = lats.reshape(1, -1, 1)  # Shape: (1, lat, 1)

        gradient = calculate_meridional_gradient(field, lats, lat_axis=1)

        assert gradient.shape == field.shape

        # For a linear field, gradient should be approximately constant
        # Convert degrees to radians and account for Earth radius
        expected_gradient = 1.0 / (np.pi / 180.0 * EARTH_RADIUS)
        assert np.allclose(gradient[0, 1:-1, 0], expected_gradient, rtol=0.1)

    def test_calculate_meridional_gradient_climate_data(self, sample_climate_data):
        """Test meridional gradient with realistic climate data."""
        temp = sample_climate_data["temperature"]

        # Take a single time slice
        temp_slice = temp.isel(time=0)

        gradient = calculate_meridional_gradient(
            temp_slice.values, temp_slice.lat.values, lat_axis=0
        )

        assert gradient.shape == temp_slice.shape
        assert np.all(np.isfinite(gradient))

        # Temperature gradient should be generally negative (decreasing toward poles)
        # Check that the mean gradient has the right sign in the tropics
        tropical_mask = (temp_slice.lat >= -30) & (temp_slice.lat <= 30)
        tropical_gradient = gradient[tropical_mask, :]
        assert np.mean(tropical_gradient) < 0  # Generally decreasing northward

    def test_calculate_meridional_gradient_different_axes(self):
        """Test meridional gradient calculation with different axis positions."""
        # Create 3D field (time, lat, lon)
        time, lat, lon = 5, 37, 72
        field = np.random.randn(time, lat, lon)
        latitudes = np.linspace(-90, 90, lat)

        # Test with lat_axis=1
        gradient1 = calculate_meridional_gradient(field, latitudes, lat_axis=1)
        assert gradient1.shape == field.shape

        # Test with different field shape (lat, lon, time)
        field_transposed = field.transpose(1, 2, 0)
        gradient2 = calculate_meridional_gradient(
            field_transposed, latitudes, lat_axis=0
        )
        assert gradient2.shape == field_transposed.shape


class TestZonalGradient:
    """Test zonal (longitudinal) gradient calculation."""

    def test_calculate_zonal_gradient_basic(self):
        """Test basic zonal gradient calculation."""
        lats = np.array([0, 30, 60])  # Different latitudes
        lons = np.linspace(0, 360, 73)  # 5-degree intervals

        # Create a field that varies linearly with longitude at each latitude
        field = np.zeros((3, 73))
        for i, lat in enumerate(lats):
            field[i, :] = lons  # Linear variation with longitude

        gradient = calculate_zonal_gradient(field, lons, lats, lon_axis=1, lat_axis=0)

        assert gradient.shape == field.shape

        # Check that gradient varies with latitude (due to cosine factor)
        gradient_equator = np.mean(gradient[0, 10:-10])  # Equator
        gradient_mid = np.mean(gradient[1, 10:-10])  # 30°N
        gradient_high = np.mean(gradient[2, 10:-10])  # 60°N

        # Gradient should decrease as we move away from equator
        assert gradient_equator > gradient_mid > gradient_high

    def test_calculate_zonal_gradient_3d(self):
        """Test zonal gradient for 3D data."""
        time, lat, lon = 10, 19, 36
        field = np.random.randn(time, lat, lon)
        latitudes = np.linspace(-90, 90, lat)
        longitudes = np.linspace(0, 355, lon)

        gradient = calculate_zonal_gradient(
            field, longitudes, latitudes, lon_axis=-1, lat_axis=-2
        )

        assert gradient.shape == field.shape
        assert np.all(np.isfinite(gradient))

    def test_calculate_zonal_gradient_4d(self):
        """Test zonal gradient for 4D data."""
        time, level, lat, lon = 5, 10, 19, 36
        field = np.random.randn(time, level, lat, lon)
        latitudes = np.linspace(-90, 90, lat)
        longitudes = np.linspace(0, 355, lon)

        gradient = calculate_zonal_gradient(
            field, longitudes, latitudes, lon_axis=-1, lat_axis=-2
        )

        assert gradient.shape == field.shape
        assert np.all(np.isfinite(gradient))

    def test_calculate_zonal_gradient_climate_data(self, sample_climate_data):
        """Test zonal gradient with realistic climate data."""
        temp = sample_climate_data["temperature"]

        # Take a single time slice
        temp_slice = temp.isel(time=0)

        gradient = calculate_zonal_gradient(
            temp_slice.values,
            temp_slice.lon.values,
            temp_slice.lat.values,
            lon_axis=1,
            lat_axis=0,
        )

        assert gradient.shape == temp_slice.shape
        assert np.all(np.isfinite(gradient))


class TestVerticalGradient:
    """Test vertical (pressure level) gradient calculation."""

    def test_calculate_vertical_gradient_basic(self):
        """Test basic vertical gradient calculation."""
        # Standard atmospheric pressure levels (decreasing)
        pressure = np.array([1000, 850, 700, 500, 300, 200, 100]) * 100  # Convert to Pa

        # Create a temperature profile (decreasing with height)
        temp_profile = np.array([288, 280, 270, 255, 230, 220, 210])
        field = temp_profile.reshape(1, -1, 1, 1)  # (time, level, lat, lon)

        gradient = calculate_vertical_gradient(field, pressure, pressure_axis=1)

        assert gradient.shape == field.shape

        # Temperature should generally decrease with decreasing pressure (increasing height)
        # So gradient w.r.t. pressure should be positive
        assert np.mean(gradient[0, 1:-1, 0, 0]) > 0

    def test_calculate_vertical_gradient_multidimensional(self):
        """Test vertical gradient for multidimensional data."""
        time, level, lat, lon = 12, 17, 73, 144
        field = np.random.randn(time, level, lat, lon)
        pressure = np.linspace(100000, 10000, level)  # Decreasing pressure

        gradient = calculate_vertical_gradient(field, pressure, pressure_axis=1)

        assert gradient.shape == field.shape
        assert np.all(np.isfinite(gradient))

    def test_calculate_vertical_gradient_climate_data(self, sample_climate_data):
        """Test vertical gradient with realistic pressure levels."""
        temp = sample_climate_data["temperature"]

        # Add a pressure dimension
        pressure_levels = np.array([1000, 850, 700, 500, 300, 200]) * 100  # Pa
        temp_3d = np.tile(
            temp.values[:, np.newaxis, :, :], (1, len(pressure_levels), 1, 1)
        )

        # Add realistic vertical temperature variation
        for i, p in enumerate(pressure_levels):
            temp_3d[:, i, :, :] *= (p / 100000) ** 0.2  # Simple lapse rate

        gradient = calculate_vertical_gradient(
            temp_3d, pressure_levels, pressure_axis=1
        )

        assert gradient.shape == temp_3d.shape
        assert np.all(np.isfinite(gradient))


class TestGradientsIntegration:
    """Integration tests for gradients module."""

    def test_gradients_with_realistic_data(self, sample_climate_data, sample_2d_field):
        """Test all gradient functions with realistic data."""
        temp_2d = sample_2d_field

        # Test meridional gradient
        grad_lat = calculate_meridional_gradient(
            temp_2d.values, temp_2d.lat.values, lat_axis=0
        )

        # Test zonal gradient
        grad_lon = calculate_zonal_gradient(
            temp_2d.values,
            temp_2d.lon.values,
            temp_2d.lat.values,
            lon_axis=1,
            lat_axis=0,
        )

        assert grad_lat.shape == temp_2d.shape
        assert grad_lon.shape == temp_2d.shape
        assert np.all(np.isfinite(grad_lat))
        assert np.all(np.isfinite(grad_lon))

        # Gradients should have reasonable magnitudes
        assert np.abs(np.mean(grad_lat)) < 1e-3  # K/m
        assert np.abs(np.mean(grad_lon)) < 1e-3  # K/m

    def test_gradients_error_handling(self):
        """Test comprehensive error handling."""
        field = np.random.randn(10, 20)
        coordinates = np.linspace(0, 100, 15)  # Wrong size

        # Test dimension mismatch
        with pytest.raises(ValueError):
            calculate_gradient(field, coordinates, axis=1)

        # Test invalid axis
        valid_coordinates = np.linspace(0, 100, 20)
        with pytest.raises(IndexError):
            calculate_gradient(field, valid_coordinates, axis=5)  # Invalid axis

    def test_gradients_consistency(self):
        """Test consistency between different gradient functions."""
        # Create simple field
        lats = np.linspace(-45, 45, 19)
        lons = np.linspace(0, 360, 36)
        field = np.random.randn(19, 36)

        # Calculate meridional gradient using both functions
        grad1 = calculate_gradient(field, lats, axis=0, radius=EARTH_RADIUS)
        grad2 = calculate_meridional_gradient(field, lats, lat_axis=0)

        # Should be approximately equal
        assert np.allclose(grad1, grad2, rtol=1e-10)

    def test_gradients_units_consistency(self):
        """Test that gradient units are physically reasonable."""
        # Temperature field in Kelvin
        temp = 273.15 + 20 * np.random.randn(37, 72)
        lats = np.linspace(-90, 90, 37)

        gradient = calculate_meridional_gradient(temp, lats, lat_axis=0)

        # Gradient units should be K/m
        # Typical atmospheric temperature gradients are O(10^-5) to O(10^-3) K/m
        typical_magnitude = np.percentile(np.abs(gradient), 75)
        assert 1e-6 < typical_magnitude < 1e-2  # Reasonable range


# Performance tests (marked as slow)
@pytest.mark.slow
class TestGradientsPerformance:
    """Performance tests for gradients module."""

    def test_gradient_calculation_large_data(self):
        """Test gradient calculation with large datasets."""
        # Large 4D dataset
        time, level, lat, lon = 100, 50, 180, 360
        field = np.random.randn(time, level, lat, lon)
        coordinates = np.linspace(0, 357.5, lon)

        # Should complete without memory issues
        gradient = calculate_gradient(field, coordinates, axis=-1)

        assert gradient.shape == field.shape
        assert np.all(np.isfinite(gradient))

    def test_zonal_gradient_large_data(self):
        """Test zonal gradient calculation with large datasets."""
        time, level, lat, lon = 50, 30, 90, 180
        field = np.random.randn(time, level, lat, lon)
        latitudes = np.linspace(-90, 90, lat)
        longitudes = np.linspace(0, 357, lon)

        gradient = calculate_zonal_gradient(
            field, longitudes, latitudes, lon_axis=-1, lat_axis=-2
        )

        assert gradient.shape == field.shape
        assert np.all(np.isfinite(gradient))


class TestGradientsEdgeCases:
    """Test edge cases and special conditions."""

    def test_gradient_uniform_field(self):
        """Test gradient of uniform field (should be zero)."""
        field = np.ones((10, 15))
        coordinates = np.linspace(0, 100, 15)

        gradient = calculate_gradient(field, coordinates, axis=1)

        # Gradient of uniform field should be approximately zero
        assert np.allclose(gradient, 0, atol=1e-15)

    def test_gradient_linear_field(self):
        """Test gradient of linear field (should be constant)."""
        coordinates = np.linspace(0, 100, 21)
        field = coordinates.reshape(1, -1)  # Linear field

        gradient = calculate_gradient(field, coordinates, axis=1)

        # Gradient should be approximately constant (and equal to 1)
        expected_gradient = 1.0
        assert np.allclose(gradient[0, 1:-1], expected_gradient, rtol=0.01)

    def test_gradient_with_cyclic_coordinates(self):
        """Test gradient calculation with longitude coordinates (cyclic)."""
        # Longitude coordinates from 0 to 355 degrees
        lons = np.linspace(0, 355, 72)
        lats = np.array([0])  # Equator only

        # Create a field that should have continuous gradient across 0/360 boundary
        field = np.sin(np.radians(lons * 2)).reshape(1, -1)

        gradient = calculate_zonal_gradient(field, lons, lats, lon_axis=1, lat_axis=0)

        assert gradient.shape == field.shape
        assert np.all(np.isfinite(gradient))

        # The gradient should be approximately 2*cos(2*theta) scaled by Earth radius
        expected_magnitude = 2.0 / (np.pi / 180.0 * EARTH_RADIUS)
        assert np.abs(np.mean(gradient[0, 10:-10])) < expected_magnitude * 2


if __name__ == "__main__":
    # Quick test runner
    pytest.main([__file__, "-v"])
