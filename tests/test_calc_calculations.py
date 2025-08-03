"""
Tests for skyborn.calc.calculations module.

This module tests the statistical and mathematical calculation functions
in the skyborn.calc.calculations module.
"""

import pytest
import numpy as np
import xarray as xr
from numpy.testing import assert_array_almost_equal, assert_array_equal
from skyborn.calc.calculations import linear_regression


class TestLinearRegression:
    """Test linear regression functionality."""

    def test_linear_regression_numpy_arrays(self, sample_regression_data):
        """Test linear regression with numpy arrays."""
        data, predictor = sample_regression_data

        # Perform regression
        slopes, p_values = linear_regression(data, predictor)

        # Check output shapes
        assert slopes.shape == data.shape[1:]
        assert p_values.shape == data.shape[1:]

        # Check that outputs are finite
        assert np.all(np.isfinite(slopes))
        assert np.all(np.isfinite(p_values))

        # Check p-values are in valid range [0, 1]
        assert np.all(p_values >= 0)
        assert np.all(p_values <= 1)

    def test_linear_regression_xarray(self, sample_regression_data):
        """Test linear regression with xarray DataArrays."""
        data_np, predictor_np = sample_regression_data

        # Convert to xarray
        data_xr = xr.DataArray(
            data_np,
            dims=["time", "lat", "lon"],
            coords={
                "time": np.arange(data_np.shape[0]),
                "lat": np.arange(data_np.shape[1]),
                "lon": np.arange(data_np.shape[2]),
            },
        )
        predictor_xr = xr.DataArray(predictor_np, dims=["time"])

        # Test with xarray inputs
        slopes, p_values = linear_regression(data_xr, predictor_xr)

        # Should produce same results as numpy version
        slopes_np, p_values_np = linear_regression(data_np, predictor_np)

        assert_array_almost_equal(slopes, slopes_np)
        assert_array_almost_equal(p_values, p_values_np)

    def test_linear_regression_known_relationship(self):
        """Test linear regression with known relationship."""
        n_time = 100
        predictor = np.linspace(-2, 2, n_time)

        # Create data with known slope and intercept
        true_slope = 3.5
        true_intercept = 1.2
        noise_level = 0.1

        # Single grid point with known relationship
        data = np.zeros((n_time, 1, 1))
        data[:, 0, 0] = (
            true_slope * predictor
            + true_intercept
            + np.random.randn(n_time) * noise_level
        )

        slopes, p_values = linear_regression(data, predictor)

        # Check that recovered slope is close to true slope
        assert abs(slopes[0, 0] - true_slope) < 0.2

        # With strong relationship, p-value should be very small
        assert p_values[0, 0] < 0.01

    def test_linear_regression_no_relationship(self):
        """Test linear regression with no relationship (random data)."""
        n_time = 50
        predictor = np.random.randn(n_time)

        # Create random data with no relationship to predictor
        data = np.random.randn(n_time, 5, 5)

        slopes, p_values = linear_regression(data, predictor)

        # Slopes should be close to zero on average
        assert abs(np.mean(slopes)) < 0.5

        # Most p-values should be > 0.05 (not significant)
        assert np.mean(p_values > 0.05) > 0.8

    def test_linear_regression_input_validation(self):
        """Test input validation for linear regression."""
        # Test mismatched dimensions
        data = np.random.randn(50, 10, 10)
        predictor = np.random.randn(40)  # Wrong length

        with pytest.raises(ValueError, match="Number of samples in data"):
            linear_regression(data, predictor)

        # Test with 2D data (should fail)
        data_2d = np.random.randn(50, 10)
        predictor_valid = np.random.randn(50)

        with pytest.raises(ValueError):
            linear_regression(data_2d, predictor_valid)

    def test_linear_regression_edge_cases(self):
        """Test edge cases for linear regression."""
        # Test with constant predictor
        n_time = 30
        predictor = np.ones(n_time)  # Constant predictor
        data = np.random.randn(n_time, 3, 3)

        slopes, p_values = linear_regression(data, predictor)

        # With constant predictor, slopes should be near zero
        assert np.all(np.abs(slopes) < 1e-10)

        # Test with single time step
        predictor_single = np.array([1.0])
        data_single = np.random.randn(1, 2, 2)

        # This should work but produce NaN p-values
        slopes, p_values = linear_regression(data_single, predictor_single)
        assert slopes.shape == (2, 2)
        # With only one point, can't compute meaningful statistics

    def test_linear_regression_output_types(self, sample_regression_data):
        """Test that outputs are numpy arrays regardless of input type."""
        data, predictor = sample_regression_data

        slopes, p_values = linear_regression(data, predictor)

        assert isinstance(slopes, np.ndarray)
        assert isinstance(p_values, np.ndarray)

        # Test with xarray input
        data_xr = xr.DataArray(data, dims=["time", "lat", "lon"])
        predictor_xr = xr.DataArray(predictor, dims=["time"])

        slopes_xr, p_values_xr = linear_regression(data_xr, predictor_xr)

        assert isinstance(slopes_xr, np.ndarray)
        assert isinstance(p_values_xr, np.ndarray)


class TestCalculationsIntegration:
    """Integration tests for calculations module."""

    def test_calculations_with_climate_data(self, sample_climate_data):
        """Test calculations using realistic climate data."""
        temp = sample_climate_data["temperature"]

        # Create a simple index (e.g., global mean temperature)
        global_temp = temp.mean(dim=["lat", "lon"])

        # Test regression of local temperature against global temperature
        slopes, p_values = linear_regression(temp.values, global_temp.values)

        # Should get reasonable results
        assert slopes.shape == (73, 144)  # lat, lon
        assert p_values.shape == (73, 144)

        # Most locations should have positive correlation with global mean
        assert np.mean(slopes > 0) > 0.7

        # Many locations should have significant correlations
        assert np.mean(p_values < 0.05) > 0.3

    def test_calculations_error_handling(self):
        """Test comprehensive error handling."""
        # Test with wrong input types
        with pytest.raises(ValueError):
            linear_regression("not_an_array", np.array([1, 2, 3]))

        with pytest.raises(ValueError):
            linear_regression(np.array([1, 2, 3]), "not_an_array")

        # Test with incompatible shapes
        data = np.random.randn(10, 5, 5)
        predictor = np.random.randn(5, 5)  # Wrong shape

        with pytest.raises(ValueError):
            linear_regression(data, predictor)


# Performance tests (marked as slow)
@pytest.mark.slow
class TestCalculationsPerformance:
    """Performance tests for calculations module."""

    def test_linear_regression_large_data(self):
        """Test linear regression with large datasets."""
        # Large dataset
        data = np.random.randn(1000, 100, 100)
        predictor = np.random.randn(1000)

        # Should complete without memory issues
        slopes, p_values = linear_regression(data, predictor)

        assert slopes.shape == (100, 100)
        assert p_values.shape == (100, 100)
        assert np.all(np.isfinite(slopes))
        assert np.all(np.isfinite(p_values))


if __name__ == "__main__":
    # Quick test runner
    pytest.main([__file__, "-v"])
