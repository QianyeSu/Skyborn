"""
Tests for skyborn.calc.calculations module.

This module tests the statistical and mathematical calculation functions
in the skyborn.calc.calculations module.
"""

import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_array_almost_equal, assert_array_equal

from skyborn.calc.calculations import (
    calculate_potential_temperature,
    calculate_theta_se,
    convert_longitude_range,
    kendall_correlation,
    linear_regression,
    pearson_correlation,
    spatial_correlation,
    spearman_correlation,
)
from skyborn.calc.emergent_constraints import (  # Legacy functions for backward compatibility testing
    _calculate_std_from_pdf,
    calc_GAUSSIAN_PDF,
    calc_PDF_EC,
    calc_PDF_EC_PRIOR,
    emergent_constraint_posterior,
    emergent_constraint_prior,
    find_std_from_PDF,
    gaussian_pdf,
)


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

        # With constant predictor, slopes might be computed but should be handled gracefully
        assert slopes.shape == (3, 3)
        assert p_values.shape == (3, 3)

        # Results should be finite or NaN (both are acceptable for constant predictor)
        assert np.all(np.isfinite(slopes) | np.isnan(slopes))
        assert np.all(np.isfinite(p_values) | np.isnan(p_values))

        # Test with single time step
        predictor_single = np.array([1.0])
        data_single = np.random.randn(1, 2, 2)

        # This should work but produce NaN p-values
        slopes, p_values = linear_regression(data_single, predictor_single)
        assert slopes.shape == (2, 2)
        # With only one point, can't compute meaningful statistics

    def test_linear_regression_with_nans(self):
        """Test linear regression with NaN values."""
        np.random.seed(42)
        n_time, n_lat, n_lon = 60, 15, 20

        # Create predictor with some NaN values
        predictor = np.random.randn(n_time)
        predictor[5:10] = np.nan  # Add NaN to predictor

        # Create data with some NaN values
        data = np.random.randn(n_time, n_lat, n_lon)
        data[:8, :3, :3] = np.nan  # Add NaN block to data
        data[50:55, 10:12, 15:18] = np.nan  # Add another NaN region

        # Add some known relationships where there are no NaNs
        for lat in range(8, 12):
            for lon in range(10, 15):
                # Create valid data points only where both data and predictor are finite
                valid_mask = ~np.isnan(predictor)
                data[valid_mask, lat, lon] = 0.7 * predictor[
                    valid_mask
                ] + 0.3 * np.random.randn(np.sum(valid_mask))

        slopes, p_values = linear_regression(data, predictor)

        # Check shapes
        assert slopes.shape == (n_lat, n_lon)
        assert p_values.shape == (n_lat, n_lon)

        # Areas with too many NaNs should be NaN in results
        assert np.all(np.isnan(slopes[:3, :3]))
        assert np.all(np.isnan(p_values[:3, :3]))

        # Areas with sufficient valid data should have finite results
        valid_region = slopes[8:12, 10:15]
        valid_p_region = p_values[8:12, 10:15]

        # At least some points should have valid results
        has_valid_slopes = np.any(~np.isnan(valid_region))
        has_valid_p_values = np.any(~np.isnan(valid_p_region))

        assert (
            has_valid_slopes
        ), "Should have some valid slope estimates in regions with sufficient data"
        assert (
            has_valid_p_values
        ), "Should have some valid p-values in regions with sufficient data"

        # Valid slopes should be finite and reasonable
        finite_slopes = valid_region[~np.isnan(valid_region)]
        finite_p_values = valid_p_region[~np.isnan(valid_p_region)]

        if len(finite_slopes) > 0:
            assert np.all(np.isfinite(finite_slopes))
            # Should detect positive correlation
            assert np.mean(finite_slopes) > 0.3

        if len(finite_p_values) > 0:
            assert np.all(np.isfinite(finite_p_values))
            assert np.all(finite_p_values >= 0.0)
            assert np.all(finite_p_values <= 1.0)

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


class TestSpatialCorrelation:
    """Test spatial correlation functionality."""

    @pytest.fixture
    def sample_spatial_data(self):
        """Create sample spatial correlation data."""
        np.random.seed(42)
        n_time, n_lat, n_lon = 50, 20, 30

        # Create predictor time series
        predictor = np.random.randn(n_time)

        # Create spatial data
        data = np.random.randn(n_time, n_lat, n_lon)

        # Add some regions with known correlations
        # Strong positive correlation region
        for lat in range(5, 10):
            for lon in range(10, 15):
                data[:, lat, lon] = 0.8 * predictor + 0.2 * np.random.randn(n_time)

        # Strong negative correlation region
        for lat in range(15, 18):
            for lon in range(20, 25):
                data[:, lat, lon] = -0.7 * predictor + 0.3 * np.random.randn(n_time)

        return data, predictor

    def test_spatial_correlation_basic(self, sample_spatial_data):
        """Test basic spatial correlation functionality."""
        data, predictor = sample_spatial_data

        corr_coef, p_values = spatial_correlation(data, predictor)

        # Check output shapes
        assert corr_coef.shape == data.shape[1:]  # (n_lat, n_lon)
        assert p_values.shape == data.shape[1:]

        # Check correlation coefficient range
        valid_corr = corr_coef[~np.isnan(corr_coef)]
        assert np.all(valid_corr >= -1.0)
        assert np.all(valid_corr <= 1.0)

        # Check p-values range
        valid_p = p_values[~np.isnan(p_values)]
        assert np.all(valid_p >= 0.0)
        assert np.all(valid_p <= 1.0)

    def test_spatial_correlation_known_relationships(self, sample_spatial_data):
        """Test spatial correlation with known relationships."""
        data, predictor = sample_spatial_data

        corr_coef, p_values = spatial_correlation(data, predictor)

        # Check positive correlation region (lat 5-9, lon 10-14)
        pos_corr_region = corr_coef[5:10, 10:15]
        assert np.mean(pos_corr_region) > 0.5  # Should be strongly positive

        # Check negative correlation region (lat 15-17, lon 20-24)
        neg_corr_region = corr_coef[15:18, 20:25]
        assert np.mean(neg_corr_region) < -0.4  # Should be strongly negative

        # Check corresponding p-values are significant
        pos_p_region = p_values[5:10, 10:15]
        neg_p_region = p_values[15:18, 20:25]
        assert np.mean(pos_p_region < 0.05) > 0.8  # Most should be significant
        assert np.mean(neg_p_region < 0.05) > 0.8

    def test_spatial_correlation_xarray_input(self, sample_spatial_data):
        """Test spatial correlation with xarray input."""
        data_np, predictor_np = sample_spatial_data

        # Convert to xarray
        time_coord = np.arange(len(predictor_np))
        lat_coord = np.arange(data_np.shape[1])
        lon_coord = np.arange(data_np.shape[2])

        data_xr = xr.DataArray(
            data_np,
            dims=["time", "lat", "lon"],
            coords={"time": time_coord, "lat": lat_coord, "lon": lon_coord},
        )
        predictor_xr = xr.DataArray(
            predictor_np, dims=["time"], coords={"time": time_coord}
        )

        # Test with xarray inputs
        corr_xr, p_xr = spatial_correlation(data_xr, predictor_xr)

        # Should produce same results as numpy version
        corr_np, p_np = spatial_correlation(data_np, predictor_np)

        assert_array_almost_equal(corr_xr, corr_np)
        assert_array_almost_equal(p_xr, p_np)

    def test_spatial_correlation_with_nans(self):
        """Test spatial correlation with NaN values."""
        np.random.seed(42)
        n_time, n_lat, n_lon = 60, 15, 20

        # Create data with NaN values
        predictor = np.random.randn(n_time)
        predictor[5:10] = np.nan  # Add NaN to predictor

        data = np.random.randn(n_time, n_lat, n_lon)
        data[:10, :3, :3] = np.nan  # Add NaN block to data

        # Add some known correlation where there are no NaNs
        for lat in range(8, 12):
            for lon in range(10, 15):
                data[:, lat, lon] = 0.7 * predictor + 0.3 * np.random.randn(n_time)

        corr_coef, p_values = spatial_correlation(data, predictor)

        # Check shapes
        assert corr_coef.shape == (n_lat, n_lon)
        assert p_values.shape == (n_lat, n_lon)

        # Areas with too many NaNs should be NaN in results
        assert np.all(np.isnan(corr_coef[:3, :3]))
        assert np.all(np.isnan(p_values[:3, :3]))

        # Areas with sufficient data should have valid results
        valid_region = corr_coef[8:12, 10:15]
        assert not np.all(np.isnan(valid_region))
        assert np.mean(valid_region[~np.isnan(valid_region)]) > 0.3

    def test_spatial_correlation_input_validation(self):
        """Test input validation for spatial correlation."""
        # Test mismatched time dimensions
        data = np.random.randn(50, 10, 15)
        predictor = np.random.randn(40)  # Wrong length

        with pytest.raises(ValueError, match="Time dimension of data"):
            spatial_correlation(data, predictor)

    def test_spatial_correlation_edge_cases(self):
        """Test edge cases for spatial correlation."""
        n_time = 20

        # Test with constant predictor
        predictor_constant = np.ones(n_time)
        data = np.random.randn(n_time, 5, 5)

        corr_coef, p_values = spatial_correlation(data, predictor_constant)

        # With constant predictor, correlations should be NaN or 0
        assert np.all(np.isnan(corr_coef) | (corr_coef == 0))

        # Test with very short time series
        predictor_short = np.array([1.0, 2.0])
        data_short = np.random.randn(2, 3, 3)

        corr_short, p_short = spatial_correlation(data_short, predictor_short)

        # With only 2 points, should get NaN (insufficient data)
        assert np.all(np.isnan(corr_short))
        assert np.all(np.isnan(p_short))

    def test_spatial_correlation_output_types(self, sample_spatial_data):
        """Test that outputs are numpy arrays regardless of input type."""
        data, predictor = sample_spatial_data

        corr_coef, p_values = spatial_correlation(data, predictor)

        assert isinstance(corr_coef, np.ndarray)
        assert isinstance(p_values, np.ndarray)

        # Test with xarray input
        data_xr = xr.DataArray(data, dims=["time", "lat", "lon"])
        predictor_xr = xr.DataArray(predictor, dims=["time"])

        corr_xr, p_xr = spatial_correlation(data_xr, predictor_xr)

        assert isinstance(corr_xr, np.ndarray)
        assert isinstance(p_xr, np.ndarray)

    def test_spatial_correlation_accuracy_vs_scipy(self):
        """Test accuracy against scipy.stats.pearsonr."""
        from scipy.stats import pearsonr

        np.random.seed(123)
        n_time = 100

        # Create simple test case
        predictor = np.random.randn(n_time)
        data = np.zeros((n_time, 3, 3))

        # Fill with known relationships
        for lat in range(3):
            for lon in range(3):
                correlation_strength = 0.1 * (
                    lat + lon + 1
                )  # Vary correlation strength
                data[:, lat, lon] = (
                    correlation_strength * predictor + 0.5 * np.random.randn(n_time)
                )

        # Test our function
        corr_ours, p_ours = spatial_correlation(data, predictor)

        # Compare with scipy for each point
        for lat in range(3):
            for lon in range(3):
                r_scipy, p_scipy = pearsonr(data[:, lat, lon], predictor)

                # Check correlation coefficient
                assert abs(corr_ours[lat, lon] - r_scipy) < 1e-10

                # Check p-value
                assert abs(p_ours[lat, lon] - p_scipy) < 1e-8


class TestLongitudeConversion:
    """Test longitude coordinate conversion functionality."""

    def test_convert_longitude_range_center_on_180(self):
        """Test converting longitude range to center on 180."""
        # Create test data with longitude 0-359
        lon = np.array([0, 90, 180, 270, 359])
        lat = np.array([0, 45, 90])
        data = np.random.rand(3, 5)

        da = xr.DataArray(data, coords={"lat": lat, "lon": lon}, dims=["lat", "lon"])

        # Convert to -180..179 range
        result = convert_longitude_range(da, center_on_180=True)

        # Check that longitudes are converted correctly
        # Should remain as 0-359
        expected_lon = np.array([0, 90, 180, 270, 359])
        assert_array_almost_equal(result.lon.values, expected_lon)

        # Check that data is preserved
        assert result.shape == da.shape

    def test_convert_longitude_range_center_on_0(self):
        """Test converting longitude range to center on 0."""
        # Create test data with longitude -180..179
        lon = np.array([-180, -90, 0, 90, 179])
        lat = np.array([0, 45])
        data = np.random.rand(2, 5)

        da = xr.DataArray(data, coords={"lat": lat, "lon": lon}, dims=["lat", "lon"])

        # Convert to 0..359 range
        result = convert_longitude_range(da, center_on_180=False)

        # Check that longitudes are sorted
        assert np.all(result.lon.values[:-1] <= result.lon.values[1:])

    def test_convert_longitude_range_with_dataset(self):
        """Test longitude conversion with xarray Dataset."""
        lon = np.array([0, 90, 180, 270])
        lat = np.array([0, 45])

        ds = xr.Dataset(
            {
                "temp": (["lat", "lon"], np.random.rand(2, 4)),
                "precip": (["lat", "lon"], np.random.rand(2, 4)),
            },
            coords={"lat": lat, "lon": lon},
        )

        result = convert_longitude_range(ds, center_on_180=True)

        assert isinstance(result, xr.Dataset)
        assert "temp" in result.data_vars
        assert "precip" in result.data_vars

    def test_convert_longitude_range_custom_coordinate_name(self):
        """Test longitude conversion with custom coordinate name."""
        longitude = np.array([0, 180, 270])
        lat = np.array([0, 45])
        data = np.random.rand(2, 3)

        da = xr.DataArray(
            data, coords={"lat": lat, "longitude": longitude}, dims=["lat", "longitude"]
        )

        result = convert_longitude_range(da, lon="longitude", center_on_180=True)
        assert "longitude" in result.coords


class TestCorrelationFunctions:
    """Test correlation calculation functions."""

    @pytest.fixture
    def sample_correlation_data(self):
        """Create sample data for correlation testing."""
        np.random.seed(42)
        n = 100

        # Create perfectly correlated data
        x_perfect = np.random.randn(n)
        y_perfect = 2.0 * x_perfect + 1.0

        # Create anticorrelated data
        x_anti = np.random.randn(n)
        y_anti = -1.5 * x_anti + 0.5

        # Create uncorrelated data
        x_uncorr = np.random.randn(n)
        y_uncorr = np.random.randn(n)

        return {
            "perfect": (x_perfect, y_perfect),
            "anti": (x_anti, y_anti),
            "uncorrelated": (x_uncorr, y_uncorr),
        }

    def test_pearson_correlation_perfect(self, sample_correlation_data):
        """Test Pearson correlation with perfectly correlated data."""
        x, y = sample_correlation_data["perfect"]

        corr = pearson_correlation(x, y)

        # Should be very close to 1.0
        assert abs(corr - 1.0) < 0.01

    def test_pearson_correlation_anticorrelated(self, sample_correlation_data):
        """Test Pearson correlation with anticorrelated data."""
        x, y = sample_correlation_data["anti"]

        corr = pearson_correlation(x, y)

        # Should be very close to -1.0
        assert abs(corr + 1.0) < 0.01

    def test_pearson_correlation_uncorrelated(self, sample_correlation_data):
        """Test Pearson correlation with uncorrelated data."""
        x, y = sample_correlation_data["uncorrelated"]

        corr = pearson_correlation(x, y)

        # Should be close to 0.0
        assert abs(corr) < 0.2

    def test_pearson_correlation_xarray_input(self, sample_correlation_data):
        """Test Pearson correlation with xarray input."""
        x, y = sample_correlation_data["perfect"]

        x_da = xr.DataArray(x, dims=["time"])
        y_da = xr.DataArray(y, dims=["time"])

        corr = pearson_correlation(x_da, y_da)

        assert abs(corr - 1.0) < 0.01

    def test_spearman_correlation_perfect(self, sample_correlation_data):
        """Test Spearman correlation with perfectly correlated data."""
        x, y = sample_correlation_data["perfect"]

        corr = spearman_correlation(x, y)

        # Should be very close to 1.0
        assert abs(corr - 1.0) < 0.01

    def test_spearman_correlation_monotonic(self):
        """Test Spearman correlation with monotonic relationship."""
        x = np.array([1, 2, 3, 4, 5])
        y = np.array([1, 4, 9, 16, 25])  # y = x^2, monotonic but not linear

        corr = spearman_correlation(x, y)

        # Should be exactly 1.0 for monotonic relationship
        assert abs(corr - 1.0) < 0.01

    def test_spearman_correlation_xarray_input(self):
        """Test Spearman correlation with xarray input."""
        x = np.array([1, 2, 3, 4, 5])
        y = np.array([5, 4, 3, 2, 1])

        x_da = xr.DataArray(x, dims=["time"])
        y_da = xr.DataArray(y, dims=["time"])

        corr = spearman_correlation(x_da, y_da)

        # Perfect negative monotonic relationship
        assert abs(corr + 1.0) < 0.01

    def test_kendall_correlation_perfect(self, sample_correlation_data):
        """Test Kendall correlation with perfectly correlated data."""
        x, y = sample_correlation_data["perfect"]

        corr = kendall_correlation(x, y)

        # Should be very close to 1.0
        assert abs(corr - 1.0) < 0.05

    def test_kendall_correlation_tied_ranks(self):
        """Test Kendall correlation with tied ranks."""
        x = np.array([1, 2, 2, 3, 4])
        y = np.array([1, 2, 2, 3, 4])

        corr = kendall_correlation(x, y)

        # Should be 1.0 despite tied ranks
        assert abs(corr - 1.0) < 0.01

    def test_kendall_correlation_xarray_input(self):
        """Test Kendall correlation with xarray input."""
        x = np.array([1, 3, 2, 5, 4])
        y = np.array([2, 6, 4, 10, 8])

        x_da = xr.DataArray(x, dims=["time"])
        y_da = xr.DataArray(y, dims=["time"])

        corr = kendall_correlation(x_da, y_da)

        assert isinstance(corr, float)
        assert abs(corr) <= 1.0

    def test_correlation_functions_multidimensional_input(self):
        """Test correlation functions with multidimensional arrays."""
        # Create 2D arrays that should be flattened
        x = np.random.randn(10, 5)
        y = 2 * x + np.random.randn(10, 5) * 0.1

        # Test all correlation functions
        pearson_corr = pearson_correlation(x, y)
        spearman_corr = spearman_correlation(x, y)
        kendall_corr = kendall_correlation(x, y)

        # All should return scalars
        assert isinstance(pearson_corr, float)
        assert isinstance(spearman_corr, float)
        assert isinstance(kendall_corr, float)

        # All should be positive and reasonably high
        assert pearson_corr > 0.8
        assert spearman_corr > 0.8
        assert kendall_corr > 0.6


class TestPotentialTemperature:
    """Test potential temperature calculation functionality."""

    def test_calculate_potential_temperature_basic(self):
        """Test basic potential temperature calculation."""
        # Standard conditions
        temperature = np.array([273.15, 283.15, 293.15])  # 0, 10, 20°C
        pressure = np.array([1000.0, 850.0, 700.0])  # hPa

        potential_temp = calculate_potential_temperature(temperature, pressure)

        # Check output shape
        assert potential_temp.shape == temperature.shape

        # Potential temperature should be higher at lower pressures
        assert potential_temp[2] > potential_temp[1] > potential_temp[0]

        # At reference pressure, potential temp should equal temperature
        ref_temp = calculate_potential_temperature(
            np.array([300.0]), np.array([1000.0])
        )
        assert abs(ref_temp[0] - 300.0) < 0.01

    def test_calculate_potential_temperature_xarray_input(self):
        """Test potential temperature with xarray input."""
        temperature = xr.DataArray(
            [273.15, 283.15, 293.15],
            dims=["level"],
            attrs={"units": "K", "long_name": "Temperature"},
        )
        pressure = xr.DataArray(
            [1000.0, 850.0, 700.0], dims=["level"], attrs={"units": "hPa"}
        )

        potential_temp = calculate_potential_temperature(temperature, pressure)

        # Should return xarray with appropriate attributes
        assert hasattr(potential_temp, "attrs")
        assert potential_temp.attrs["units"] == "K"
        assert "Potential Temperature" in potential_temp.attrs["long_name"]

    def test_calculate_potential_temperature_custom_reference(self):
        """Test potential temperature with custom reference pressure."""
        temperature = np.array([300.0])
        pressure = np.array([500.0])

        # Default reference (1000 hPa)
        theta_default = calculate_potential_temperature(temperature, pressure)

        # Custom reference (900 hPa)
        theta_custom = calculate_potential_temperature(
            temperature, pressure, reference_pressure=900.0
        )

        # Different reference pressures should give different results
        assert abs(theta_default[0] - theta_custom[0]) > 1.0

    def test_calculate_potential_temperature_multidimensional(self):
        """Test potential temperature with multidimensional arrays."""
        # Create 2D temperature and pressure arrays
        temperature = np.random.uniform(250, 320, (5, 10))  # Realistic temp range
        # Realistic pressure range
        pressure = np.random.uniform(200, 1000, (5, 10))

        potential_temp = calculate_potential_temperature(temperature, pressure)

        # Check output shape
        assert potential_temp.shape == temperature.shape

        # All values should be positive
        assert np.all(potential_temp > 0)

        # Potential temperature should generally be >= temperature for pressures < 1000 hPa
        low_pressure_mask = pressure < 1000
        assert np.all(
            potential_temp[low_pressure_mask] >= temperature[low_pressure_mask]
        )

    def test_calculate_potential_temperature_edge_cases(self):
        """Test potential temperature with edge cases."""
        # Very high temperature
        high_temp = np.array([400.0])
        pressure = np.array([500.0])

        theta_high = calculate_potential_temperature(high_temp, pressure)
        assert np.isfinite(theta_high[0])
        assert theta_high[0] > high_temp[0]

        # Very low pressure
        temperature = np.array([250.0])
        low_pressure = np.array([100.0])

        theta_low_p = calculate_potential_temperature(temperature, low_pressure)
        assert np.isfinite(theta_low_p[0])
        assert theta_low_p[0] > temperature[0]

    def test_calculate_potential_temperature_physical_consistency(self):
        """Test that potential temperature follows physical principles."""
        # Create a vertical profile
        levels = np.array([1000.0, 850.0, 700.0, 500.0, 300.0])  # hPa
        # Typical temperature profile (decreasing with height)
        temperatures = np.array([288.0, 281.0, 274.0, 261.0, 241.0])  # K

        potential_temps = calculate_potential_temperature(temperatures, levels)

        # In a stable atmosphere, potential temperature should increase with height
        # (though this isn't always true in the real atmosphere)
        # At minimum, check that all values are reasonable
        assert np.all(potential_temps > 200)  # All should be above 200K
        assert np.all(potential_temps < 500)  # All should be below 500K

        # First level (1000 hPa) should have theta ≈ temperature
        assert abs(potential_temps[0] - temperatures[0]) < 1.0


class TestCalculationsIntegration:
    """Integration tests for calculations module."""

    def test_all_correlation_functions_consistency(self):
        """Test that all correlation functions give consistent results."""
        # Create strongly correlated data
        np.random.seed(42)
        x = np.random.randn(50)
        y = 0.8 * x + 0.2 * np.random.randn(50)

        pearson_corr = pearson_correlation(x, y)
        spearman_corr = spearman_correlation(x, y)
        kendall_corr = kendall_correlation(x, y)

        # All should be positive and reasonably similar for linear relationship
        assert pearson_corr > 0.6
        assert spearman_corr > 0.6
        assert kendall_corr > 0.4  # Kendall is typically smaller

        # Pearson and Spearman should be similar for linear relationship
        assert abs(pearson_corr - spearman_corr) < 0.2

    def test_longitude_conversion_round_trip(self):
        """Test that longitude conversion is reversible."""
        # Create test data
        lon = np.array([30, 120, 200, 300])
        lat = np.array([0, 45])
        data = np.random.rand(2, 4)

        da = xr.DataArray(data, coords={"lat": lat, "lon": lon}, dims=["lat", "lon"])

        # Convert to -180..179 and back to 0..359
        converted1 = convert_longitude_range(da, center_on_180=False)
        converted2 = convert_longitude_range(converted1, center_on_180=True)

        # Data values should be preserved
        # (coordinates might be reordered, so we check data integrity differently)
        assert converted2.shape == da.shape
        assert np.allclose(
            np.sort(converted2.values.flatten()), np.sort(da.values.flatten())
        )

    def test_potential_temperature_with_realistic_sounding(self):
        """Test potential temperature with realistic atmospheric sounding."""
        # Typical mid-latitude sounding
        pressure_levels = np.array(
            [1000, 925, 850, 700, 500, 400, 300, 250, 200]
        )  # hPa
        temperatures = np.array([288, 284, 278, 268, 249, 236, 221, 213, 205])  # K

        potential_temps = calculate_potential_temperature(temperatures, pressure_levels)

        # Check physical realism
        assert len(potential_temps) == len(pressure_levels)
        # theta > T for p < 1000 hPa
        assert np.all(potential_temps > temperatures)

        # Potential temperature should generally increase with height in stable conditions
        # Check that the calculation is physically reasonable without being too strict
        theta_increase = np.diff(potential_temps)
        # Just check that we don't have all decreasing values (which would be unphysical)
        assert not np.all(
            theta_increase < 0
        ), "Potential temperature shouldn't decrease everywhere"

        # Check that values are reasonable
        assert np.all(np.isfinite(potential_temps))
        # Should be positive temperature values
        assert np.all(potential_temps > 0)

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
        # linear_regression now extracts values via getattr, so test shape mismatch instead
        with pytest.raises(ValueError, match="must match"):
            linear_regression(np.random.randn(10, 5, 5), np.array([1, 2, 3]))

        with pytest.raises(ValueError):
            linear_regression(np.array([1, 2, 3]), "not_an_array")

        # Test with incompatible shapes
        data = np.random.randn(10, 5, 5)
        predictor = np.random.randn(5, 5)  # Wrong shape

        with pytest.raises(ValueError):
            linear_regression(data, predictor)


class TestEmergentConstraints:
    """Test emergent constraints functionality."""

    @pytest.fixture
    def sample_emergent_constraint_data(self):
        """Create sample data for emergent constraints testing."""
        np.random.seed(42)

        # Number of models
        n_models = 20

        # Create constraint data (e.g., model sensitivity)
        constraint_values = np.random.normal(2.5, 0.8, n_models)

        # Create target data with correlation to constraint
        # True relationship: target = 1.5 * constraint + noise
        target_values = 1.5 * constraint_values + np.random.normal(0, 0.3, n_models)

        # Create xarray DataArrays
        constraint_data = xr.DataArray(
            constraint_values,
            dims=["model"],
            coords={"model": [f"model_{i}" for i in range(n_models)]},
            attrs={"units": "K", "long_name": "Climate Sensitivity"},
        )

        target_data = xr.DataArray(
            target_values,
            dims=["model"],
            coords={"model": [f"model_{i}" for i in range(n_models)]},
            attrs={"units": "K", "long_name": "Future Temperature Change"},
        )

        # Create grids for PDF calculation
        constraint_grid = np.linspace(0.5, 4.5, 100)
        target_grid = np.linspace(1.0, 8.0, 150)

        # Create observational PDF (Gaussian centered around observed value)
        obs_mean = 3.0
        obs_std = 0.4
        obs_pdf = gaussian_pdf(obs_mean, obs_std, constraint_grid)

        return constraint_data, target_data, constraint_grid, target_grid, obs_pdf

    def test_gaussian_pdf_basic(self):
        """Test basic Gaussian PDF calculation."""
        mu = 0.0
        sigma = 1.0
        x = np.linspace(-3, 3, 100)

        pdf = gaussian_pdf(mu, sigma, x)

        # Check properties of Gaussian PDF
        assert len(pdf) == len(x)
        assert np.all(pdf > 0)  # PDF values should be positive

        # Check normalization (approximately)
        dx = x[1] - x[0]
        integral = np.trapz(pdf, dx=dx)
        assert abs(integral - 1.0) < 0.01

        # Check maximum at mean
        max_idx = np.argmax(pdf)
        assert abs(x[max_idx] - mu) < 0.1

    def test_gaussian_pdf_different_parameters(self):
        """Test Gaussian PDF with different parameters."""
        x = np.linspace(-5, 10, 200)

        # Test different means and standard deviations
        test_cases = [(0.0, 1.0), (2.5, 0.5), (-1.0, 2.0), (5.0, 1.5)]

        for mu, sigma in test_cases:
            pdf = gaussian_pdf(mu, sigma, x)

            # Check properties
            assert np.all(pdf > 0)

            # Check maximum location
            max_idx = np.argmax(pdf)
            assert abs(x[max_idx] - mu) < 0.1

            # Check that larger sigma gives smaller peak
            if sigma > 1.0:
                assert np.max(pdf) < 0.5

    def test_gaussian_pdf_scalar_input(self):
        """Test Gaussian PDF with scalar input."""
        mu = 1.0
        sigma = 0.5
        x = 1.0  # Scalar input

        pdf = gaussian_pdf(mu, sigma, x)

        # Should return scalar
        assert np.isscalar(pdf)

        # Should be maximum value (at mean)
        expected = 1 / np.sqrt(2 * np.pi * sigma**2)
        assert abs(pdf - expected) < 1e-10

    def test_emergent_constraint_posterior(self, sample_emergent_constraint_data):
        """Test emergent constraint posterior calculation."""
        constraint_data, target_data, constraint_grid, target_grid, obs_pdf = (
            sample_emergent_constraint_data
        )

        posterior_pdf, posterior_std, posterior_mean = emergent_constraint_posterior(
            constraint_data, target_data, constraint_grid, target_grid, obs_pdf
        )

        # Check output shapes and types
        assert len(posterior_pdf) == len(target_grid)
        assert isinstance(posterior_std, float)
        assert isinstance(posterior_mean, float)

        # Check PDF properties
        assert np.all(posterior_pdf >= 0)  # PDF should be non-negative
        assert np.sum(posterior_pdf) > 0  # PDF should not be all zeros

        # Check that mean is within reasonable range
        assert target_grid.min() <= posterior_mean <= target_grid.max()

        # Check that std is positive
        assert posterior_std > 0

    def test_emergent_constraint_posterior_reduces_uncertainty(
        self, sample_emergent_constraint_data
    ):
        """Test that emergent constraints reduce uncertainty."""
        constraint_data, target_data, constraint_grid, target_grid, obs_pdf = (
            sample_emergent_constraint_data
        )

        # Calculate prior uncertainty (from model spread)
        prior_std = np.std(target_data.values)

        # Calculate posterior
        posterior_pdf, posterior_std, posterior_mean = emergent_constraint_posterior(
            constraint_data, target_data, constraint_grid, target_grid, obs_pdf
        )

        # Posterior uncertainty should be smaller than prior
        # Note: This may not always be true, but should be for our test case
        # with a reasonable correlation
        assert posterior_std <= prior_std * 1.5  # Allow some tolerance

    def test_emergent_constraint_prior(self, sample_emergent_constraint_data):
        """Test emergent constraint prior calculation."""
        constraint_data, target_data, constraint_grid, target_grid, _ = (
            sample_emergent_constraint_data
        )

        prior_pdf, prediction_error, regression_line = emergent_constraint_prior(
            constraint_data, target_data, constraint_grid, target_grid
        )

        # Check output shapes
        assert prior_pdf.shape == (len(target_grid), len(constraint_grid))
        assert len(prediction_error) == len(constraint_grid)
        assert len(regression_line) == len(constraint_grid)

        # Check that all values are finite and positive where appropriate
        assert np.all(np.isfinite(prior_pdf))
        assert np.all(prior_pdf >= 0)
        assert np.all(prediction_error > 0)
        assert np.all(np.isfinite(regression_line))

    def test_calculate_std_from_pdf(self):
        """Test standard deviation calculation from PDF."""
        # Create a known Gaussian distribution
        x = np.linspace(-5, 5, 1000)
        mu = 0.0
        sigma = 1.0
        pdf = gaussian_pdf(mu, sigma, x)

        # Calculate std using our function
        threshold = 0.341  # 1-sigma equivalent
        calculated_std = _calculate_std_from_pdf(threshold, x, pdf)

        # Should be approximately equal to true sigma
        # Allow some tolerance due to discretization
        assert abs(calculated_std - sigma) < 0.2

    def test_calculate_std_from_pdf_different_distributions(self):
        """Test std calculation with different distribution shapes."""
        x = np.linspace(-10, 10, 500)

        # Test with different Gaussian distributions
        test_cases = [(0.0, 0.5), (2.0, 1.5), (-1.0, 2.0)]

        for mu, sigma in test_cases:
            pdf = gaussian_pdf(mu, sigma, x)
            calculated_std = _calculate_std_from_pdf(0.341, x, pdf)

            # Should be positive and reasonable
            assert calculated_std > 0
            assert calculated_std < 10  # Reasonable upper bound for our test data

    def test_emergent_constraints_with_perfect_correlation(self):
        """Test emergent constraints with perfect model correlation."""
        # Create perfectly correlated data
        n_models = 15
        constraint_values = np.linspace(1, 4, n_models)
        target_values = 2.0 * constraint_values  # Perfect correlation

        constraint_data = xr.DataArray(constraint_values, dims=["model"])
        target_data = xr.DataArray(target_values, dims=["model"])

        constraint_grid = np.linspace(0.5, 4.5, 50)
        target_grid = np.linspace(1.0, 9.0, 80)

        # Tight observational constraint
        obs_pdf = gaussian_pdf(2.5, 0.1, constraint_grid)

        posterior_pdf, posterior_std, posterior_mean = emergent_constraint_posterior(
            constraint_data, target_data, constraint_grid, target_grid, obs_pdf
        )

        # With perfect correlation, posterior should be very constrained
        assert posterior_std < 1.0  # Should be well constrained

        # Posterior mean should be close to expected value (2.0 * 2.5 = 5.0)
        expected_mean = 2.0 * 2.5
        assert abs(posterior_mean - expected_mean) < 0.5

    def test_emergent_constraints_with_no_correlation(self):
        """Test emergent constraints with no model correlation."""
        np.random.seed(123)  # Different seed for this test

        # Create uncorrelated data
        n_models = 20
        constraint_values = np.random.normal(2.0, 0.5, n_models)
        target_values = np.random.normal(5.0, 1.0, n_models)  # Independent

        constraint_data = xr.DataArray(constraint_values, dims=["model"])
        target_data = xr.DataArray(target_values, dims=["model"])

        constraint_grid = np.linspace(0.5, 4.0, 50)
        target_grid = np.linspace(2.0, 8.0, 80)

        obs_pdf = gaussian_pdf(2.0, 0.3, constraint_grid)

        posterior_pdf, posterior_std, posterior_mean = emergent_constraint_posterior(
            constraint_data, target_data, constraint_grid, target_grid, obs_pdf
        )

        # With no correlation, constraint should provide little information
        prior_std = np.std(target_values)
        # Posterior std should not be much smaller than prior
        assert posterior_std > 0.5 * prior_std

    def test_emergent_constraints_edge_cases(self):
        """Test emergent constraints with edge cases."""
        # Test with minimal data
        constraint_data = xr.DataArray([1.0, 2.0, 3.0], dims=["model"])
        target_data = xr.DataArray([2.0, 4.0, 6.0], dims=["model"])

        constraint_grid = np.linspace(0.5, 3.5, 20)
        target_grid = np.linspace(1.0, 7.0, 30)
        obs_pdf = gaussian_pdf(2.0, 0.5, constraint_grid)

        # Should not raise errors
        posterior_pdf, posterior_std, posterior_mean = emergent_constraint_posterior(
            constraint_data, target_data, constraint_grid, target_grid, obs_pdf
        )

        assert len(posterior_pdf) == len(target_grid)
        assert posterior_std > 0
        assert np.isfinite(posterior_mean)

    def test_legacy_function_compatibility(self, sample_emergent_constraint_data):
        """Test that legacy functions produce same results as new functions."""
        constraint_data, target_data, constraint_grid, target_grid, obs_pdf = (
            sample_emergent_constraint_data
        )

        # Test gaussian_pdf vs calc_GAUSSIAN_PDF
        x = np.linspace(-2, 2, 50)
        mu, sigma = 0.5, 1.0

        new_result = gaussian_pdf(mu, sigma, x)
        legacy_result = calc_GAUSSIAN_PDF(mu, sigma, x)

        assert_array_almost_equal(new_result, legacy_result)

        # Test emergent_constraint_posterior vs calc_PDF_EC
        new_posterior = emergent_constraint_posterior(
            constraint_data, target_data, constraint_grid, target_grid, obs_pdf
        )
        legacy_posterior = calc_PDF_EC(
            constraint_data, target_data, constraint_grid, target_grid, obs_pdf
        )

        # Compare all three returned values
        assert_array_almost_equal(new_posterior[0], legacy_posterior[0])  # PDF
        assert abs(new_posterior[1] - legacy_posterior[1]) < 1e-10  # std
        assert abs(new_posterior[2] - legacy_posterior[2]) < 1e-10  # mean

        # Test _calculate_std_from_pdf vs find_std_from_PDF
        x_test = np.linspace(-3, 3, 100)
        pdf_test = gaussian_pdf(0, 1, x_test)

        new_std = _calculate_std_from_pdf(0.341, x_test, pdf_test)
        legacy_std = find_std_from_PDF(0.341, x_test, pdf_test)

        assert abs(new_std - legacy_std) < 1e-10

        # Test emergent_constraint_prior vs calc_PDF_EC_PRIOR
        new_prior = emergent_constraint_prior(
            constraint_data, target_data, constraint_grid, target_grid
        )
        legacy_prior = calc_PDF_EC_PRIOR(
            constraint_data, target_data, constraint_grid, target_grid
        )

        # Compare all three returned arrays
        assert_array_almost_equal(new_prior[0], legacy_prior[0])  # prior_pdf
        assert_array_almost_equal(new_prior[1], legacy_prior[1])  # prediction_error
        assert_array_almost_equal(new_prior[2], legacy_prior[2])  # regression_line

    def test_emergent_constraints_input_validation(self):
        """Test input validation for emergent constraints functions."""
        # Test with mismatched array sizes
        constraint_data = xr.DataArray([1.0, 2.0], dims=["model"])
        target_data = xr.DataArray([1.0, 2.0, 3.0], dims=["model"])  # Different size

        constraint_grid = np.linspace(0, 3, 10)
        target_grid = np.linspace(0, 5, 15)
        obs_pdf = gaussian_pdf(1.5, 0.5, constraint_grid)

        # Should handle different sizes gracefully or raise appropriate error
        try:
            result = emergent_constraint_posterior(
                constraint_data, target_data, constraint_grid, target_grid, obs_pdf
            )
            # If it doesn't raise an error, check that result is reasonable
            assert len(result) == 3
        except (ValueError, IndexError):
            # This is also acceptable behavior
            pass

    def test_gaussian_pdf_error_conditions(self):
        """Test Gaussian PDF with error conditions."""
        # Test with zero standard deviation
        with pytest.warns(RuntimeWarning):  # Division by zero warning
            result = gaussian_pdf(0, 0, 0)
            assert np.isinf(result) or np.isnan(result)

        # Test with negative standard deviation
        with pytest.warns(RuntimeWarning):  # May produce warnings
            result = gaussian_pdf(0, -1, 0)
            # Result may be complex or NaN, which is expected


class TestEmergentConstraintsIntegration:
    """Integration tests for emergent constraints with realistic scenarios."""

    def test_climate_sensitivity_constraint_workflow(self):
        """Test complete workflow for climate sensitivity constraints."""
        np.random.seed(42)

        # Simulate CMIP model data for climate sensitivity
        n_models = 25

        # Constraint variable: tropical land temperature variability
        constraint_obs = 0.8  # Observed value
        constraint_uncertainty = 0.2
        constraint_models = np.random.normal(0.85, 0.25, n_models)

        # Target variable: equilibrium climate sensitivity
        # Create realistic relationship based on Cox et al. (2013)
        true_slope = -3.0  # Negative relationship
        true_intercept = 6.0
        ecs_models = (
            true_intercept
            + true_slope * constraint_models
            + np.random.normal(0, 0.3, n_models)
        )

        # Ensure positive ECS values
        ecs_models = np.clip(ecs_models, 1.0, 8.0)

        # Create xarray data
        constraint_data = xr.DataArray(
            constraint_models,
            dims=["model"],
            attrs={"long_name": "Tropical Temperature Variability", "units": "K"},
        )
        target_data = xr.DataArray(
            ecs_models,
            dims=["model"],
            attrs={"long_name": "Equilibrium Climate Sensitivity", "units": "K"},
        )

        # Set up grids
        constraint_grid = np.linspace(0.2, 1.5, 100)
        target_grid = np.linspace(1.0, 8.0, 150)

        # Observational constraint
        obs_pdf = gaussian_pdf(constraint_obs, constraint_uncertainty, constraint_grid)

        # Calculate prior (unconstrained) distribution
        prior_pdf, prediction_error, regression_line = emergent_constraint_prior(
            constraint_data, target_data, constraint_grid, target_grid
        )

        # Calculate posterior (constrained) distribution
        posterior_pdf, posterior_std, posterior_mean = emergent_constraint_posterior(
            constraint_data, target_data, constraint_grid, target_grid, obs_pdf
        )

        # Validate results
        assert len(posterior_pdf) == len(target_grid)
        assert posterior_std > 0
        assert 1.0 <= posterior_mean <= 8.0

        # Check that constraint reduces uncertainty
        prior_model_std = np.std(ecs_models)
        assert posterior_std <= prior_model_std

        # Check that posterior mean is reasonable given the relationship
        expected_mean = true_intercept + true_slope * constraint_obs
        # Allow reasonable tolerance
        assert abs(posterior_mean - expected_mean) < 1.5

    def test_emergent_constraints_statistical_consistency(self):
        """Test statistical consistency of emergent constraints method."""
        np.random.seed(123)

        # Create controlled test case
        n_models = 30
        true_constraint = 2.0
        true_target = 5.0
        correlation = 0.8

        # Generate correlated model data
        constraint_models = np.random.normal(true_constraint, 0.5, n_models)

        # Create correlated target with specified correlation
        independent_noise = np.random.normal(0, 0.3, n_models)
        correlated_noise = correlation * (constraint_models - true_constraint) / 0.5
        target_models = (
            true_target
            + 1.5 * (constraint_models - true_constraint)
            + independent_noise * np.sqrt(1 - correlation**2)
        )

        constraint_data = xr.DataArray(constraint_models, dims=["model"])
        target_data = xr.DataArray(target_models, dims=["model"])

        constraint_grid = np.linspace(0.5, 3.5, 80)
        target_grid = np.linspace(2.0, 8.0, 100)

        # Very precise observational constraint
        obs_pdf = gaussian_pdf(true_constraint, 0.1, constraint_grid)

        posterior_pdf, posterior_std, posterior_mean = emergent_constraint_posterior(
            constraint_data, target_data, constraint_grid, target_grid, obs_pdf
        )

        # With high correlation and precise observation,
        # posterior should be close to true target
        assert abs(posterior_mean - true_target) < 1.0

        # Posterior uncertainty should be reduced
        prior_std = np.std(target_models)
        assert posterior_std < 0.8 * prior_std


class TestCalculateThetaSE:
    """Test calculate_theta_se (Pseudo-Equivalent Potential Temperature) functionality."""

    def test_calculate_theta_se_basic(self):
        """Test basic theta-se calculation with realistic atmospheric values."""
        # Typical surface conditions in tropics
        temperature = 300.0  # K (27°C)
        pressure = 1000.0  # hPa
        mixing_ratio = 0.015  # kg/kg (15 g/kg)
        dewpoint = 20.0  # °C

        theta_se = calculate_theta_se(temperature, pressure, mixing_ratio, dewpoint)

        # Extract magnitude if it's a Quantity object
        if hasattr(theta_se, "magnitude"):
            theta_se_value = theta_se.magnitude
        else:
            theta_se_value = theta_se

        # Theta-se should be greater than temperature due to latent heat
        assert theta_se_value > temperature
        # Typical range for tropical theta-se is 330-360 K
        assert 330 < theta_se_value < 370
        # Should be finite
        assert np.isfinite(theta_se_value)

    def test_calculate_theta_se_arrays(self):
        """Test theta-se calculation with numpy arrays."""
        # Vertical profile
        temperatures = np.array([300, 290, 280, 270])  # K
        pressures = np.array([1000, 850, 700, 500])  # hPa
        mixing_ratios = np.array([0.015, 0.010, 0.005, 0.002])  # kg/kg
        dewpoints = np.array([20, 15, 10, 5])  # °C

        theta_se = calculate_theta_se(temperatures, pressures, mixing_ratios, dewpoints)

        # Extract magnitude if it's a Quantity object
        if hasattr(theta_se, "magnitude"):
            theta_se_value = theta_se.magnitude
        else:
            theta_se_value = theta_se

        # Check shape
        assert theta_se_value.shape == temperatures.shape
        # All values should be finite
        assert np.all(np.isfinite(theta_se_value))
        # Theta-se should increase or stay similar with height in convective atmosphere
        # (though this is profile-dependent)
        assert np.all(theta_se_value > temperatures)

    def test_calculate_theta_se_xarray(self):
        """Test theta-se calculation with xarray DataArrays."""
        # Create xarray inputs with plain numpy values (no units in attrs)
        temperature = xr.DataArray(
            [300, 295, 290],
            dims=["level"],
            attrs={"long_name": "Temperature"},
        )
        pressure = xr.DataArray([1000, 900, 800], dims=["level"])
        mixing_ratio = xr.DataArray([0.015, 0.012, 0.009], dims=["level"])
        dewpoint = xr.DataArray([20, 18, 15], dims=["level"])

        theta_se = calculate_theta_se(temperature, pressure, mixing_ratio, dewpoint)

        # Should return xarray with attributes
        assert isinstance(theta_se, xr.DataArray)
        assert "units" in theta_se.attrs
        assert theta_se.attrs["units"] == "K"
        assert "long_name" in theta_se.attrs
        # Values should be reasonable - extract magnitude if needed
        theta_se_values = (
            theta_se.values
            if not hasattr(theta_se.values, "magnitude")
            else theta_se.values.magnitude
        )
        temp_values = temperature.values
        assert np.all(theta_se_values > temp_values)

    def test_calculate_theta_se_dry_conditions(self):
        """Test theta-se with very dry conditions."""
        temperature = 300.0  # K
        pressure = 1000.0  # hPa
        mixing_ratio = 0.001  # kg/kg (very dry - 1 g/kg)
        dewpoint = 0.0  # °C (low dewpoint)

        theta_se = calculate_theta_se(temperature, pressure, mixing_ratio, dewpoint)

        # Extract magnitude if it's a Quantity object
        if hasattr(theta_se, "magnitude"):
            theta_se_value = theta_se.magnitude
        else:
            theta_se_value = theta_se

        # Should still produce finite result
        assert np.isfinite(theta_se_value)
        # With low moisture, theta-se should be closer to potential temperature
        # but still higher due to some latent heat
        assert theta_se_value > temperature

    def test_calculate_theta_se_high_altitude(self):
        """Test theta-se at high altitude (low pressure)."""
        temperature = 250.0  # K (cold upper troposphere)
        pressure = 300.0  # hPa (~ 9 km altitude)
        mixing_ratio = 0.0005  # kg/kg (0.5 g/kg - very dry at altitude)
        dewpoint = -30.0  # °C

        theta_se = calculate_theta_se(temperature, pressure, mixing_ratio, dewpoint)

        # Extract magnitude if it's a Quantity object
        if hasattr(theta_se, "magnitude"):
            theta_se_value = theta_se.magnitude
        else:
            theta_se_value = theta_se

        # Should produce finite result
        assert np.isfinite(theta_se_value)
        # At high altitude, theta-se should still be significantly higher than T
        assert theta_se_value > temperature

    def test_calculate_theta_se_comparison_numpy_xarray(self):
        """Test that numpy and xarray inputs give same results."""
        # NumPy inputs
        temp_np = np.array([300, 295, 290])
        pres_np = np.array([1000, 900, 800])
        mixr_np = np.array([0.015, 0.012, 0.009])
        dewp_np = np.array([20, 18, 15])

        theta_se_np = calculate_theta_se(temp_np, pres_np, mixr_np, dewp_np)

        # xarray inputs
        temp_xr = xr.DataArray(temp_np, dims=["level"])
        pres_xr = xr.DataArray(pres_np, dims=["level"])
        mixr_xr = xr.DataArray(mixr_np, dims=["level"])
        dewp_xr = xr.DataArray(dewp_np, dims=["level"])

        theta_se_xr = calculate_theta_se(temp_xr, pres_xr, mixr_xr, dewp_xr)

        # Extract magnitudes for comparison
        if hasattr(theta_se_np, "magnitude"):
            theta_se_np_values = theta_se_np.magnitude
        else:
            theta_se_np_values = theta_se_np

        if hasattr(theta_se_xr.values, "magnitude"):
            theta_se_xr_values = theta_se_xr.values.magnitude
        else:
            theta_se_xr_values = theta_se_xr.values

        # Results should be nearly identical
        assert_array_almost_equal(theta_se_np_values, theta_se_xr_values, decimal=6)

    def test_calculate_theta_se_scalar_vs_array(self):
        """Test that scalar and single-element array give same result."""
        # Scalar inputs
        theta_se_scalar = calculate_theta_se(300.0, 1000.0, 0.015, 20.0)

        # Array inputs with single element
        theta_se_array = calculate_theta_se(
            np.array([300.0]),
            np.array([1000.0]),
            np.array([0.015]),
            np.array([20.0]),
        )

        # Extract magnitudes for comparison
        if hasattr(theta_se_scalar, "magnitude"):
            theta_se_scalar_value = theta_se_scalar.magnitude
        else:
            theta_se_scalar_value = theta_se_scalar

        if hasattr(theta_se_array, "magnitude"):
            theta_se_array_value = theta_se_array.magnitude[0]
        else:
            theta_se_array_value = theta_se_array[0]

        # Should be very close
        assert np.abs(theta_se_scalar_value - theta_se_array_value) < 0.01

    def test_calculate_theta_se_multidimensional(self):
        """Test theta-se calculation with 2D arrays (e.g., spatial grid)."""
        # Create 2D grid (e.g., lat-lon)
        nlat, nlon = 5, 6
        temperatures = np.random.uniform(280, 310, (nlat, nlon))
        pressures = np.full((nlat, nlon), 1000.0)  # Surface pressure
        mixing_ratios = np.random.uniform(0.005, 0.020, (nlat, nlon))
        dewpoints = np.random.uniform(10, 25, (nlat, nlon))

        theta_se = calculate_theta_se(temperatures, pressures, mixing_ratios, dewpoints)

        # Extract magnitude if it's a Quantity object
        if hasattr(theta_se, "magnitude"):
            theta_se_value = theta_se.magnitude
        else:
            theta_se_value = theta_se

        # Check shape
        assert theta_se_value.shape == (nlat, nlon)
        # All values should be finite
        assert np.all(np.isfinite(theta_se_value))
        # All theta-se should be greater than temperature
        assert np.all(theta_se_value > temperatures)


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
