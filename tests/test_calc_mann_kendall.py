"""
Tests for Mann-Kendall trend analysis module.

This module tests both the numpy and xarray implementations of Mann-Kendall
trend detection algorithms.
"""

import numpy as np
import pytest
from skyborn.calc.mann_kendall import (
    mann_kendall_test,
    mann_kendall_multidim_numpy,
    trend_analysis,
)

try:
    import xarray as xr
    from skyborn.calc.mann_kendall import mann_kendall_xarray

    HAS_XARRAY = True
except ImportError:
    HAS_XARRAY = False


class TestMannKendall:
    """Test suite for Mann-Kendall functionality."""

    def test_basic_trend_detection(self):
        """Test basic trend detection with known trends."""
        # Create data with known positive trend
        n = 50
        x = np.arange(n)
        y = 2 * x + np.random.randn(n) * 0.5  # Strong positive trend

        result = mann_kendall_test(y)

        assert result["h"] == True, "Should detect significant trend"
        assert result["trend"] > 0, "Should detect positive trend"
        assert 0 < result["p"] < 0.05, "Should have significant p-value"
        assert result["z"] > 1.96, "Should have significant z-score"

    def test_no_trend_detection(self):
        """Test detection of no trend in random data."""
        # Create random data with no trend
        np.random.seed(42)
        y = np.random.randn(50)

        result = mann_kendall_test(y)

        # With random data, should generally not detect trend
        assert abs(result["trend"]) < 1, "Trend should be small for random data"
        assert result["p"] > 0.01, "P-value should be large for random data"

    def test_negative_trend(self):
        """Test detection of negative trends."""
        n = 50
        x = np.arange(n)
        y = -1.5 * x + np.random.randn(n) * 0.3  # Strong negative trend

        result = mann_kendall_test(y)

        assert result["h"] == True, "Should detect significant trend"
        assert result["trend"] < 0, "Should detect negative trend"
        assert result["z"] < -1.96, "Should have significant negative z-score"

    def test_missing_values(self):
        """Test handling of missing values."""
        # Create data with NaN values
        n = 50
        x = np.arange(n)
        y = 2 * x + np.random.randn(n) * 0.5
        y[10:15] = np.nan  # Add missing values

        result = mann_kendall_test(y)

        # Should still detect trend despite missing values
        assert not np.isnan(result["trend"]), "Should handle missing values"
        assert result["h"] == True, "Should detect trend despite missing values"

    def test_insufficient_data(self):
        """Test behavior with insufficient data."""
        # Too few data points
        y = np.array([1, 2])

        result = mann_kendall_test(y)

        assert np.isnan(result["trend"]), "Should return NaN for insufficient data"
        assert result["h"] == False, "Should not detect trend with insufficient data"

    def test_modified_mann_kendall(self):
        """Test modified Mann-Kendall for autocorrelated data."""
        # Create autocorrelated data with trend
        n = 100
        trend = np.linspace(0, 10, n)

        # Add autocorrelated noise
        noise = np.random.randn(n)
        for i in range(1, n):
            noise[i] += 0.7 * noise[i - 1]  # AR(1) process

        y = trend + noise

        # Compare standard vs modified test
        result_standard = mann_kendall_test(y, modified=False)
        result_modified = mann_kendall_test(y, modified=True)

        # Modified test should be more conservative (larger p-value)
        assert (
            result_modified["p"] >= result_standard["p"]
        ), "Modified test should be more conservative"

    def test_multidimensional_numpy(self):
        """Test multidimensional numpy implementation."""
        # Create 3D data: (time, lat, lon)
        time_steps, nlat, nlon = 50, 10, 15

        # Create spatial pattern with different trends
        trends = np.random.randn(nlat, nlon) * 0.1
        data = np.zeros((time_steps, nlat, nlon))

        for t in range(time_steps):
            data[t] = trends * t + np.random.randn(nlat, nlon) * 0.5

        results = mann_kendall_multidim_numpy(data, time_axis=0)

        # Check output shapes
        assert results["trend"].shape == (nlat, nlon)
        assert results["h"].shape == (nlat, nlon)
        assert results["p"].shape == (nlat, nlon)

        # Check that detected trends are reasonable
        detected_trends = results["trend"][results["h"]]
        if len(detected_trends) > 0:
            # Should be correlated with input trends where significant
            significant_mask = results["h"]
            input_trends_sig = trends[significant_mask]
            output_trends_sig = results["trend"][significant_mask]

            # Should have same sign for most significant trends
            same_sign = np.sign(input_trends_sig) == np.sign(output_trends_sig)
            assert np.mean(same_sign) > 0.7, "Should detect correct trend signs"

    def test_different_time_axes(self):
        """Test with time along different axes."""
        # Create data with time along axis 1
        nlat, time_steps, nlon = 8, 30, 12
        data = np.random.randn(nlat, time_steps, nlon)

        # Add trend along time axis (axis 1)
        for t in range(time_steps):
            data[:, t, :] += t * 0.1

        results = mann_kendall_multidim_numpy(data, time_axis=1)

        assert results["trend"].shape == (nlat, nlon)
        # Should detect positive trends
        assert np.mean(results["trend"][results["h"]]) > 0

    def test_chunked_processing(self):
        """Test chunked processing for memory efficiency."""
        # Large spatial dimensions
        time_steps, nlat, nlon = 20, 100, 150
        data = np.random.randn(time_steps, nlat, nlon)

        # Test with small chunk size
        results = mann_kendall_multidim_numpy(data, time_axis=0, chunk_size=1000)

        assert results["trend"].shape == (nlat, nlon)
        assert np.all(np.isfinite(results["trend"]) | np.isnan(results["trend"]))

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not available")
    def test_xarray_interface(self):
        """Test xarray interface."""
        # Create xarray DataArray
        time = np.arange(40)
        lat = np.linspace(-60, 60, 8)
        lon = np.linspace(0, 360, 12, endpoint=False)

        # Create data with spatial trend pattern
        data = np.random.randn(len(time), len(lat), len(lon))
        for i, t in enumerate(time):
            data[i] += t * 0.05  # Add trend

        da = xr.DataArray(
            data,
            dims=["time", "lat", "lon"],
            coords={"time": time, "lat": lat, "lon": lon},
        )

        result = mann_kendall_xarray(da, dim="time")

        # Check output
        assert isinstance(result, xr.Dataset)
        assert "trend" in result.data_vars
        assert "h" in result.data_vars
        assert "p" in result.data_vars

        # Check coordinates
        assert "lat" in result.coords
        assert "lon" in result.coords
        assert "time" not in result.coords  # Should be removed

        # Check shapes
        assert result["trend"].shape == (len(lat), len(lon))

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not available")
    def test_xarray_with_missing_data(self):
        """Test xarray interface with missing data."""
        time = np.arange(30)
        lat = np.linspace(-30, 30, 5)
        lon = np.linspace(0, 180, 6)

        data = np.random.randn(len(time), len(lat), len(lon))
        data[10:15, 1:3, 2:4] = np.nan  # Add missing region

        da = xr.DataArray(
            data,
            dims=["time", "lat", "lon"],
            coords={"time": time, "lat": lat, "lon": lon},
        )

        result = mann_kendall_xarray(da, dim="time", use_dask=False)

        # Should handle missing data gracefully
        assert isinstance(result, xr.Dataset)
        # Missing regions should have NaN or False for significance
        missing_region = result["h"].isel(lat=slice(1, 3), lon=slice(2, 4))
        # Should not crash and produce reasonable results

    def test_trend_analysis_unified_interface(self):
        """Test unified trend_analysis interface."""
        # Test with numpy array
        data_np = np.random.randn(40, 5, 8)
        for t in range(40):
            data_np[t] += t * 0.02

        result_np = trend_analysis(data_np, time_axis=0)
        assert isinstance(result_np, dict)
        assert "trend" in result_np

        if HAS_XARRAY:
            # Test with xarray
            da = xr.DataArray(
                data_np,
                dims=["time", "y", "x"],
                coords={"time": np.arange(40), "y": np.arange(5), "x": np.arange(8)},
            )

            result_xr = trend_analysis(da, time_axis="time")
            assert isinstance(result_xr, xr.Dataset)
            assert "trend" in result_xr.data_vars

    def test_different_methods(self):
        """Test different slope calculation methods."""
        # Create data with clear trend
        n = 50
        x = np.arange(n)
        y = 2 * x + np.random.randn(n) * 0.5

        # Test different methods
        result_theil = mann_kendall_test(y, method="theilslopes")
        result_linreg = mann_kendall_test(y, method="linregress")
        result_auto = mann_kendall_test(y, method="auto")

        # All should detect positive trend
        assert result_theil["trend"] > 0
        assert result_linreg["trend"] > 0
        assert result_auto["trend"] > 0

        # Results should be similar but not identical
        assert abs(result_theil["trend"] - result_linreg["trend"]) < 1.0

    def test_performance_comparison(self):
        """Test performance comparison between implementations."""
        # Create moderately large dataset
        time_steps, nlat, nlon = 100, 50, 60
        data = np.random.randn(time_steps, nlat, nlon)

        # Add some spatial trends
        for i in range(nlat):
            for j in range(nlon):
                trend_strength = (i - nlat // 2) * (j - nlon // 2) * 0.001
                data[:, i, j] += np.arange(time_steps) * trend_strength

        # Test numpy implementation
        import time

        start_time = time.time()
        results_np = mann_kendall_multidim_numpy(data, time_axis=0)
        numpy_time = time.time() - start_time

        # Verify results are reasonable
        assert results_np["trend"].shape == (nlat, nlon)
        assert np.any(results_np["h"])  # Should detect some trends

        print(f"Numpy implementation time: {numpy_time:.2f} seconds")
        print(f"Points processed: {nlat * nlon}")
        print(f"Processing rate: {(nlat * nlon) / numpy_time:.0f} points/second")


if __name__ == "__main__":
    # Run basic tests
    test_suite = TestMannKendall()

    print("Running Mann-Kendall tests...")

    test_suite.test_basic_trend_detection()
    print("✓ Basic trend detection test passed")

    test_suite.test_no_trend_detection()
    print("✓ No trend detection test passed")

    test_suite.test_negative_trend()
    print("✓ Negative trend test passed")

    test_suite.test_missing_values()
    print("✓ Missing values test passed")

    test_suite.test_multidimensional_numpy()
    print("✓ Multidimensional numpy test passed")

    test_suite.test_performance_comparison()
    print("✓ Performance test completed")

    if HAS_XARRAY:
        test_suite.test_xarray_interface()
        print("✓ xarray interface test passed")

    print("\nAll Mann-Kendall tests passed!")
