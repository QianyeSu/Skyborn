"""
Tests for Lanczos filter module.
"""

import numpy as np
import pytest

from skyborn.calc.filter import (
    lanczos_bandpass,
    lanczos_filter,
    lanczos_highpass,
    lanczos_lowpass,
    lanczos_weights,
)

# Try importing xarray for xarray tests
try:
    import xarray as xr

    HAS_XARRAY = True
except ImportError:
    HAS_XARRAY = False


class TestLanczosWeights:
    """Test weight calculation."""

    def test_basic_lowpass_weights(self):
        """Test basic low-pass weight generation."""
        weights = lanczos_weights(cutoff_freq=0.1, window=61, pass_type="low")

        assert len(weights) == 61
        assert weights.dtype == np.float64
        # Low-pass weights should sum to 1
        np.testing.assert_allclose(weights.sum(), 1.0, rtol=1e-10)
        # Should be symmetric
        np.testing.assert_array_almost_equal(weights[:30], weights[-30:][::-1])

    def test_basic_highpass_weights(self):
        """Test basic high-pass weight generation."""
        weights = lanczos_weights(cutoff_freq=0.1, window=61, pass_type="high")

        assert len(weights) == 61
        # High-pass weights should sum to ~0
        np.testing.assert_allclose(weights.sum(), 0.0, atol=1e-10)
        # Should be symmetric
        np.testing.assert_array_almost_equal(weights[:30], weights[-30:][::-1])

    def test_weight_normalization(self):
        """Test that weights are properly normalized."""
        # Low-pass should sum to 1
        weights_lp = lanczos_weights(0.05, 81, "low")
        np.testing.assert_allclose(weights_lp.sum(), 1.0)

        # High-pass should sum to 0
        weights_hp = lanczos_weights(0.05, 81, "high")
        np.testing.assert_allclose(weights_hp.sum(), 0.0, atol=1e-10)

    def test_invalid_window(self):
        """Test that even windows are rejected."""
        with pytest.raises(ValueError, match="window size must be odd"):
            lanczos_weights(0.1, window=60, pass_type="low")

    def test_small_window(self):
        """Test that small windows are rejected."""
        with pytest.raises(ValueError, match="window size must be >= 3"):
            lanczos_weights(0.1, window=1, pass_type="low")

    def test_invalid_cutoff_freq(self):
        """Test that invalid cutoff frequencies are rejected."""
        with pytest.raises(ValueError, match="cutoff_freq must be in"):
            lanczos_weights(0.6, window=61, pass_type="low")

        with pytest.raises(ValueError, match="cutoff_freq must be in"):
            lanczos_weights(-0.1, window=61, pass_type="low")

    def test_invalid_pass_type(self):
        """Test that invalid pass types are rejected."""
        with pytest.raises(ValueError, match="pass_type must be"):
            lanczos_weights(0.1, window=61, pass_type="bandpass")


class TestLanczosFilter1D:
    """Test 1D filtering on numpy arrays."""

    def test_lowpass_constant_signal(self):
        """Test that low-pass filter preserves constant signal."""
        data = np.ones(1000)
        filtered = lanczos_lowpass(data, cutoff_freq=0.1, window=61)

        # Should preserve constant (DC component)
        np.testing.assert_allclose(filtered, 1.0, rtol=1e-6)

    def test_highpass_constant_signal(self):
        """Test that high-pass filter removes constant signal."""
        data = np.ones(1000) * 5.0
        filtered = lanczos_highpass(data, cutoff_freq=0.1, window=61)

        # Should remove constant (DC component)
        np.testing.assert_allclose(filtered, 0.0, atol=1e-10)

    def test_lowpass_removes_high_freq(self):
        """Test that low-pass removes high-frequency noise."""
        # Create signal: low-freq sine + high-freq noise
        t = np.linspace(0, 10 * np.pi, 1000)
        signal = np.sin(t)  # Low frequency
        noise = 0.5 * np.sin(20 * t)  # High frequency
        data = signal + noise

        # Filter should remove high-freq noise
        filtered = lanczos_lowpass(data, cutoff_freq=0.05, window=81)

        # Filtered should be closer to pure signal (check in middle to avoid edges)
        center = slice(100, 900)
        signal_error = np.abs(data[center] - signal[center]).mean()
        filtered_error = np.abs(filtered[center] - signal[center]).mean()
        assert filtered_error < signal_error * 0.3  # Much less error

    def test_highpass_removes_trend(self):
        """Test that high-pass removes low-frequency trend."""
        # Create signal: linear trend + high-freq oscillation
        t = np.linspace(0, 10, 1000)
        trend = 2 * t + 5  # Linear trend
        oscillation = 0.5 * np.sin(50 * t)  # High frequency
        data = trend + oscillation

        # Filter should remove trend
        filtered = lanczos_highpass(data, cutoff_freq=0.05, window=81)

        # Mean of high-passed signal should be near zero
        assert np.abs(filtered[100:900].mean()) < 0.1
        # Should preserve oscillation amplitude (roughly)
        assert 0.3 < np.abs(filtered[100:900]).mean() < 0.7

    def test_bandpass_extracts_frequency_band(self):
        """Test that band-pass extracts specific frequency range."""
        # Create signal with 3 frequency components
        t = np.linspace(0, 100 * np.pi, 5000)
        low_freq = np.sin(t)  # f ~ 0.01
        mid_freq = 2 * np.sin(10 * t)  # f ~ 0.1
        high_freq = 0.5 * np.sin(50 * t)  # f ~ 0.5
        data = low_freq + mid_freq + high_freq

        # Extract mid-frequency band
        filtered = lanczos_bandpass(data, freq_low=0.05, freq_high=0.15, window=121)

        # Check that mid-frequency dominates in filtered signal
        # (crude test - just verify amplitude is reasonable)
        center = slice(500, 4500)
        assert 1.0 < filtered[center].std() < 3.0

    def test_boundary_reflect(self):
        """Test reflect boundary handling."""
        data = np.arange(100, dtype=float)
        filtered = lanczos_lowpass(data, cutoff_freq=0.1, window=21, boundary="reflect")

        # Should produce smooth result without edge artifacts
        assert np.all(np.isfinite(filtered))

    def test_boundary_periodic(self):
        """Test periodic boundary handling."""
        # Create periodic signal
        t = np.linspace(0, 2 * np.pi, 100, endpoint=False)
        data = np.sin(t)
        filtered = lanczos_lowpass(
            data, cutoff_freq=0.1, window=21, boundary="periodic"
        )

        # Should be periodic: start ≈ end
        np.testing.assert_allclose(filtered[0], filtered[-1], atol=0.1)

    def test_auto_window_size(self):
        """Test automatic window size determination."""
        data = np.random.randn(1000)

        # Without specifying window
        filtered = lanczos_lowpass(data, cutoff_freq=0.1, window=None)

        assert np.all(np.isfinite(filtered[100:900]))


class TestLanczosFilter2D:
    """Test 2D filtering on numpy arrays."""

    def test_2d_lowpass_constant_field(self):
        """Test that 2D low-pass preserves constant field."""
        field = np.ones((100, 200)) * 3.0
        filtered = lanczos_lowpass(field, cutoff_freq=0.1, window=21, dim=(0, 1))

        # Should preserve constant
        np.testing.assert_allclose(filtered, 3.0, rtol=1e-6)

    def test_2d_separable_filtering(self):
        """Test that 2D filter is separable (same result as two 1D filters)."""
        # Create 2D field with noise
        np.random.seed(42)
        field = np.random.randn(80, 120)

        cutoff = 0.1
        window = 21

        # Apply 2D filter
        filtered_2d = lanczos_lowpass(
            field, cutoff_freq=cutoff, window=window, dim=(0, 1)
        )

        # Apply two 1D filters sequentially
        filtered_1d = lanczos_lowpass(field, cutoff_freq=cutoff, window=window, dim=0)
        filtered_1d = lanczos_lowpass(
            filtered_1d, cutoff_freq=cutoff, window=window, dim=1
        )

        # Should be identical (separable filtering)
        np.testing.assert_allclose(filtered_2d, filtered_1d, rtol=1e-10)

    def test_2d_spatial_smoothing(self):
        """Test spatial smoothing on checkerboard pattern."""
        # Create checkerboard pattern (high spatial frequency)
        x, y = np.meshgrid(np.arange(100), np.arange(100))
        checkerboard = ((x // 5 + y // 5) % 2).astype(float)

        # Smooth with low-pass filter
        smoothed = lanczos_lowpass(
            checkerboard, cutoff_freq=0.05, window=31, dim=(0, 1)
        )

        # Smoothed version should have much less variance
        assert smoothed.std() < checkerboard.std() * 0.5


@pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
class TestLanczosFilterXarray:
    """Test filtering on xarray DataArrays."""

    def test_xarray_basic_filtering(self):
        """Test basic xarray filtering with dimension names."""
        # Create simple xarray
        time = np.arange(1000)
        data = np.sin(2 * np.pi * time / 50) + np.random.randn(1000) * 0.1
        da = xr.DataArray(data, coords={"time": time}, dims=["time"])

        # Filter along time dimension
        filtered = lanczos_lowpass(da, cutoff_freq=0.05, window=61, dim="time")

        assert isinstance(filtered, xr.DataArray)
        assert filtered.dims == da.dims
        assert "lanczos_filtered" in filtered.attrs

    def test_xarray_multidim_filtering(self):
        """Test multi-dimensional xarray filtering."""
        # Create 3D data (time, lat, lon)
        nt, nlat, nlon = 200, 50, 100
        data = np.random.randn(nt, nlat, nlon)
        da = xr.DataArray(
            data,
            coords={
                "time": np.arange(nt),
                "lat": np.linspace(-90, 90, nlat),
                "lon": np.linspace(0, 360, nlon),
            },
            dims=["time", "lat", "lon"],
        )

        # Filter along multiple dimensions
        filtered = lanczos_lowpass(
            da,
            cutoff_freq={"time": 0.05, "lat": 0.1},
            window={"time": 31, "lat": 21},
            dim=["time", "lat"],
        )

        assert isinstance(filtered, xr.DataArray)
        assert filtered.shape == da.shape
        assert "lanczos_cutoff_freq" in filtered.attrs

    def test_xarray_preserves_coords(self):
        """Test that coordinates are preserved."""
        # Create data with rich coordinates
        time = np.arange(365)
        lat = np.linspace(-90, 90, 50)
        lon = np.linspace(0, 360, 100)
        data = np.random.randn(365, 50, 100)

        da = xr.DataArray(
            data,
            coords={
                "time": time,
                "lat": lat,
                "lon": lon,
            },
            dims=["time", "lat", "lon"],
            attrs={"units": "K", "description": "temperature"},
        )

        filtered = lanczos_lowpass(da, cutoff_freq=0.05, window=31, dim="time")

        # Check coordinates preserved
        assert np.array_equal(filtered.time, da.time)
        assert np.array_equal(filtered.lat, da.lat)
        assert np.array_equal(filtered.lon, da.lon)
        # Original attrs should be preserved
        assert filtered.attrs["units"] == "K"


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_small_data_array(self):
        """Test filtering on small array."""
        data = np.array([1, 2, 3, 4, 5], dtype=float)
        # Window should adapt to data size
        filtered = lanczos_lowpass(data, cutoff_freq=0.2, window=3)
        assert len(filtered) == len(data)

    def test_single_value_array(self):
        """Test filtering on single-value array."""
        data = np.array([5.0])
        filtered = lanczos_lowpass(data, cutoff_freq=0.1, window=3)
        assert filtered.shape == data.shape

    def test_nan_handling(self):
        """Test behavior with NaN values."""
        data = np.ones(100)
        data[50] = np.nan

        # Filter propagates NaNs
        filtered = lanczos_lowpass(data, cutoff_freq=0.1, window=21)
        # Near the NaN, filtered values will also be NaN
        assert np.isnan(filtered[50])

    def test_negative_axis_indexing(self):
        """Test negative axis indexing."""
        data = np.random.randn(50, 100, 200)

        # Filter on last axis using negative index
        filtered = lanczos_lowpass(data, cutoff_freq=0.1, window=21, dim=-1)

        assert filtered.shape == data.shape

    def test_bandpass_invalid_frequencies(self):
        """Test that invalid band-pass frequencies are rejected."""
        data = np.random.randn(100)

        with pytest.raises(ValueError, match="freq_low.*must be.*freq_high"):
            lanczos_bandpass(data, freq_low=0.3, freq_high=0.1)


class TestNumericalAccuracy:
    """Test numerical accuracy and precision."""

    def test_impulse_response(self):
        """Test impulse response matches filter weights."""
        # Create impulse signal
        impulse = np.zeros(1000)
        impulse[500] = 1.0

        # Filter with specific weights
        window = 61
        cutoff = 0.1
        filtered = lanczos_lowpass(impulse, cutoff_freq=cutoff, window=window)

        # Get expected weights
        weights = lanczos_weights(cutoff, window, "low")

        # Response should match weights (centered at impulse position)
        half_window = window // 2
        response = filtered[500 - half_window : 500 + half_window + 1]

        np.testing.assert_allclose(response, weights, rtol=1e-6)

    def test_linearity(self):
        """Test that filter is linear (scaling property)."""
        data = np.random.randn(500)

        filtered1 = lanczos_lowpass(data, cutoff_freq=0.1, window=41)
        filtered2 = lanczos_lowpass(3.0 * data, cutoff_freq=0.1, window=41)

        # Scaling property: L[a*x] = a*L[x]
        np.testing.assert_allclose(filtered2, 3.0 * filtered1, rtol=1e-10)

    def test_additivity(self):
        """Test that filter is additive."""
        data1 = np.random.randn(500)
        data2 = np.random.randn(500)

        filtered1 = lanczos_lowpass(data1, cutoff_freq=0.1, window=41)
        filtered2 = lanczos_lowpass(data2, cutoff_freq=0.1, window=41)
        filtered_sum = lanczos_lowpass(data1 + data2, cutoff_freq=0.1, window=41)

        # Additivity: L[x+y] = L[x] + L[y]
        np.testing.assert_allclose(filtered_sum, filtered1 + filtered2, rtol=1e-10)
