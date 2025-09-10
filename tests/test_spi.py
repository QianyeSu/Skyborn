"""
Tests for the Standardized Precipitation Index (SPI) module.

This module tests both the core SPI calculation functions and the 
xarray integration interface.
"""

import numpy as np
import pytest
import sys
import os

# Add src to path for testing
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

try:
    import xarray as xr
    HAS_XARRAY = True
except ImportError:
    HAS_XARRAY = False
    xr = None

from numpy.testing import assert_array_almost_equal, assert_allclose
from scipy import stats

from skyborn.calc.spi.core import (
    standardized_precipitation_index,
    spi,
    _gamma_fit_vectorized,
    _calculate_spi_values,
    _rolling_sum
)

if HAS_XARRAY:
    from skyborn.calc.spi.xarray import spi_xarray, spi_dataset


class TestSPICore:
    """Test core SPI calculation functions."""
    
    def test_rolling_sum_basic(self):
        """Test basic rolling sum functionality."""
        data = np.array([1, 2, 3, 4, 5])
        
        # 3-period rolling sum
        result = _rolling_sum(data, window=3, axis=0)
        expected = np.array([np.nan, np.nan, 6, 9, 12])  # [1+2+3, 2+3+4, 3+4+5]
        
        # Compare non-NaN values
        valid_mask = ~np.isnan(expected)
        assert_allclose(result[valid_mask], expected[valid_mask])
        
        # Check NaN positions
        assert np.all(np.isnan(result[~valid_mask]))
    
    def test_rolling_sum_2d(self):
        """Test rolling sum with 2D data."""
        # Create 2D data (time, space)
        data = np.array([
            [1, 10],
            [2, 20], 
            [3, 30],
            [4, 40]
        ])
        
        result = _rolling_sum(data, window=2, axis=0)
        
        # Expected: first row should be NaN, then rolling sums
        expected = np.array([
            [np.nan, np.nan],
            [3, 30],      # [1+2, 10+20]
            [5, 50],      # [2+3, 20+30]
            [7, 70]       # [3+4, 30+40]
        ])
        
        # Test valid values
        valid_mask = ~np.isnan(expected)
        assert_allclose(result[valid_mask], expected[valid_mask])
    
    def test_rolling_sum_window_1(self):
        """Test rolling sum with window=1 (should return original data)."""
        data = np.array([1, 2, 3, 4, 5])
        result = _rolling_sum(data, window=1, axis=0)
        assert_array_almost_equal(result, data)
    
    def test_gamma_fit_vectorized(self):
        """Test gamma distribution fitting."""
        # Create synthetic gamma-distributed data
        np.random.seed(42)
        true_alpha, true_scale = 2.0, 3.0
        
        # Generate data for 2 spatial points
        data = np.random.gamma(true_alpha, true_scale, size=(100, 2))
        
        alpha, beta = _gamma_fit_vectorized(data, axis=0)
        
        # Check that fitted parameters are reasonable
        assert alpha.shape == (2,)
        assert beta.shape == (2,)
        
        # Parameters should be positive and finite
        assert np.all(alpha > 0)
        assert np.all(beta > 0)
        assert np.all(np.isfinite(alpha))
        assert np.all(np.isfinite(beta))
        
        # Fitted parameters should be reasonably close to true values
        # (allowing for some variation due to random sampling)
        assert np.all(np.abs(alpha - true_alpha) < 1.0)
        assert np.all(np.abs(beta - true_scale) < 2.0)
    
    def test_gamma_fit_with_zeros(self):
        """Test gamma fitting with zero precipitation values."""
        # Create data with some zeros
        data = np.array([
            [0, 1, 2, 0, 3, 4, 0, 5],  # Mix of zeros and positive values
            [1, 2, 3, 4, 5, 6, 7, 8]   # All positive values
        ]).T
        
        alpha, beta = _gamma_fit_vectorized(data, axis=0)
        
        # Should still produce reasonable fits
        assert alpha.shape == (2,)
        assert beta.shape == (2,)
        assert np.all(alpha > 0)
        assert np.all(beta > 0)
        assert np.all(np.isfinite(alpha))
        assert np.all(np.isfinite(beta))
    
    def test_spi_calculation_basic(self):
        """Test basic SPI calculation."""
        # Create synthetic precipitation data
        np.random.seed(42)
        # Generate gamma-distributed precipitation for 5 years (60 months)
        precip = np.random.gamma(2, 2, size=(60,))
        
        # Calculate SPI with 1-month time scale
        spi_values = standardized_precipitation_index(precip, time_scale=1, axis=0)
        
        # Check basic properties
        assert spi_values.shape == precip.shape
        assert np.isfinite(spi_values).sum() > 40  # Most values should be finite
        
        # SPI should have approximately mean=0, std=1 for sufficiently long series
        valid_spi = spi_values[np.isfinite(spi_values)]
        if len(valid_spi) > 30:
            assert abs(np.mean(valid_spi)) < 0.5  # Mean should be close to 0
            assert abs(np.std(valid_spi) - 1.0) < 0.5  # Std should be close to 1
    
    def test_spi_multidimensional(self):
        """Test SPI calculation with multi-dimensional data."""
        # Create 3D precipitation data (time, lat, lon)
        np.random.seed(42)
        precip = np.random.gamma(2, 2, size=(48, 5, 4))  # 4 years, 5x4 grid
        
        # Calculate 3-month SPI
        spi_values = standardized_precipitation_index(precip, time_scale=3, axis=0)
        
        # Check output shape
        assert spi_values.shape == precip.shape
        
        # Check that we have reasonable number of finite values
        # (some will be NaN due to rolling window at beginning)
        n_finite = np.isfinite(spi_values).sum()
        expected_min_finite = 40 * 5 * 4  # At least 40 time steps for each grid point
        assert n_finite >= expected_min_finite
    
    def test_spi_different_time_scales(self):
        """Test SPI calculation with different time scales."""
        np.random.seed(42)
        precip = np.random.gamma(2, 2, size=(120,))  # 10 years monthly data
        
        # Test different time scales
        for time_scale in [1, 3, 6, 12]:
            spi_values = standardized_precipitation_index(precip, time_scale=time_scale, axis=0)
            
            assert spi_values.shape == precip.shape
            
            # Check that we lose some data at the beginning due to rolling window
            if time_scale > 1:
                # First (time_scale - 1) values should be NaN
                assert np.all(np.isnan(spi_values[:time_scale-1]))
            
            # Later values should be mostly finite
            later_values = spi_values[time_scale+10:]  # Skip initial period
            finite_ratio = np.isfinite(later_values).mean()
            assert finite_ratio > 0.8  # At least 80% should be finite
    
    def test_spi_alias(self):
        """Test that spi is an alias for standardized_precipitation_index."""
        np.random.seed(42)
        precip = np.random.gamma(2, 2, size=(60,))
        
        spi1 = standardized_precipitation_index(precip, time_scale=3)
        spi2 = spi(precip, time_scale=3)
        
        assert_array_almost_equal(spi1, spi2)
    
    def test_invalid_inputs(self):
        """Test handling of invalid inputs."""
        precip = np.random.gamma(2, 2, size=(60,))
        
        # Invalid distribution
        with pytest.raises(ValueError, match="only 'gamma' distribution is supported"):
            standardized_precipitation_index(precip, distribution='normal')
        
        # Invalid time scale
        with pytest.raises(ValueError, match="Time scale must be >= 1"):
            standardized_precipitation_index(precip, time_scale=0)


@pytest.mark.skipif(not HAS_XARRAY, reason="xarray not available")
class TestSPIXarray:
    """Test xarray interface for SPI calculations."""
    
    def test_spi_xarray_basic(self):
        """Test basic xarray SPI calculation."""
        # Create sample xarray DataArray
        np.random.seed(42)
        time = range(60)  # 5 years monthly
        lat = np.linspace(-30, 30, 5)
        lon = np.linspace(0, 360, 4, endpoint=False)
        
        precip_data = np.random.gamma(2, 2, size=(60, 5, 4))
        precip = xr.DataArray(
            precip_data,
            coords={'time': time, 'lat': lat, 'lon': lon},
            dims=['time', 'lat', 'lon'],
            attrs={'units': 'mm', 'long_name': 'precipitation'}
        )
        
        # Calculate SPI
        spi_result = spi_xarray(precip, time_scale=3)
        
        # Check output properties
        assert isinstance(spi_result, xr.DataArray)
        assert spi_result.shape == precip.shape
        assert spi_result.dims == precip.dims
        
        # Check coordinates are preserved
        for coord_name in ['time', 'lat', 'lon']:
            assert coord_name in spi_result.coords
            assert_array_almost_equal(spi_result.coords[coord_name], precip.coords[coord_name])
        
        # Check attributes
        assert spi_result.attrs['units'] == '1'
        assert 'Standardized Precipitation Index' in spi_result.attrs['long_name']
        assert spi_result.attrs['spi_time_scale'] == 3
    
    def test_spi_xarray_time_dim_detection(self):
        """Test automatic time dimension detection."""
        np.random.seed(42)
        
        # Test with different time dimension names
        for time_name in ['time', 'Time', 't']:
            coords = {time_name: range(36), 'space': range(5)}
            dims = [time_name, 'space']
            
            precip = xr.DataArray(
                np.random.gamma(2, 2, size=(36, 5)),
                coords=coords,
                dims=dims
            )
            
            spi_result = spi_xarray(precip, time_scale=1)
            assert spi_result.shape == precip.shape
    
    def test_spi_xarray_explicit_time_dim(self):
        """Test explicit time dimension specification."""
        np.random.seed(42)
        precip = xr.DataArray(
            np.random.gamma(2, 2, size=(5, 36)),  # space, time
            coords={'space': range(5), 'month': range(36)},
            dims=['space', 'month']
        )
        
        # Specify time dimension explicitly
        spi_result = spi_xarray(precip, time_scale=3, time_dim='month')
        assert spi_result.shape == precip.shape
    
    def test_spi_dataset(self):
        """Test SPI calculation for entire dataset with multiple time scales."""
        np.random.seed(42)
        time = range(60)
        lat = np.linspace(-30, 30, 3)
        lon = np.linspace(0, 360, 4, endpoint=False)
        
        # Create dataset with precipitation
        precip_data = np.random.gamma(2, 2, size=(60, 3, 4))
        dataset = xr.Dataset({
            'precipitation': xr.DataArray(
                precip_data,
                coords={'time': time, 'lat': lat, 'lon': lon},
                dims=['time', 'lat', 'lon']
            )
        })
        
        # Calculate SPI for multiple time scales
        spi_ds = spi_dataset(dataset, 'precipitation', time_scales=[1, 3, 6])
        
        # Check output
        assert isinstance(spi_ds, xr.Dataset)
        assert 'spi_1m' in spi_ds.data_vars
        assert 'spi_3m' in spi_ds.data_vars
        assert 'spi_6m' in spi_ds.data_vars
        
        # Check shapes
        for var_name in ['spi_1m', 'spi_3m', 'spi_6m']:
            assert spi_ds[var_name].shape == (60, 3, 4)
    
    def test_invalid_precipitation_variable(self):
        """Test error handling for invalid precipitation variable."""
        dataset = xr.Dataset({'temperature': xr.DataArray([1, 2, 3])})
        
        with pytest.raises(ValueError, match="Precipitation variable 'precip' not found"):
            spi_dataset(dataset, 'precip')
    
    def test_non_xarray_input(self):
        """Test error handling for non-xarray input."""
        with pytest.raises(TypeError, match="precipitation must be an xarray DataArray"):
            spi_xarray(np.array([1, 2, 3]))


class TestSPIIntegration:
    """Integration tests for SPI calculations."""
    
    def test_spi_known_values(self):
        """Test SPI calculation with known values for validation."""
        # Create a simple test case with known drought/wet patterns
        # Alternate between low and high precipitation
        precip = np.array([0.1, 5.0, 0.1, 5.0, 0.1, 5.0] * 10)  # 60 months
        
        spi_values = standardized_precipitation_index(precip, time_scale=1, axis=0)
        
        # Check that low precipitation periods have negative SPI
        # and high precipitation periods have positive SPI
        # (after sufficient data for distribution fitting)
        
        valid_indices = np.where(np.isfinite(spi_values))[0]
        if len(valid_indices) > 20:  # Need sufficient data
            # Low precip indices (even positions in pattern)
            low_precip_mask = np.array([(i % 6) in [0, 2, 4] for i in valid_indices])
            high_precip_mask = np.array([(i % 6) in [1, 3, 5] for i in valid_indices])
            
            if np.any(low_precip_mask) and np.any(high_precip_mask):
                low_spi = spi_values[valid_indices[low_precip_mask]]
                high_spi = spi_values[valid_indices[high_precip_mask]]
                
                # Most low precipitation should have negative SPI
                assert np.mean(low_spi < 0) > 0.6
                # Most high precipitation should have positive SPI  
                assert np.mean(high_spi > 0) > 0.6
    
    def test_spi_consistency_across_scales(self):
        """Test that SPI is consistent across different time scales."""
        np.random.seed(42)
        precip = np.random.gamma(2, 2, size=(120,))  # 10 years
        
        spi_1m = standardized_precipitation_index(precip, time_scale=1)
        spi_3m = standardized_precipitation_index(precip, time_scale=3)
        spi_12m = standardized_precipitation_index(precip, time_scale=12)
        
        # Longer time scales should be smoother (less variance)
        valid_1m = spi_1m[np.isfinite(spi_1m)]
        valid_3m = spi_3m[np.isfinite(spi_3m)]
        valid_12m = spi_12m[np.isfinite(spi_12m)]
        
        if len(valid_1m) > 30 and len(valid_3m) > 30 and len(valid_12m) > 30:
            var_1m = np.var(valid_1m)
            var_3m = np.var(valid_3m)
            var_12m = np.var(valid_12m)
            
            # Longer time scales should generally have less variance
            # (though this is not strictly guaranteed for all datasets)
            assert var_12m <= var_1m * 1.5  # Allow some flexibility
    
    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not available")
    def test_xarray_numpy_consistency(self):
        """Test that xarray and numpy interfaces give same results."""
        np.random.seed(42)
        precip_data = np.random.gamma(2, 2, size=(60, 3, 4))
        
        # Calculate with numpy interface
        spi_numpy = standardized_precipitation_index(precip_data, time_scale=3, axis=0)
        
        # Calculate with xarray interface
        precip_xr = xr.DataArray(
            precip_data,
            coords={'time': range(60), 'lat': range(3), 'lon': range(4)},
            dims=['time', 'lat', 'lon']
        )
        spi_xarray_result = spi_xarray(precip_xr, time_scale=3)
        
        # Results should be very close
        assert_allclose(spi_numpy, spi_xarray_result.values, rtol=1e-10)


if __name__ == "__main__":
    # Run basic tests if executed directly
    test_core = TestSPICore()
    test_core.test_rolling_sum_basic()
    test_core.test_spi_calculation_basic()
    print("Basic SPI tests passed!")
    
    if HAS_XARRAY:
        test_xr = TestSPIXarray()
        test_xr.test_spi_xarray_basic()
        print("Xarray SPI tests passed!")
    else:
        print("Xarray not available, skipping xarray tests")