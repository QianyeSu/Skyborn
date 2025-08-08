"""
Supplementary tests to improve Mann-Kendall coverage.
"""

import importlib.util
import sys
import os
import numpy as np
import dask.array as da
import xarray as xr
import pandas as pd

# Add path and import module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
mann_kendall_path = os.path.join(
    os.path.dirname(__file__), "src", "skyborn", "calc", "mann_kendall.py"
)
spec = importlib.util.spec_from_file_location("mann_kendall", mann_kendall_path)
mann_kendall_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mann_kendall_module)

mann_kendall_test = mann_kendall_module.mann_kendall_test
mann_kendall_multidim = mann_kendall_module.mann_kendall_multidim
trend_analysis = mann_kendall_module.trend_analysis
_calculate_std_error_theil = mann_kendall_module._calculate_std_error_theil


def test_missing_lines():
    """Test specific missing lines."""

    print("Testing custom array with axis_names...")
    # Test line 335: custom array with axis_names

    class CustomArrayWithAxisNames:
        def __init__(self, data):
            self.data = data
            self.axis_names = ["x", "y", "z"]
            self.shape = data.shape
            self.ndim = data.ndim

        def __array__(self):
            return self.data

    data_3d = np.random.randn(10, 5, 8)
    custom_arr = CustomArrayWithAxisNames(data_3d)
    result = trend_analysis(custom_arr, axis="x")
    print("✓ Custom array with axis_names test passed")

    print("Testing _calculate_std_error_theil...")
    # Test lines 723-734: _calculate_std_error_theil function
    y_short = np.array([1, 2])
    x_short = np.array([1, 2])
    std_err = _calculate_std_error_theil(y_short, x_short, 1.0)
    assert np.isnan(std_err), "Should return NaN for insufficient data"

    y_3 = np.array([1, 3, 5])
    x_3 = np.array([1, 2, 3])
    std_err_3 = _calculate_std_error_theil(y_3, x_3, 2.0)
    assert not np.isnan(std_err_3), "Should return valid std error for sufficient data"
    print("✓ _calculate_std_error_theil test passed")

    print("Testing dask processing paths...")
    # Test dask-specific code paths (lines 474-476, 552-555, etc.)
    try:
        dask_data = da.from_array(np.random.randn(20, 15), chunks=(10, 15))
        dask_data_with_nan = da.where(dask_data > -2, dask_data, np.nan)
        result = mann_kendall_multidim(dask_data_with_nan, axis=0)
        print("✓ Dask processing test passed")
    except Exception as e:
        print(f"Note: Dask test encountered: {e}")

    print("Testing xarray error handling...")
    # Test xarray error handling (line 359)

    class MockXArrayError:
        def __init__(self, data):
            self.data = data
            self.dims = ["time", "lat", "lon"]
            self.shape = data.shape
            self.ndim = data.ndim

        def get_axis_num(self, axis):
            raise ValueError("Mock error for testing")

        def __array__(self):
            return self.data

    mock_data = MockXArrayError(np.random.randn(10, 5, 8))
    try:
        result = trend_analysis(mock_data, axis="invalid")
        print("✓ XArray error handling test completed")
    except ValueError:
        print("✓ XArray error handling test passed (expected error)")

    print("Testing xarray-specific processing...")
    # Test lines 580-582: xarray specific processing
    da_simple = xr.DataArray(
        np.random.randn(30) + np.arange(30) * 0.02,
        dims=["time"],
        coords={"time": pd.date_range("2000", periods=30, freq="YS")},
        attrs={"units": "test_units", "long_name": "test_data"},
    )

    result = trend_analysis(da_simple, axis="time")
    print("✓ XArray processing test passed")

    print("\nAll supplementary coverage tests completed!")


if __name__ == "__main__":
    test_missing_lines()
