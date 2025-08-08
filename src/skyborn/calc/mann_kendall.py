"""
Mann-Kendall trend analysis for Skyborn package.

This module provides optimized implementations of the Mann-Kendall test
for trend detection in time series data, supporting both xarray and numpy arrays.
"""

import numpy as np
import scipy.stats as stats
from typing import Tuple, Union, Optional, Dict, Any
import warnings

try:
    import xarray as xr

    HAS_XARRAY = True
except ImportError:
    HAS_XARRAY = False
    xr = None

try:
    import dask.array as da
    from dask.diagnostics import ProgressBar

    HAS_DASK = True
except ImportError:
    HAS_DASK = False


def mann_kendall_test(
    data: Union[np.ndarray, xr.DataArray],
    alpha: float = 0.05,
    method: str = "auto",
    modified: bool = False,
) -> Dict[str, Union[float, bool]]:
    """
    Perform Mann-Kendall test for trend detection on 1D time series.

    Parameters
    ----------
    data : array-like
        1D time series data
    alpha : float, default 0.05
        Significance level for hypothesis testing
    method : str, default 'auto'
        Method for calculating slope ('theilslopes', 'linregress', 'auto')
    modified : bool, default False
        Use modified Mann-Kendall test (Yue and Wang, 2004) to account for autocorrelation

    Returns
    -------
    result : dict
        Dictionary containing:
        - 'trend': slope of the trend
        - 'h': hypothesis test result (True if significant trend exists)
        - 'p': p-value
        - 'z': normalized test statistic
        - 'tau': Kendall's tau
        - 'std_error': standard error of the trend
    """
    # Handle both numpy and xarray inputs
    if hasattr(data, "values"):
        values = data.values
    else:
        values = np.asarray(data)

    # Remove NaN values
    valid_mask = np.isfinite(values)
    y = values[valid_mask]
    x = np.arange(len(values))[valid_mask]
    n = len(y)

    if n < 3:
        warnings.warn("Need at least 3 data points for Mann-Kendall test")
        return {
            "trend": np.nan,
            "h": False,
            "p": np.nan,
            "z": np.nan,
            "tau": np.nan,
            "std_error": np.nan,
        }

    # Calculate Mann-Kendall score (S)
    S = _calculate_mk_score(y)

    # Calculate variance of S
    var_s = _calculate_mk_variance(y, n, modified)

    # Calculate normalized test statistic
    if S > 0:
        z = (S - 1) / np.sqrt(var_s)
    elif S == 0:
        z = 0
    else:
        z = (S + 1) / np.sqrt(var_s)

    # Calculate p-value and hypothesis test
    p_value = 2 * (1 - stats.norm.cdf(abs(z)))
    h = abs(z) > stats.norm.ppf(1 - alpha / 2)

    # Calculate slope
    if method == "auto":
        method = "theilslopes" if n > 100 else "linregress"

    if method == "theilslopes":
        slope, intercept = stats.mstats.theilslopes(y, x)[:2]
    elif method == "linregress":
        slope, intercept = stats.linregress(x, y)[:2]
    else:
        raise ValueError(f"Unknown method: {method}")

    # Calculate Kendall's tau
    tau = S / (0.5 * n * (n - 1))

    # Calculate standard error
    y_detrend = y - (x * slope + intercept)
    std_error = np.std(y_detrend) / np.sqrt(n)

    return {
        "trend": slope,
        "h": h,
        "p": p_value,
        "z": z,
        "tau": tau,
        "std_error": std_error,
    }


def _calculate_mk_score(y: np.ndarray) -> int:
    """Calculate Mann-Kendall score S."""
    n = len(y)
    S = 0

    for i in range(n - 1):
        for j in range(i + 1, n):
            S += np.sign(y[j] - y[i])

    return S


def _calculate_mk_variance(y: np.ndarray, n: int, modified: bool = False) -> float:
    """Calculate variance of Mann-Kendall score."""
    # Count tied groups
    unique_vals, counts = np.unique(y, return_counts=True)
    tied_groups = counts[counts > 1]

    # Basic variance calculation
    if len(tied_groups) == 0:
        var_s = n * (n - 1) * (2 * n + 5) / 18
    else:
        var_s = (
            n * (n - 1) * (2 * n + 5)
            - np.sum(tied_groups * (tied_groups - 1) * (2 * tied_groups + 5))
        ) / 18

    # Apply modification for autocorrelation if requested
    if modified:
        # Yue and Wang (2004) modification
        x = np.arange(n)
        slope, intercept = stats.linregress(x, y)[:2]
        y_detrend = y - (x * slope + intercept)

        # Calculate autocorrelation
        autocorr = _calculate_autocorrelation(y_detrend)

        # Calculate effective sample size multiplier
        n_star = 1 + 2 * np.sum([(1 - k / n) * autocorr[k] for k in range(1, n)])
        var_s *= n_star

    return var_s


def _calculate_autocorrelation(y: np.ndarray) -> np.ndarray:
    """Calculate autocorrelation function."""
    n = len(y)
    y_centered = y - np.mean(y)

    # Full autocorrelation using numpy correlate
    full_corr = np.correlate(y_centered, y_centered, mode="full")
    autocorr = full_corr[n - 1 :] / full_corr[n - 1]  # Normalize

    return autocorr


def mann_kendall_multidim_numpy(
    data: np.ndarray,
    time_axis: int = 0,
    alpha: float = 0.05,
    method: str = "auto",
    modified: bool = False,
    chunk_size: Optional[int] = None,
) -> Dict[str, np.ndarray]:
    """
    Optimized numpy-based Mann-Kendall test for multidimensional arrays.

    This implementation is typically faster than xarray-based versions for large arrays
    because it avoids xarray overhead and uses optimized numpy operations.

    Parameters
    ----------
    data : np.ndarray
        Input array with time series along one axis
    time_axis : int, default 0
        Axis along which to compute trends
    alpha : float, default 0.05
        Significance level
    method : str, default 'auto'
        Slope calculation method
    modified : bool, default False
        Use modified Mann-Kendall test
    chunk_size : int, optional
        Process data in chunks to manage memory usage

    Returns
    -------
    result : dict
        Dictionary with result arrays for each statistic
    """
    # Move time axis to the front
    data = np.moveaxis(data, time_axis, 0)
    time_steps, *spatial_shape = data.shape

    # Reshape to 2D: (time, space)
    data_2d = data.reshape(time_steps, -1)
    n_points = data_2d.shape[1]

    # Initialize output arrays
    results = {
        "trend": np.full(spatial_shape, np.nan),
        "h": np.zeros(spatial_shape, dtype=bool),
        "p": np.full(spatial_shape, np.nan),
        "z": np.full(spatial_shape, np.nan),
        "tau": np.full(spatial_shape, np.nan),
        "std_error": np.full(spatial_shape, np.nan),
    }

    # Flatten result arrays for efficient assignment
    results_flat = {key: val.ravel() for key, val in results.items()}

    # Determine chunk size
    if chunk_size is None:
        chunk_size = min(10000, n_points)  # Process up to 10k points at once

    # Process in chunks
    for start_idx in range(0, n_points, chunk_size):
        end_idx = min(start_idx + chunk_size, n_points)
        chunk_data = data_2d[:, start_idx:end_idx]

        # Vectorized Mann-Kendall calculation for chunk
        chunk_results = _vectorized_mk_test(
            chunk_data, alpha=alpha, method=method, modified=modified
        )

        # Store results
        for key, values in chunk_results.items():
            results_flat[key][start_idx:end_idx] = values

    # Reshape results back to original spatial shape
    for key in results:
        results[key] = results_flat[key].reshape(spatial_shape)

    return results


def _vectorized_mk_test(
    data_chunk: np.ndarray,
    alpha: float = 0.05,
    method: str = "auto",
    modified: bool = False,
) -> Dict[str, np.ndarray]:
    """
    Vectorized Mann-Kendall test for a chunk of time series.

    Uses optimized numpy operations to process multiple time series simultaneously.
    """
    time_steps, n_series = data_chunk.shape
    x = np.arange(time_steps)

    # Initialize result arrays
    results = {
        "trend": np.full(n_series, np.nan),
        "h": np.zeros(n_series, dtype=bool),
        "p": np.full(n_series, np.nan),
        "z": np.full(n_series, np.nan),
        "tau": np.full(n_series, np.nan),
        "std_error": np.full(n_series, np.nan),
    }

    # Find valid time series (with enough non-NaN values)
    valid_counts = np.sum(np.isfinite(data_chunk), axis=0)
    valid_series = valid_counts >= 3

    if not np.any(valid_series):
        return results

    # Process valid series
    valid_data = data_chunk[:, valid_series]
    n_valid = np.sum(valid_series)

    # Calculate Mann-Kendall scores vectorized
    S_values = np.zeros(n_valid)
    for i in range(n_valid):
        y = valid_data[:, i]
        finite_mask = np.isfinite(y)
        if np.sum(finite_mask) < 3:
            continue
        y_clean = y[finite_mask]
        S_values[i] = _calculate_mk_score(y_clean)

    # Calculate slopes and other statistics
    for i, series_idx in enumerate(np.where(valid_series)[0]):
        y = data_chunk[:, series_idx]
        finite_mask = np.isfinite(y)
        n_finite = np.sum(finite_mask)

        if n_finite < 3:
            continue

        y_clean = y[finite_mask]
        x_clean = x[finite_mask]

        # Mann-Kendall test
        mk_result = mann_kendall_test(
            y_clean, alpha=alpha, method=method, modified=modified
        )

        # Store results
        for key, value in mk_result.items():
            results[key][series_idx] = value

    return results


if HAS_XARRAY:

    def mann_kendall_xarray(
        data: xr.DataArray,
        dim: str = "time",
        alpha: float = 0.05,
        method: str = "auto",
        modified: bool = False,
        use_dask: bool = True,
    ) -> xr.Dataset:
        """
        Mann-Kendall test for xarray DataArray.

        Parameters
        ----------
        data : xr.DataArray
            Input data array
        dim : str, default 'time'
            Time dimension name
        alpha : float, default 0.05
            Significance level
        method : str, default 'auto'
            Slope calculation method
        modified : bool, default False
            Use modified Mann-Kendall test
        use_dask : bool, default True
            Use dask for computation if available

        Returns
        -------
        result : xr.Dataset
            Dataset containing trend analysis results
        """
        # Get time axis
        time_axis = data.get_axis_num(dim)

        if use_dask and HAS_DASK and hasattr(data.data, "chunks"):
            # Use dask for computation
            results = _dask_mann_kendall(data, dim, alpha, method, modified)
        else:
            # Use numpy implementation
            results = mann_kendall_multidim_numpy(
                data.values,
                time_axis=time_axis,
                alpha=alpha,
                method=method,
                modified=modified,
            )

        # Create output coordinates (all dims except time)
        coords = {k: v for k, v in data.coords.items() if k != dim}

        # Create Dataset
        ds = xr.Dataset(
            data_vars={
                "trend": (list(coords.keys()), results["trend"]),
                "h": (list(coords.keys()), results["h"]),
                "p": (list(coords.keys()), results["p"]),
                "z": (list(coords.keys()), results["z"]),
                "tau": (list(coords.keys()), results["tau"]),
                "std_error": (list(coords.keys()), results["std_error"]),
            },
            coords=coords,
            attrs={
                "title": "Mann-Kendall Trend Analysis",
                "alpha": alpha,
                "method": method,
                "modified": modified,
            },
        )

        return ds

    def _dask_mann_kendall(
        data: xr.DataArray, dim: str, alpha: float, method: str, modified: bool
    ) -> Dict[str, np.ndarray]:
        """Use dask for Mann-Kendall computation."""
        import dask.array as da

        # Apply along time axis using dask
        def mk_wrapper(arr, axis=0):
            return mann_kendall_multidim_numpy(
                arr, time_axis=axis, alpha=alpha, method=method, modified=modified
            )

        # Get time axis
        time_axis = data.get_axis_num(dim)

        # Apply function
        with ProgressBar():
            result_dict = da.apply_along_axis(
                mk_wrapper,
                time_axis,
                data.data,
                dtype=object,
                shape=data.shape[:time_axis] + data.shape[time_axis + 1 :],
            ).compute()

        return result_dict


# Convenience function for easy access
def trend_analysis(
    data: Union[np.ndarray, xr.DataArray],
    time_axis: Union[int, str] = 0,
    alpha: float = 0.05,
    method: str = "auto",
    modified: bool = False,
    **kwargs,
) -> Union[Dict[str, np.ndarray], xr.Dataset]:
    """
    Unified interface for Mann-Kendall trend analysis.

    Automatically chooses the best implementation based on input type.
    """
    if HAS_XARRAY and isinstance(data, xr.DataArray):
        return mann_kendall_xarray(
            data, dim=time_axis, alpha=alpha, method=method, modified=modified, **kwargs
        )
    else:
        if isinstance(time_axis, str):
            raise ValueError("String time_axis only supported for xarray DataArray")
        return mann_kendall_multidim_numpy(
            data,
            time_axis=time_axis,
            alpha=alpha,
            method=method,
            modified=modified,
            **kwargs,
        )
