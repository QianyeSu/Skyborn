"""
Mann-Kendall trend analysis for Skyborn package.

This module provides optimized implementations of the Mann-Kendall test
for trend detection in time series data, supporting both xarray and numpy arrays.

Examples
--------
Basic usage with 1D time series:

>>> import numpy as np
>>> from skyborn.calc import mann_kendall_test
>>>
>>> # Create time series with trend
>>> x = np.arange(100)
>>> y = 0.5 * x + np.random.randn(100)
>>>
>>> # Perform Mann-Kendall test
>>> result = mann_kendall_test(y)
>>> print(f"Trend: {result['trend']:.3f}")
>>> print(f"Significant: {result['h']}")

Multidimensional trend analysis:

>>> import numpy as np
>>> from skyborn.calc import trend_analysis
>>>
>>> # Create 3D data (time, lat, lon)
>>> data = np.random.randn(100, 50, 80) + np.linspace(0, 2, 100)[:, None, None]
>>>
>>> # Analyze trends along time axis (axis 0)
>>> results = trend_analysis(data, axis=0)
>>> print(f"Trend shape: {results['trend'].shape}")  # (50, 80)

Flexible dimension specification with xarray:

>>> import xarray as xr
>>> from skyborn.calc import trend_analysis
>>>
>>> # Create sample data
>>> time = pd.date_range('2000', periods=120, freq='M')
>>> data = xr.DataArray(
...     np.random.randn(120, 10, 15) + np.linspace(0, 3, 120)[:, None, None],
...     coords={'time': time, 'lat': range(10), 'lon': range(15)},
...     dims=['time', 'lat', 'lon']
... )
>>>
>>> # Analyze trends along any dimension
>>> results = trend_analysis(data, dim='time')
>>> print(results.trend.values.shape)  # (10, 15)
"""

import warnings
from importlib import import_module
from importlib import util as importlib_util
from importlib.machinery import EXTENSION_SUFFIXES
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union

import dask.array as da
import numpy as np
import scipy.stats as stats
import xarray as xr
from dask.diagnostics import ProgressBar

__all__ = [
    "mann_kendall_test",
    "mann_kendall_multidim",
    "mann_kendall_xarray",
    "trend_analysis",
    "mk_test",
    "mk_multidim",
]


def _load_mk_backend():
    """Best-effort loader for the compiled Mann-Kendall backend."""
    candidate_names = []
    if __package__:
        candidate_names.append(f"{__package__}.mann_kendall_core.mann_kendall_core")
    candidate_names.append("skyborn.calc.mann_kendall_core.mann_kendall_core")

    for module_name in candidate_names:
        try:
            return import_module(module_name)
        except Exception:
            continue

    backend_dir = Path(__file__).resolve().parent / "mann_kendall_core"
    for suffix in EXTENSION_SUFFIXES:
        candidate = backend_dir / f"mann_kendall_core{suffix}"
        if not candidate.exists():
            continue

        try:
            spec = importlib_util.spec_from_file_location(
                "mann_kendall_core", candidate
            )
            if spec is None or spec.loader is None:
                continue
            backend = importlib_util.module_from_spec(spec)
            spec.loader.exec_module(backend)
            return backend
        except Exception:
            continue

    return None


_mk_backend = _load_mk_backend()
_mk_score_var_batch_clean_backend = getattr(
    _mk_backend, "mk_score_var_batch_clean", None
)
_mk_score_var_sen_batch_clean_backend = getattr(
    _mk_backend, "mk_score_var_sen_batch_clean", None
)
_sen_slope_batch_clean_backend = getattr(_mk_backend, "sen_slope_batch_clean", None)


def _backend_available() -> bool:
    """Return whether all compiled clean-series kernels are available."""
    return (
        _mk_score_var_batch_clean_backend is not None
        and _mk_score_var_sen_batch_clean_backend is not None
        and _sen_slope_batch_clean_backend is not None
    )


def _as_fortran_float64_2d(data_2d: np.ndarray) -> np.ndarray:
    """Convert clean 2D data to the float64 Fortran layout expected by f2py."""
    array = np.asarray(data_2d, dtype=np.float64)
    if array.ndim != 2:
        raise ValueError("Expected a 2D array for Mann-Kendall clean-series kernels.")
    if array.flags.f_contiguous:
        return array
    return np.asfortranarray(array)


def _require_backend_function(function, function_name: str):
    """Return a compiled backend entry point or raise a clear import error."""
    if function is None:
        raise ImportError(
            "skyborn.calc.mann_kendall requires the compiled "
            f"mann_kendall_core entry '{function_name}'. "
            "Install a prebuilt wheel or build the Skyborn extensions first."
        )
    return function


def _mk_score_var_batch_clean(
    data_2d: np.ndarray, modified: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
    """Run the compiled clean 2D score/variance kernel."""
    clean_data = np.asarray(data_2d, dtype=np.float64)
    backend_function = _require_backend_function(
        _mk_score_var_batch_clean_backend, "mk_score_var_batch_clean"
    )
    s_values, var_values = backend_function(
        _as_fortran_float64_2d(clean_data), int(modified)
    )
    return np.asarray(s_values, dtype=np.float64), np.asarray(
        var_values, dtype=np.float64
    )


def _sen_slope_batch_clean(data_2d: np.ndarray) -> np.ndarray:
    """Run the compiled clean 2D Theil-Sen slope kernel."""
    clean_data = np.asarray(data_2d, dtype=np.float64)
    backend_function = _require_backend_function(
        _sen_slope_batch_clean_backend, "sen_slope_batch_clean"
    )
    slopes = backend_function(_as_fortran_float64_2d(clean_data))
    return np.asarray(slopes, dtype=np.float64)


def _mk_score_var_sen_batch_clean(
    data_2d: np.ndarray, modified: bool = False
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Run the compiled clean 2D combined score/variance/slope kernel."""
    clean_data = np.asarray(data_2d, dtype=np.float64)
    backend_function = _require_backend_function(
        _mk_score_var_sen_batch_clean_backend, "mk_score_var_sen_batch_clean"
    )
    s_values, var_values, slopes = backend_function(
        _as_fortran_float64_2d(clean_data), int(modified)
    )
    return (
        np.asarray(s_values, dtype=np.float64),
        np.asarray(var_values, dtype=np.float64),
        np.asarray(slopes, dtype=np.float64),
    )


def mann_kendall_test(
    data: Union[np.ndarray, "xr.DataArray"],
    alpha: float = 0.05,
    method: str = "theilslopes",
    modified: bool = False,
) -> Dict[str, Union[float, bool]]:
    """
    Perform Mann-Kendall test for trend detection on 1D time series.

    The Mann-Kendall test is a nonparametric test for detecting monotonic trends
    in time series data. It makes no assumptions about the distribution of the data.

    Parameters
    ----------
    data : array-like
        1D time series data. Can be numpy array or xarray DataArray.
    alpha : float, default 0.05
        Significance level for hypothesis testing (Type I error probability).
    method : str, default 'theilslopes'
        Method for calculating slope:
        - 'theilslopes': Theil-Sen slope estimator (robust, recommended)
        - 'linregress': Linear regression slope (faster but less robust)
    modified : bool, default False
        Use modified Mann-Kendall test (Yue and Wang, 2004) to account for
        serial autocorrelation in the data.

    Returns
    -------
    result : dict
        Dictionary containing trend analysis results:

        - 'trend' : float
            Slope of the trend (units per time step)
        - 'h' : bool
            Hypothesis test result. True if significant trend exists at
            the specified alpha level.
        - 'p' : float
            Two-tailed p-value of the test
        - 'z' : float
            Normalized test statistic (z-score)
        - 'tau' : float
            Kendall's tau correlation coefficient
        - 'std_error' : float
            Standard error of the detrended residuals

    Notes
    -----
    The Mann-Kendall test statistic S is calculated as:

    .. math::
        S = \\sum_{i=1}^{n-1} \\sum_{j=i+1}^{n} \\text{sign}(x_j - x_i)

    The test assumes that:
    - Data are independent (or modified=True for autocorrelated data)
    - Data come from the same distribution
    - Missing values are handled by removal

    References
    ----------
    Mann, H. B. (1945). Nonparametric tests against trend. Econometrica, 13(3), 245-259.
    Kendall, M. G. (1948). Rank correlation methods. Griffin, London.
    Yue, S., & Wang, C. (2004). The Mann-Kendall test modified by effective sample size
    to detect trend in serially correlated hydrological series. Water resources
    management, 18(3), 201-218.

    Examples
    --------
    >>> import numpy as np
    >>> from skyborn.calc import mann_kendall_test
    >>>
    >>> # Generate trend data
    >>> t = np.arange(50)
    >>> data = 0.3 * t + np.random.normal(0, 0.5, 50)
    >>>
    >>> # Test for trend
    >>> result = mann_kendall_test(data)
    >>> print(f"Trend slope: {result['trend']:.3f}")
    >>> print(f"P-value: {result['p']:.3f}")
    >>> print(f"Significant trend: {result['h']}")
    """
    if hasattr(data, "values"):
        values = data.values
    else:
        values = np.asarray(data)

    if values.ndim > 1:
        raise ValueError(
            "mann_kendall_test only accepts 1D data. Use trend_analysis for multidimensional data."
        )

    valid_mask = np.isfinite(values)
    y = np.asarray(values[valid_mask], dtype=np.float64)
    x = np.arange(len(y), dtype=np.float64)
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

    if method == "theilslopes":
        S_values, var_values, slopes = _mk_score_var_sen_batch_clean(
            y.reshape(-1, 1), modified=modified
        )
        S = int(np.rint(S_values[0]))
        var_s = float(var_values[0])
        slope = float(slopes[0])
        intercept = np.median(y) - slope * np.median(x)
    else:
        S = _calculate_mk_score(y)
        var_s = _calculate_mk_variance(y, n, modified)

    if S > 0:
        z = (S - 1) / np.sqrt(var_s)
    elif S == 0:
        z = 0
    else:
        z = (S + 1) / np.sqrt(var_s)

    p_value = 2 * (1 - stats.norm.cdf(abs(z)))
    h = abs(z) > stats.norm.ppf(1 - alpha / 2)

    if method == "linregress":
        slope, intercept = stats.linregress(x, y)[:2]
    elif method != "theilslopes":
        raise ValueError(f"Unknown method: {method}. Use 'theilslopes' or 'linregress'")

    tau = S / (0.5 * n * (n - 1))
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
    s_values, _ = _mk_score_var_batch_clean(
        np.asarray(y, dtype=np.float64).reshape(-1, 1)
    )
    return int(np.rint(s_values[0]))


def _calculate_mk_variance(y: np.ndarray, n: int, modified: bool = False) -> float:
    """Calculate variance of Mann-Kendall score."""
    _, variances = _mk_score_var_batch_clean(
        np.asarray(y, dtype=np.float64).reshape(-1, 1), modified=modified
    )
    return float(variances[0])


def mann_kendall_multidim(
    data: np.ndarray,
    axis: Union[int, str] = 0,
    alpha: float = 0.05,
    method: str = "theilslopes",
    modified: bool = False,
    chunk_size: Optional[int] = None,
    dim: Optional[Union[int, str]] = None,
) -> Dict[str, np.ndarray]:
    """
    Optimized numpy-based Mann-Kendall test for multidimensional arrays.

    This implementation is typically faster than xarray-based versions for large arrays
    because it avoids xarray overhead and uses optimized numpy operations.

    Supports various input dimensions:
    - 1D: (time,) - single time series
    - 2D: (time, ensemble) or (ensemble, time) - ensemble time series
    - 3D: (time, lat, lon) - spatial grid
    - ND: arbitrary multidimensional arrays

    Parameters
    ----------
    data : np.ndarray
        Input array with time series along one axis
    axis : int or str, default 0
        Axis along which to compute trends. Can be integer index or dimension name.
        For numpy arrays with string names, the array must have named axes.
    alpha : float, default 0.05
        Significance level
    method : str, default 'theilslopes'
        Slope calculation method ('theilslopes' or 'linregress')
    modified : bool, default False
        Use modified Mann-Kendall test for autocorrelated data
    chunk_size : int, optional
        Process data in chunks to manage memory usage
    dim : int or str, optional
        Alternative name for axis parameter (XArray style). Takes precedence over axis.

    Returns
    -------
    result : dict
        Dictionary with result arrays for each statistic:
        - 'trend': slope values
        - 'h': significance test results (boolean)
        - 'p': p-values
        - 'z': z-scores
        - 'tau': Kendall's tau values
        - 'std_error': standard errors

    Notes
    -----
    Parameter priority: dim > axis
    Both axis and dim can accept integer indices or string dimension names.
    For numpy arrays, string names require the array to have a 'dims' attribute or similar.
    """

    def _resolve_axis(arr: Any, axis_param: Union[int, str]) -> int:
        """Resolve axis parameter to integer index."""
        if isinstance(axis_param, str):
            # For numpy arrays, we need to check if it has dimension names
            if hasattr(arr, "dims"):
                # XArray-like object
                return arr.get_axis_num(axis_param)
            elif hasattr(arr, "axis_names") and arr.axis_names:
                # Custom numpy array with axis names
                return arr.axis_names.index(axis_param)
            else:
                # Plain numpy array - assume common dimension names
                common_time_names = [
                    "time",
                    "t",
                    "year",
                    "years",
                    "month",
                    "months",
                    "day",
                    "days",
                    "hour",
                    "hours",
                ]
                if axis_param.lower() in common_time_names:
                    # Assume time is the first dimension by default
                    return 0
                else:
                    raise ValueError(
                        f"Cannot resolve string axis '{axis_param}' for numpy array without dimension names. "
                        f"Use integer index or ensure array has 'dims' attribute."
                    )
        else:
            return axis_param

    # Handle parameter aliases with priority: dim > axis
    if dim is not None:
        actual_time_axis = _resolve_axis(data, dim)
    else:
        actual_time_axis = _resolve_axis(data, axis)
    # Handle 1D input
    if data.ndim == 1:
        return mann_kendall_test(data, alpha=alpha, method=method, modified=modified)

    # Get original shape for result reshaping
    original_shape = data.shape

    # Move time axis to the front
    data = np.moveaxis(data, actual_time_axis, 0)
    time_steps = data.shape[0]

    # Handle different input shapes
    if data.ndim == 2:
        # Already in the right format: (time, n_series)
        data_2d = data
        spatial_shape = (data.shape[1],)
    else:
        # Reshape to 2D: (time, space)
        spatial_shape = data.shape[1:]
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
        # Estimate memory usage and set reasonable chunk size
        max_memory_mb = 200  # 200MB limit
        bytes_per_element = 8
        max_chunk_size = (max_memory_mb * 1024 * 1024) // (
            time_steps * time_steps * bytes_per_element
        )
        chunk_size = min(max_chunk_size, n_points, 10000)
        chunk_size = max(chunk_size, 100)  # Minimum chunk size

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
    method: str = "theilslopes",
    modified: bool = False,
) -> Dict[str, np.ndarray]:
    """
    Truly vectorized Mann-Kendall test for a chunk of time series.

    Uses optimized numpy operations to process multiple time series simultaneously.
    This is a major performance improvement over the previous loop-based implementation.

    Handles missing data (NaN) gracefully by:
    1. Full vectorization for series without NaN
    2. Individual processing for series with NaN
    3. Automatic skipping of series with insufficient data
    """
    time_steps, n_series = data_chunk.shape
    x = np.arange(time_steps, dtype=np.float64)

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

    # Get valid indices for efficient assignment
    valid_indices = np.where(valid_series)[0]

    # Separate series with complete data from those with missing values
    no_nan_mask = valid_counts[valid_series] == time_steps
    no_nan_indices = valid_indices[no_nan_mask]

    # Process series without NaN using full vectorization.
    # NumPy returns a scalar directly; Dask returns a deferred scalar.
    n_clean_series = no_nan_mask.sum()
    if hasattr(n_clean_series, "compute"):
        n_clean_series = n_clean_series.compute()

    if n_clean_series > 0:
        if hasattr(data_chunk, "compute"):
            clean_data_computed = data_chunk[:, no_nan_indices].compute()
        else:
            clean_data_computed = np.asarray(data_chunk[:, no_nan_indices])

        clean_data_computed = np.asarray(clean_data_computed, dtype=np.float64)
        if method == "theilslopes":
            S_values, var_s_values, slopes = _mk_score_var_sen_batch_clean(
                clean_data_computed, modified=modified
            )
        else:
            S_values, var_s_values = _mk_score_var_batch_clean(
                clean_data_computed, modified=modified
            )

        # Vectorized z-score calculation
        z_values = np.zeros_like(S_values, dtype=np.float64)
        positive_mask = S_values > 0
        negative_mask = S_values < 0
        z_values[positive_mask] = (S_values[positive_mask] - 1) / np.sqrt(
            var_s_values[positive_mask]
        )
        z_values[negative_mask] = (S_values[negative_mask] + 1) / np.sqrt(
            var_s_values[negative_mask]
        )

        p_values = 2 * (1 - stats.norm.cdf(np.abs(z_values)))
        h_values = np.abs(z_values) > stats.norm.ppf(1 - alpha / 2)

        if method != "theilslopes":
            slopes = _vectorized_linregress_slopes(clean_data_computed, x)

        n = time_steps
        tau_values = S_values / (0.5 * n * (n - 1))

        if method == "theilslopes":
            std_errors = _vectorized_std_error_theil(clean_data_computed, x, slopes)
        else:
            std_errors = _vectorized_std_error_linregress(
                clean_data_computed, x, slopes
            )

        results["trend"][no_nan_indices] = slopes
        results["h"][no_nan_indices] = h_values
        results["p"][no_nan_indices] = p_values
        results["z"][no_nan_indices] = z_values
        results["tau"][no_nan_indices] = tau_values
        results["std_error"][no_nan_indices] = std_errors

    # Handle series with NaN individually (but batch them for efficiency)
    nan_mask = valid_counts[valid_indices] != time_steps
    nan_indices = valid_indices[nan_mask]

    # Convert to regular numpy if it's a Dask array to avoid indexing issues
    if hasattr(data_chunk, "compute"):
        data_for_nan = data_chunk.compute()
        # Also convert indices to regular numpy for iteration
        if hasattr(nan_indices, "compute"):
            nan_indices = nan_indices.compute()
    else:
        data_for_nan = data_chunk

    for series_idx in nan_indices:
        y = data_for_nan[:, series_idx]
        finite_mask = np.isfinite(y)

        y_clean = y[finite_mask]

        # Mann-Kendall test for individual series with NaN
        try:
            mk_result = mann_kendall_test(
                y_clean, alpha=alpha, method=method, modified=modified
            )

            # Store results
            for key, value in mk_result.items():
                results[key][series_idx] = value
        except Exception:
            # Skip problematic series - results remain NaN
            continue

    return results


def _vectorized_linregress_slopes(data_2d: np.ndarray, x: np.ndarray) -> np.ndarray:
    """
    Vectorized linear regression slope calculation.

    Much faster than calling stats.linregress for each series individually.
    """
    time_steps, n_series = data_2d.shape

    # Center the data
    x_mean = np.mean(x)
    y_mean = np.mean(data_2d, axis=0)  # Shape: (n_series,)

    # Calculate numerator and denominator for slope
    x_centered = x - x_mean  # Shape: (time_steps,)
    y_centered = data_2d - y_mean  # Shape: (time_steps, n_series)

    # Vectorized slope calculation: sum(xy) / sum(x²)
    numerator = np.sum(
        x_centered[:, np.newaxis] * y_centered, axis=0
    )  # Shape: (n_series,)
    denominator = np.sum(x_centered**2)  # Scalar

    slopes = numerator / denominator

    return slopes


def _vectorized_std_error_theil(
    data_2d: np.ndarray, x: np.ndarray, slopes: np.ndarray
) -> np.ndarray:
    """
    Vectorized standard error calculation for Theil-Sen slopes.
    """
    time_steps, n_series = data_2d.shape

    # Calculate median intercepts
    intercepts = np.median(data_2d - x[:, np.newaxis] * slopes, axis=0)

    # Calculate residuals
    y_pred = x[:, np.newaxis] * slopes + intercepts  # Shape: (time_steps, n_series)
    residuals = data_2d - y_pred

    # Standard error estimation
    std_errors = np.std(residuals, axis=0) / np.sqrt(time_steps)

    return std_errors


def _vectorized_std_error_linregress(
    data_2d: np.ndarray, x: np.ndarray, slopes: np.ndarray
) -> np.ndarray:
    """
    Vectorized residual standard error calculation for linear regression.

    This mirrors the scalar `mann_kendall_test` implementation, which returns
    the standard error of the detrended residuals rather than the regression
    slope standard error.
    """
    time_steps, n_series = data_2d.shape

    # Calculate intercepts
    x_mean = np.mean(x)
    y_mean = np.mean(data_2d, axis=0)
    intercepts = y_mean - slopes * x_mean

    # Calculate residuals
    y_pred = x[:, np.newaxis] * slopes + intercepts
    residuals = data_2d - y_pred

    # Match the scalar API: residual spread scaled by sqrt(n)
    std_errors = np.std(residuals, axis=0) / np.sqrt(time_steps)

    return std_errors


def mann_kendall_xarray(
    data,  # xr.DataArray
    dim: str = "time",
    alpha: float = 0.05,
    method: str = "theilslopes",
    modified: bool = False,
    use_dask: bool = True,
):  # -> xr.Dataset
    """
    Mann-Kendall test for xarray DataArray with intelligent dimension handling.

    Automatically handles:
    - 1D time series (pure time dimension)
    - Multi-dimensional data with preserved dimension order
    - Scalar dimensions (size=1) from sel() operations
    - Correct dimension ordering in output

    Parameters
    ----------
    data : xr.DataArray
        Input data array
    dim : str, default 'time'
        Time dimension name to analyze along
    alpha : float, default 0.05
        Significance level
    method : str, default 'theilslopes'
        Slope calculation method
    modified : bool, default False
        Use modified Mann-Kendall test
    use_dask : bool, default True
        Use dask for computation if available

    Returns
    -------
    result : xr.Dataset
        Dataset containing trend analysis results with preserved dimension order

    Examples
    --------
    >>> # 3D data (time, lat, lon) -> (lat, lon)
    >>> ds = mann_kendall_xarray(data_3d, dim='time')
    >>>
    >>> # 3D data (year, ensemble, lat) -> (ensemble, lat)
    >>> ds = mann_kendall_xarray(data_3d, dim='year')
    >>>
    >>> # Handle scalar dimension from sel
    >>> subset = data.sel(lat=0)  # lat becomes scalar
    >>> ds = mann_kendall_xarray(subset, dim='time')  # Works correctly
    """
    # Get time axis
    time_axis = data.get_axis_num(dim)

    if use_dask and hasattr(data.data, "chunks"):
        # Use dask for computation
        results = _dask_mann_kendall(data, dim, alpha, method, modified)
    else:
        # Use numpy implementation
        results = mann_kendall_multidim(
            data.values,
            axis=time_axis,
            alpha=alpha,
            method=method,
            modified=modified,
        )

    # Smart dimension handling: preserve order and filter scalar dimensions
    # Get all dimensions except the analyzed dimension, in original order
    output_dims = []
    output_coords = {}

    for d in data.dims:
        if d == dim:
            # Skip the analyzed dimension
            continue

        # Check if this is a true dimension (size > 1) or scalar (size == 1)
        if data.sizes[d] > 1:
            output_dims.append(d)
            output_coords[d] = data.coords[d]
        # Scalar dimensions (size == 1) are excluded from output

    # Handle different cases based on output dimensions
    if not output_dims:
        # Pure scalar output: no non-time dimensions or all are size=1
        # Return scalar dataset
        ds = xr.Dataset(
            data_vars={
                "trend": (
                    [],
                    (
                        results["trend"]
                        if np.isscalar(results["trend"])
                        else results["trend"].item()
                    ),
                ),
                "h": (
                    [],
                    results["h"] if np.isscalar(results["h"]) else results["h"].item(),
                ),
                "p": (
                    [],
                    results["p"] if np.isscalar(results["p"]) else results["p"].item(),
                ),
                "z": (
                    [],
                    results["z"] if np.isscalar(results["z"]) else results["z"].item(),
                ),
                "tau": (
                    [],
                    (
                        results["tau"]
                        if np.isscalar(results["tau"])
                        else results["tau"].item()
                    ),
                ),
                "std_error": (
                    [],
                    (
                        results["std_error"]
                        if np.isscalar(results["std_error"])
                        else results["std_error"].item()
                    ),
                ),
            },
            coords={},
            attrs={
                "title": "Mann-Kendall Trend Analysis",
                "alpha": alpha,
                "method": method,
                "modified": modified,
                "input_dims": str(data.dims),
                "analyzed_dim": dim,
            },
        )
    else:
        # Multi-dimensional output: preserve original dimension order
        ds = xr.Dataset(
            data_vars={
                "trend": (output_dims, results["trend"]),
                "h": (output_dims, results["h"]),
                "p": (output_dims, results["p"]),
                "z": (output_dims, results["z"]),
                "tau": (output_dims, results["tau"]),
                "std_error": (output_dims, results["std_error"]),
            },
            coords=output_coords,
            attrs={
                "title": "Mann-Kendall Trend Analysis",
                "alpha": alpha,
                "method": method,
                "modified": modified,
                "input_dims": str(data.dims),
                "analyzed_dim": dim,
            },
        )

    return ds


def _dask_mann_kendall(
    data, dim: str, alpha: float, method: str, modified: bool  # xr.DataArray
):  # -> Dict[str, np.ndarray]
    """Use dask map_blocks for Mann-Kendall computation."""

    # Get time axis
    time_axis = data.get_axis_num(dim)

    # Convert to dask array if not already
    if hasattr(data, "data") and isinstance(data.data, da.Array):
        dask_data = data.data
    else:
        # Create dask array with reasonable chunking
        chunks = list(data.shape)
        chunks[time_axis] = -1  # Keep time dimension unchunked

        # Chunk spatial dimensions reasonably
        for i, size in enumerate(chunks):
            if i != time_axis and size > 100:
                chunks[i] = min(50, size // 2)

        dask_data = da.from_array(data.values, chunks=chunks)

    def mk_block_processor(block):
        """Process a block of data with Mann-Kendall analysis."""
        # Move time axis to front
        block = np.moveaxis(block, time_axis, 0)
        time_steps = block.shape[0]

        # Get spatial shape and flatten
        spatial_shape = block.shape[1:]
        spatial_size = np.prod(spatial_shape)
        block_2d = block.reshape(time_steps, spatial_size)

        # Initialize output arrays
        output_shape = (6,) + spatial_shape  # 6 output variables
        output = np.full(output_shape, np.nan)
        output_2d = output.reshape(6, spatial_size)

        # Process each spatial point
        for i in range(spatial_size):
            series = block_2d[:, i]

            # Skip if insufficient data
            if len(series) < 3 or np.sum(np.isfinite(series)) < 3:
                continue

            try:
                # Perform Mann-Kendall test
                result = mann_kendall_test(
                    series, alpha=alpha, method=method, modified=modified
                )

                # Store results: [trend, h, p, z, tau, std_error]
                output_2d[0, i] = result["trend"]
                output_2d[1, i] = float(result["h"])
                output_2d[2, i] = result["p"]
                output_2d[3, i] = result["z"]
                output_2d[4, i] = result["tau"]
                output_2d[5, i] = result["std_error"]

            except Exception:
                continue  # Skip problematic series

        return output.reshape(output_shape)

    # Determine output chunks
    output_chunks = list(dask_data.chunks)
    output_chunks[time_axis] = (6,)  # Replace time with 6 output variables

    # Apply Mann-Kendall analysis using map_blocks
    with ProgressBar():
        result_array = da.map_blocks(
            mk_block_processor,
            dask_data,
            dtype=float,
            chunks=output_chunks,
            drop_axis=time_axis,
            new_axis=time_axis,
        ).compute()

    # Move result axis to front and extract variables
    result_array = np.moveaxis(result_array, time_axis, 0)

    # Create result dictionary
    result_dict = {
        "trend": result_array[0],
        "h": result_array[1].astype(bool),
        "p": result_array[2],
        "z": result_array[3],
        "tau": result_array[4],
        "std_error": result_array[5],
    }

    return result_dict


# Convenience function for easy access
def trend_analysis(
    data: Any,
    axis: Union[int, str] = 0,
    alpha: float = 0.05,
    method: str = "theilslopes",
    modified: bool = False,
    dim: Optional[Union[int, str]] = None,
    **kwargs,
) -> Any:
    """
    Unified interface for Mann-Kendall trend analysis.

    Automatically chooses the best implementation based on input type.

    Parameters
    ----------
    data : np.ndarray or xr.DataArray
        Input data array
    axis : int or str, default 0
        Axis along which to compute trends. Can be integer index or dimension name.
    alpha : float, default 0.05
        Significance level
    method : str, default 'theilslopes'
        Slope calculation method
    modified : bool, default False
        Use modified Mann-Kendall test
    dim : int or str, optional
        Alternative name for axis parameter (XArray style). Takes precedence over axis.
    **kwargs
        Additional arguments passed to underlying functions

    Notes
    -----
    Parameter priority: dim > axis
    Both axis and dim support integer indices and string dimension names.
    """
    # Handle parameter aliases with priority: dim > axis
    if dim is not None:
        actual_axis = dim
    else:
        actual_axis = axis

    if isinstance(data, xr.DataArray):
        return mann_kendall_xarray(
            data,
            dim=actual_axis,
            alpha=alpha,
            method=method,
            modified=modified,
            **kwargs,
        )
    else:
        return mann_kendall_multidim(
            data,
            axis=actual_axis,
            alpha=alpha,
            method=method,
            modified=modified,
            **kwargs,
        )


# Aliases for backward compatibility and convenience
mk_test = mann_kendall_test  # Short alias
mk_multidim = mann_kendall_multidim  # Short alias for multidimensional
