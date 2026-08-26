"""
Lanczos filter implementation for multi-dimensional data.

This module provides Lanczos filtering for 1D, 2D, and N-D arrays with
optional Fortran acceleration for high performance.
"""

from __future__ import annotations

import warnings
from typing import Literal, Sequence

import numpy as np

try:
    import xarray as xr

    HAS_XARRAY = True
except ImportError:
    HAS_XARRAY = False
    xr = None

# Import Fortran backend (required)
from . import _lanczos_core

__all__ = [
    "lanczos_filter",
    "lanczos_lowpass",
    "lanczos_highpass",
    "lanczos_bandpass",
    "lanczos_weights",
]


def lanczos_weights(
    cutoff_freq: float,
    window: int,
    pass_type: Literal["low", "high"] = "low",
) -> np.ndarray:
    """
    Compute Lanczos filter weights (coefficients).

    Parameters
    ----------
    cutoff_freq : float
        Cutoff frequency in normalized units (0 to 0.5).
        For example, for daily data with 30-day cutoff: cutoff_freq = 1/30
    window : int
        Number of filter weights (must be odd).
        Larger window gives sharper cutoff but more edge effects.
        Typical values: 61, 81, 121
    pass_type : {'low', 'high'}, default='low'
        Type of filter:
        - 'low': low-pass (smoothing, removes high frequencies)
        - 'high': high-pass (removes trends, keeps high frequencies)

    Returns
    -------
    weights : ndarray
        Filter coefficients of length `window`, normalized to sum to 1 (for low-pass)
        or sum to 0 (for high-pass)

    Examples
    --------
    >>> # 30-day low-pass filter with 61-point window
    >>> weights = lanczos_weights(cutoff_freq=1/30, window=61, pass_type='low')
    >>> len(weights)
    61
    >>> weights.sum()  # Low-pass sums to 1
    1.0

    Notes
    -----
    The Lanczos filter is a windowed sinc filter:
    - Ideal low-pass filter: sinc(2*pi*fc*t)
    - Windowed by Lanczos sigma factor: sinc(pi*t/T)
    where fc is cutoff frequency and T is half-window length

    References
    ----------
    Duchon, C. E. (1979). Lanczos Filtering in One and Two Dimensions.
        Journal of Applied Meteorology, 18(8), 1016-1022.
    """
    # Map pass_type to integer for Fortran
    pass_type_int = 1 if pass_type == "low" else 2 if pass_type == "high" else 0

    if pass_type_int == 0:
        raise ValueError(f"pass_type must be 'low' or 'high', got {pass_type}")

    # Call Fortran routine
    weights, ierr = _lanczos_core.compute_lanczos_weights(
        cutoff_freq, window, pass_type_int
    )

    if ierr != 0:
        # Error already raised by C wrapper
        raise RuntimeError(
            f"Fortran compute_lanczos_weights failed with error code {ierr}"
        )

    return weights


def _apply_1d_filter(
    data: np.ndarray,
    weights: np.ndarray,
    axis: int = -1,
    boundary: Literal["reflect", "periodic", "constant"] = "reflect",
    fill_value: float = np.nan,
) -> np.ndarray:
    """
    Apply 1D convolution with filter weights along specified axis.

    Parameters
    ----------
    data : ndarray
        Input data array
    weights : ndarray
        1D filter weights
    axis : int, default=-1
        Axis along which to apply filter
    boundary : {'reflect', 'periodic', 'constant'}, default='reflect'
        Boundary handling mode
    fill_value : float, default=np.nan
        Fill value for 'constant' boundary mode

    Returns
    -------
    filtered : ndarray
        Filtered data with same shape as input
    """
    if boundary == "constant":
        raise NotImplementedError(
            "constant boundary mode not yet implemented in Fortran backend"
        )

    # Move axis to last position
    data = np.moveaxis(data, axis, -1)
    original_shape = data.shape

    # Flatten all dimensions except the last one
    data_2d = data.reshape(-1, data.shape[-1])
    n_series = data_2d.shape[0]
    n_points = data_2d.shape[1]

    # Prepare output
    filtered_2d = np.empty_like(data_2d)

    # Ensure data is Fortran-contiguous and float64
    weights = np.asfortranarray(weights, dtype=np.float64)

    # Apply filter to each series
    if boundary == "reflect":
        for i in range(n_series):
            data_series = np.asfortranarray(data_2d[i, :], dtype=np.float64)
            filtered_2d[i, :] = _lanczos_core.apply_filter_1d(data_series, weights)
    elif boundary == "periodic":
        for i in range(n_series):
            data_series = np.asfortranarray(data_2d[i, :], dtype=np.float64)
            filtered_2d[i, :] = _lanczos_core.apply_filter_1d_periodic(
                data_series, weights
            )
    else:
        raise ValueError(f"Unknown boundary mode: {boundary}")

    # Reshape back to original shape
    filtered = filtered_2d.reshape(original_shape)

    # Move axis back to original position
    filtered = np.moveaxis(filtered, -1, axis)

    return filtered


def _apply_2d_filter(
    data: np.ndarray,
    weights_x: np.ndarray,
    weights_y: np.ndarray,
    axes: tuple[int, int] = (-2, -1),
    boundary: Literal["reflect", "periodic", "constant"] = "reflect",
    fill_value: float = np.nan,
) -> np.ndarray:
    """
    Apply 2D separable convolution.

    The 2D filter is applied as two successive 1D filters (separable filtering).

    Parameters
    ----------
    data : ndarray
        Input data array
    weights_x : ndarray
        1D filter weights for first dimension
    weights_y : ndarray
        1D filter weights for second dimension
    axes : tuple of two ints, default=(-2, -1)
        Axes along which to apply filter
    boundary : str
        Boundary handling mode
    fill_value : float
        Fill value for constant boundary

    Returns
    -------
    filtered : ndarray
        Filtered data
    """
    # Apply filter along first axis
    filtered = _apply_1d_filter(
        data, weights_x, axis=axes[0], boundary=boundary, fill_value=fill_value
    )

    # Apply filter along second axis
    filtered = _apply_1d_filter(
        filtered, weights_y, axis=axes[1], boundary=boundary, fill_value=fill_value
    )

    return filtered


def lanczos_filter(
    data,
    cutoff_freq: float | dict[str, float],
    window: int | dict[str, int] | None = None,
    dim: str | int | Sequence[str | int] = -1,
    pass_type: Literal["low", "high"] = "low",
    boundary: Literal["reflect", "periodic", "constant"] = "reflect",
    fill_value: float = np.nan,
):
    """
    Apply Lanczos filter to data along specified dimension(s).

    Supports 1D, 2D, and multi-dimensional filtering with automatic handling
    of numpy arrays and xarray DataArrays.

    Parameters
    ----------
    data : array-like or xr.DataArray
        Input data to filter. Can be numpy array or xarray DataArray.
    cutoff_freq : float or dict
        Cutoff frequency in normalized units (0 to 0.5).
        - For 1D: single float (e.g., 1/30 for 30-timestep cutoff)
        - For 2D/nD: dict mapping dimension names to frequencies
          (e.g., {'time': 1/30, 'lat': 1/500})
    window : int, dict, or None
        Filter window length (number of weights, must be odd).
        - For 1D: single int
        - For 2D/nD: dict mapping dimension names to window sizes
        - If None: automatically set to 6/cutoff_freq + 1 (odd)
    dim : str, int, or sequence, default=-1
        Dimension(s) along which to filter.
        - For numpy: int or sequence of ints (axis indices)
        - For xarray: str or sequence of strs (dimension names)
    pass_type : {'low', 'high'}, default='low'
        Type of filter:
        - 'low': low-pass (smoothing)
        - 'high': high-pass (detrending)
    boundary : {'reflect', 'periodic', 'constant'}, default='reflect'
        Boundary handling:
        - 'reflect': mirror at edges (good for most cases)
        - 'periodic': wrap around (for global/cyclic data)
        - 'constant': pad with fill_value
    fill_value : float, default=np.nan
        Value for 'constant' boundary mode

    Returns
    -------
    filtered : same type as input
        Filtered data with same shape and type as input

    Examples
    --------
    >>> # 1D: Remove high-frequency noise from time series
    >>> import numpy as np
    >>> data = np.random.randn(1000)
    >>> smoothed = lanczos_filter(data, cutoff_freq=0.1, window=61)

    >>> # 2D: Spatial smoothing of a field
    >>> field = np.random.randn(180, 360)  # lat x lon
    >>> smooth_field = lanczos_filter(
    ...     field,
    ...     cutoff_freq={'lat': 1/500, 'lon': 1/500},
    ...     dim=['lat', 'lon']
    ... )

    >>> # xarray: Filter along time dimension
    >>> import xarray as xr
    >>> ds = xr.tutorial.open_dataset('air_temperature')
    >>> filtered_temp = lanczos_filter(
    ...     ds['air'],
    ...     cutoff_freq=1/30,
    ...     dim='time'
    ... )

    Notes
    -----
    For 2D/multi-dimensional filtering, the filter is applied separably
    (successively along each dimension), which is computationally efficient
    and produces isotropic results for equal cutoff frequencies.

    The default window size is chosen as 6/cutoff_freq, which provides
    a good balance between sharp cutoff and minimal edge effects.

    References
    ----------
    Duchon, C. E. (1979). Lanczos Filtering in One and Two Dimensions.
        Journal of Applied Meteorology, 18(8), 1016-1022.
    """
    # Check if xarray
    is_xarray = HAS_XARRAY and isinstance(data, xr.DataArray)

    if is_xarray:
        return _lanczos_filter_xarray(
            data, cutoff_freq, window, dim, pass_type, boundary, fill_value
        )
    else:
        return _lanczos_filter_numpy(
            data, cutoff_freq, window, dim, pass_type, boundary, fill_value
        )


def _lanczos_filter_numpy(
    data: np.ndarray,
    cutoff_freq: float | dict,
    window: int | dict | None,
    dim: int | Sequence[int],
    pass_type: str,
    boundary: str,
    fill_value: float,
) -> np.ndarray:
    """Lanczos filter for numpy arrays."""
    data = np.asarray(data, dtype=np.float64)

    # Handle single dimension case
    if isinstance(dim, int):
        dim = [dim]
    elif not isinstance(dim, (list, tuple)):
        dim = [dim]

    # Normalize negative axes
    ndim = data.ndim
    dim = [ax if ax >= 0 else ndim + ax for ax in dim]

    # Handle cutoff_freq and window
    if isinstance(cutoff_freq, dict):
        raise ValueError(
            "For numpy arrays, cutoff_freq must be a single float when filtering "
            "multiple dimensions. Use xarray DataArray for dimension-specific frequencies."
        )

    if isinstance(window, dict):
        raise ValueError(
            "For numpy arrays, window must be a single int when filtering "
            "multiple dimensions. Use xarray DataArray for dimension-specific windows."
        )

    # Auto-determine window if not provided
    if window is None:
        window = int(6.0 / cutoff_freq)
        if window % 2 == 0:
            window += 1

    # Compute weights
    weights = lanczos_weights(cutoff_freq, window, pass_type)

    # Apply filter along each dimension
    filtered = data.copy()
    for axis in dim:
        filtered = _apply_1d_filter(filtered, weights, axis, boundary, fill_value)

    return filtered


def _lanczos_filter_xarray(
    data: "xr.DataArray",
    cutoff_freq: float | dict,
    window: int | dict | None,
    dim: str | Sequence[str],
    pass_type: str,
    boundary: str,
    fill_value: float,
) -> "xr.DataArray":
    """Lanczos filter for xarray DataArrays."""
    # Handle single dimension case
    if isinstance(dim, str):
        dim = [dim]

    # Convert cutoff_freq and window to dicts if needed
    if not isinstance(cutoff_freq, dict):
        cutoff_freq = {d: cutoff_freq for d in dim}

    if window is None:
        window = {}
        for d in dim:
            w = int(6.0 / cutoff_freq[d])
            if w % 2 == 0:
                w += 1
            window[d] = w
    elif not isinstance(window, dict):
        window = {d: window for d in dim}

    # Apply filter along each dimension
    filtered = data.copy()
    for d in dim:
        if d not in cutoff_freq:
            raise ValueError(f"Dimension '{d}' not found in cutoff_freq")

        # Get axis index
        axis = data.dims.index(d)

        # Compute weights for this dimension
        weights = lanczos_weights(cutoff_freq[d], window[d], pass_type)

        # Apply filter
        filtered.values = _apply_1d_filter(
            filtered.values, weights, axis, boundary, fill_value
        )

    # Preserve attributes
    filtered.attrs = data.attrs.copy()
    filtered.attrs["lanczos_filtered"] = True
    filtered.attrs["lanczos_cutoff_freq"] = cutoff_freq
    filtered.attrs["lanczos_window"] = window
    filtered.attrs["lanczos_pass_type"] = pass_type

    return filtered


def lanczos_lowpass(
    data,
    cutoff_freq: float | dict,
    window: int | dict | None = None,
    dim: str | int | Sequence[str | int] = -1,
    **kwargs,
):
    """
    Apply Lanczos low-pass filter (smoothing, removes high frequencies).

    This is a convenience wrapper for lanczos_filter with pass_type='low'.

    Parameters
    ----------
    data : array-like or xr.DataArray
        Input data to filter
    cutoff_freq : float or dict
        Cutoff frequency (normalized, 0 to 0.5)
    window : int, dict, or None
        Filter window length (must be odd)
    dim : str, int, or sequence
        Dimension(s) along which to filter
    **kwargs
        Additional arguments passed to lanczos_filter

    Returns
    -------
    filtered : same type as input
        Low-pass filtered data

    Examples
    --------
    >>> # Smooth daily data with 30-day window
    >>> smoothed = lanczos_lowpass(daily_data, cutoff_freq=1/30)

    See Also
    --------
    lanczos_filter : General Lanczos filtering function
    lanczos_highpass : High-pass filtering
    """
    return lanczos_filter(data, cutoff_freq, window, dim, pass_type="low", **kwargs)


def lanczos_highpass(
    data,
    cutoff_freq: float | dict,
    window: int | dict | None = None,
    dim: str | int | Sequence[str | int] = -1,
    **kwargs,
):
    """
    Apply Lanczos high-pass filter (removes trends, keeps high frequencies).

    This is a convenience wrapper for lanczos_filter with pass_type='high'.

    Parameters
    ----------
    data : array-like or xr.DataArray
        Input data to filter
    cutoff_freq : float or dict
        Cutoff frequency (normalized, 0 to 0.5)
    window : int, dict, or None
        Filter window length (must be odd)
    dim : str, int, or sequence
        Dimension(s) along which to filter
    **kwargs
        Additional arguments passed to lanczos_filter

    Returns
    -------
    filtered : same type as input
        High-pass filtered data

    Examples
    --------
    >>> # Remove seasonal cycle (periods > 365 days)
    >>> anomalies = lanczos_highpass(daily_data, cutoff_freq=1/365)

    See Also
    --------
    lanczos_filter : General Lanczos filtering function
    lanczos_lowpass : Low-pass filtering
    """
    return lanczos_filter(data, cutoff_freq, window, dim, pass_type="high", **kwargs)


def lanczos_bandpass(
    data,
    freq_low: float,
    freq_high: float,
    window: int | None = None,
    dim: str | int | Sequence[str | int] = -1,
    **kwargs,
):
    """
    Apply Lanczos band-pass filter (keep only frequencies in range).

    Band-pass filtering is achieved by subtracting a low-pass filtered version
    from a high-pass filtered version.

    Parameters
    ----------
    data : array-like or xr.DataArray
        Input data to filter
    freq_low : float
        Lower cutoff frequency (normalized, 0 to 0.5)
    freq_high : float
        Upper cutoff frequency (normalized, 0 to 0.5)
        Must be > freq_low
    window : int or None
        Filter window length (must be odd)
        If None, determined from freq_low
    dim : str, int, or sequence
        Dimension(s) along which to filter
    **kwargs
        Additional arguments passed to lanczos_filter

    Returns
    -------
    filtered : same type as input
        Band-pass filtered data

    Examples
    --------
    >>> # Extract intraseasonal variability (30-90 days)
    >>> iso = lanczos_bandpass(daily_data, freq_low=1/90, freq_high=1/30)

    >>> # Extract ENSO signal (2-7 years from daily data)
    >>> enso = lanczos_bandpass(sst_daily, freq_low=1/(7*365), freq_high=1/(2*365))

    See Also
    --------
    lanczos_filter : General Lanczos filtering function
    lanczos_lowpass : Low-pass filtering
    lanczos_highpass : High-pass filtering
    """
    if freq_low >= freq_high:
        raise ValueError(f"freq_low ({freq_low}) must be < freq_high ({freq_high})")

    # Auto-determine window from lower frequency (needs longer window)
    if window is None:
        window = int(6.0 / freq_low)
        if window % 2 == 0:
            window += 1

    # Apply two filters: high-pass at freq_low, low-pass at freq_high
    # Band-pass = high-pass(freq_low) - high-pass(freq_high)
    # Or equivalently: low-pass(freq_high) - low-pass(freq_low)

    highpassed = lanczos_highpass(data, freq_low, window, dim, **kwargs)
    lowpassed_hp = lanczos_lowpass(highpassed, freq_high, window, dim, **kwargs)

    return lowpassed_hp
