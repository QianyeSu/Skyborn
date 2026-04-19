"""Independent partial Mann-Kendall interfaces.

This module intentionally stays separate from ``core.py`` because partial
Mann-Kendall uses two input series: a response variable and a covariate.
That changes the missing-value semantics and the broadcast rules compared with
the single-series MK families implemented in the core module.
"""

import warnings
from typing import Any, Dict, Optional, Union

import numpy as np
import scipy.stats as stats
import xarray as xr

from .bindings import _score_variance_batch, _sen_slope_batch
from .regression import (
    _linregress_slope_batch,
    _linregress_std_error_batch,
    _theil_std_error_batch,
)

__all__ = [
    "partial_mann_kendall_test",
    "partial_mann_kendall_multidim",
    "partial_mann_kendall_xarray",
    "partial_test",
    "partial_multidim",
]


def _variance_without_ties(sample_size: int) -> float:
    """Return the no-tie Mann-Kendall variance for one sample size."""
    return sample_size * (sample_size - 1) * (2 * sample_size + 5) / 18.0


def _validate_partial_method(method: str) -> str:
    """Validate and normalize the response-slope method."""
    normalized = method.lower()
    if normalized not in {"theilslopes", "linregress"}:
        raise ValueError(f"Unknown method: {method}. Use 'theilslopes' or 'linregress'")
    return normalized


def _empty_partial_result(n_total: int, n_valid: int) -> Dict[str, Union[float, bool]]:
    """Return the partial MK NaN result payload for insufficient samples."""
    return {
        "trend": np.nan,
        "h": False,
        "p": np.nan,
        "z": np.nan,
        "tau": np.nan,
        "s": np.nan,
        "var_s": np.nan,
        "intercept": np.nan,
        "std_error": np.nan,
        "n": float(n_valid),
        "n_dropped": float(n_total - n_valid),
    }


def _partial_k_stat_batch(
    response_2d: np.ndarray, covariate_2d: np.ndarray
) -> np.ndarray:
    """Return the partial-MK K statistic for each dense column pair."""
    n_time, n_series = response_2d.shape
    k_values = np.zeros(n_series, dtype=np.float64)

    for start in range(n_time - 1):
        dx = response_2d[start:, :] - response_2d[start, :]
        dy = covariate_2d[start:, :] - covariate_2d[start, :]
        k_values += np.sum(np.sign(dx * dy), axis=0)

    return k_values


def _partial_rank_cross_sum_batch(
    response_2d: np.ndarray, covariate_2d: np.ndarray
) -> np.ndarray:
    """Return ``sum(R_x * R_y)`` for each dense column pair."""
    n_time = response_2d.shape[0]
    rank_x = np.empty_like(response_2d, dtype=np.float64)
    rank_y = np.empty_like(covariate_2d, dtype=np.float64)

    for row in range(n_time):
        rank_x[row, :] = (
            n_time
            + 1.0
            + np.sum(np.sign(response_2d[row : row + 1, :] - response_2d), axis=0)
        ) / 2.0
        rank_y[row, :] = (
            n_time
            + 1.0
            + np.sum(np.sign(covariate_2d[row : row + 1, :] - covariate_2d), axis=0)
        ) / 2.0

    return np.sum(rank_x * rank_y, axis=0)


def _partial_statistics_batch(
    response_2d: np.ndarray, covariate_2d: np.ndarray
) -> Dict[str, np.ndarray]:
    """Return partial-MK statistics for a clean dense 2D batch."""
    response = np.asarray(response_2d, dtype=np.float64)
    covariate = np.asarray(covariate_2d, dtype=np.float64)

    if response.ndim != 2 or covariate.ndim != 2:
        raise ValueError("Partial MK batch helpers expect two 2D arrays.")
    if response.shape != covariate.shape:
        raise ValueError("Response and covariate batches must have the same shape.")

    n_time = response.shape[0]
    if n_time < 3:
        raise ValueError("Partial MK batch helpers need at least 3 time steps.")

    response_scores, _ = _score_variance_batch(response, modified=False)
    covariate_scores, _ = _score_variance_batch(covariate, modified=False)

    k_values = _partial_k_stat_batch(response, covariate)
    rank_cross_sum = _partial_rank_cross_sum_batch(response, covariate)
    sigma = (k_values + 4.0 * rank_cross_sum - n_time * (n_time + 1.0) ** 2) / 3.0

    variance_no_ties = _variance_without_ties(n_time)
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = sigma / variance_no_ties

    s_values = response_scores - rho * covariate_scores
    var_values = (1.0 - rho**2) * variance_no_ties
    tau_values = response_scores / (0.5 * n_time * (n_time - 1))

    return {
        "s": np.asarray(s_values, dtype=np.float64),
        "var_s": np.asarray(var_values, dtype=np.float64),
        "tau": np.asarray(tau_values, dtype=np.float64),
    }


def _theil_slope_intercept_scalar(
    response_1d: np.ndarray,
) -> tuple[float, float, np.ndarray, np.ndarray]:
    """Return Theil-Sen slope/intercept using the response's original time index."""
    response = np.asarray(response_1d, dtype=np.float64)
    valid = np.isfinite(response)
    x = np.arange(response.size, dtype=np.float64)[valid]
    y = response[valid]
    n = y.size

    if n < 2:
        return np.nan, np.nan, x, y

    slopes = np.empty(n * (n - 1) // 2, dtype=np.float64)
    offset = 0
    for start in range(n - 1):
        count = n - start - 1
        slopes[offset : offset + count] = (y[start + 1 :] - y[start]) / (
            x[start + 1 :] - x[start]
        )
        offset += count

    slope = float(np.nanmedian(slopes))
    intercept = float(np.nanmedian(y) - np.median(x) * slope)
    return slope, intercept, x, y


def _response_trend_scalar(
    response_1d: np.ndarray, method: str
) -> tuple[float, float, float]:
    """Return response-only slope, intercept, and std_error for one 1D series."""
    normalized_method = _validate_partial_method(method)
    response = np.asarray(response_1d, dtype=np.float64)
    valid = np.isfinite(response)
    x = np.arange(response.size, dtype=np.float64)[valid]
    y = response[valid]

    if y.size < 2:
        return np.nan, np.nan, np.nan

    if normalized_method == "theilslopes":
        slope, intercept, _, _ = _theil_slope_intercept_scalar(response)
        std_error = float(
            _theil_std_error_batch(
                y.reshape(-1, 1), x, np.asarray([slope], dtype=np.float64)
            )[0]
        )
        return slope, intercept, std_error

    slope, intercept = stats.linregress(x, y)[:2]
    std_error = float(
        _linregress_std_error_batch(
            y.reshape(-1, 1), x, np.asarray([slope], dtype=np.float64)
        )[0]
    )
    return float(slope), float(intercept), std_error


def _response_trend_batch(
    response_2d: np.ndarray, method: str
) -> Dict[str, np.ndarray]:
    """Return response-only slope, intercept, and std_error for a clean 2D batch."""
    normalized_method = _validate_partial_method(method)
    response = np.asarray(response_2d, dtype=np.float64)
    if response.ndim != 2:
        raise ValueError("Expected a 2D array for response trend batch estimation.")

    x = np.arange(response.shape[0], dtype=np.float64)

    if normalized_method == "theilslopes":
        slopes = _sen_slope_batch(response)
        intercepts = np.median(response, axis=0) - np.median(x) * slopes
        std_errors = _theil_std_error_batch(response, x, slopes)
    else:
        slopes = _linregress_slope_batch(response, x)
        intercepts = np.mean(response, axis=0) - slopes * np.mean(x)
        std_errors = _linregress_std_error_batch(response, x, slopes)

    return {
        "trend": np.asarray(slopes, dtype=np.float64),
        "intercept": np.asarray(intercepts, dtype=np.float64),
        "std_error": np.asarray(std_errors, dtype=np.float64),
    }


def partial_mann_kendall_test(
    data: np.ndarray,
    covariate: np.ndarray,
    alpha: float = 0.05,
    method: str = "theilslopes",
) -> Dict[str, Union[float, bool]]:
    """Perform the partial Mann-Kendall test for one response/covariate pair."""
    normalized_method = _validate_partial_method(method)

    response = np.asarray(data, dtype=np.float64)
    cov = np.asarray(covariate, dtype=np.float64)
    if response.ndim != 1 or cov.ndim != 1:
        raise ValueError("partial_mann_kendall_test expects two 1D arrays.")
    if response.shape != cov.shape:
        raise ValueError("Response and covariate must have the same shape.")

    joint_valid = np.isfinite(response) & np.isfinite(cov)
    n_valid = int(np.sum(joint_valid))
    if n_valid < 3:
        warnings.warn("Need at least 3 joint-valid data points for partial MK test")
        return _empty_partial_result(response.size, n_valid)

    stats_values = _partial_statistics_batch(
        response[joint_valid].reshape(-1, 1),
        cov[joint_valid].reshape(-1, 1),
    )
    s_value = float(stats_values["s"][0])
    var_s = float(stats_values["var_s"][0])
    tau = float(stats_values["tau"][0])

    with np.errstate(divide="ignore", invalid="ignore"):
        z = s_value / np.sqrt(var_s)
    p_value = 2.0 * (1.0 - stats.norm.cdf(abs(z)))
    h = bool(abs(z) > stats.norm.ppf(1.0 - alpha / 2.0)) if np.isfinite(z) else False

    slope, intercept, std_error = _response_trend_scalar(response, normalized_method)

    return {
        "trend": slope,
        "h": h,
        "p": float(p_value),
        "z": float(z),
        "tau": tau,
        "s": s_value,
        "var_s": var_s,
        "intercept": intercept,
        "std_error": std_error,
        "n": float(n_valid),
        "n_dropped": float(response.size - n_valid),
    }


def _broadcast_numpy_covariate(
    data: np.ndarray, covariate: np.ndarray, axis: int
) -> tuple[np.ndarray, np.ndarray]:
    """Broadcast NumPy covariate input onto the response array."""
    response = np.asarray(data, dtype=np.float64)
    cov = np.asarray(covariate, dtype=np.float64)
    time_axis = int(axis)

    response_time_first = np.moveaxis(response, time_axis, 0)

    if cov.ndim == 1:
        if cov.shape[0] != response.shape[time_axis]:
            raise ValueError(
                "1D covariate must have the same length as the response time axis."
            )
        cov_time_first = cov.reshape((cov.shape[0],) + (1,) * (response.ndim - 1))
    elif cov.shape == response.shape:
        cov_time_first = np.moveaxis(cov, time_axis, 0)
    else:
        raise ValueError(
            "For NumPy inputs, covariate must be 1D with the same time length or "
            "have the same shape as data."
        )

    return response_time_first, np.broadcast_to(
        cov_time_first, response_time_first.shape
    )


def partial_mann_kendall_multidim(
    data: np.ndarray,
    covariate: np.ndarray,
    axis: Union[int, str] = 0,
    alpha: float = 0.05,
    method: str = "theilslopes",
    chunk_size: Optional[int] = None,
    dim: Optional[Union[int, str]] = None,
) -> Dict[str, np.ndarray]:
    """Apply the partial Mann-Kendall test along one axis of a NumPy array."""
    if dim is not None:
        axis = dim
    if isinstance(axis, str):
        raise ValueError("String axis names are only supported for xarray inputs.")

    normalized_method = _validate_partial_method(method)
    response_time_first, cov_time_first = _broadcast_numpy_covariate(
        np.asarray(data, dtype=np.float64),
        np.asarray(covariate, dtype=np.float64),
        int(axis),
    )

    if response_time_first.ndim == 1:
        return partial_mann_kendall_test(
            response_time_first,
            cov_time_first,
            alpha=alpha,
            method=normalized_method,
        )

    time_steps = response_time_first.shape[0]
    spatial_shape = response_time_first.shape[1:]
    n_points = int(np.prod(spatial_shape))

    response_2d = response_time_first.reshape(time_steps, n_points)
    cov_2d = cov_time_first.reshape(time_steps, n_points)

    results_flat = {
        "trend": np.full(n_points, np.nan, dtype=np.float64),
        "h": np.zeros(n_points, dtype=bool),
        "p": np.full(n_points, np.nan, dtype=np.float64),
        "z": np.full(n_points, np.nan, dtype=np.float64),
        "tau": np.full(n_points, np.nan, dtype=np.float64),
        "s": np.full(n_points, np.nan, dtype=np.float64),
        "var_s": np.full(n_points, np.nan, dtype=np.float64),
        "intercept": np.full(n_points, np.nan, dtype=np.float64),
        "std_error": np.full(n_points, np.nan, dtype=np.float64),
        "n": np.zeros(n_points, dtype=np.float64),
        "n_dropped": np.full(n_points, float(time_steps), dtype=np.float64),
    }

    if chunk_size is None:
        max_memory_mb = 200
        bytes_per_series = max(1, time_steps) * 8 * 6
        chunk_size = max(
            100, min(n_points, (max_memory_mb * 1024 * 1024) // bytes_per_series)
        )

    for start_idx in range(0, n_points, chunk_size):
        end_idx = min(start_idx + chunk_size, n_points)
        response_chunk = response_2d[:, start_idx:end_idx]
        cov_chunk = cov_2d[:, start_idx:end_idx]

        joint_valid_counts = np.sum(
            np.isfinite(response_chunk) & np.isfinite(cov_chunk), axis=0
        )
        response_valid_counts = np.sum(np.isfinite(response_chunk), axis=0)
        valid_series = joint_valid_counts >= 3
        if not np.any(valid_series):
            continue

        local_indices = np.where(valid_series)[0]
        clean_local = local_indices[joint_valid_counts[local_indices] == time_steps]

        if clean_local.size > 0:
            response_clean = np.asarray(
                response_chunk[:, clean_local], dtype=np.float64
            )
            cov_clean = np.asarray(cov_chunk[:, clean_local], dtype=np.float64)

            partial_stats = _partial_statistics_batch(response_clean, cov_clean)
            response_trend = _response_trend_batch(response_clean, normalized_method)

            with np.errstate(divide="ignore", invalid="ignore"):
                z_values = partial_stats["s"] / np.sqrt(partial_stats["var_s"])
            p_values = 2.0 * (1.0 - stats.norm.cdf(np.abs(z_values)))
            h_values = np.abs(z_values) > stats.norm.ppf(1.0 - alpha / 2.0)

            global_clean = start_idx + clean_local
            results_flat["trend"][global_clean] = response_trend["trend"]
            results_flat["h"][global_clean] = h_values
            results_flat["p"][global_clean] = p_values
            results_flat["z"][global_clean] = z_values
            results_flat["tau"][global_clean] = partial_stats["tau"]
            results_flat["s"][global_clean] = partial_stats["s"]
            results_flat["var_s"][global_clean] = partial_stats["var_s"]
            results_flat["intercept"][global_clean] = response_trend["intercept"]
            results_flat["std_error"][global_clean] = response_trend["std_error"]
            results_flat["n"][global_clean] = float(time_steps)
            results_flat["n_dropped"][global_clean] = 0.0

        fallback_local = local_indices[joint_valid_counts[local_indices] != time_steps]
        for offset in fallback_local:
            result = partial_mann_kendall_test(
                response_chunk[:, offset],
                cov_chunk[:, offset],
                alpha=alpha,
                method=normalized_method,
            )
            global_idx = start_idx + int(offset)
            for key, value in result.items():
                results_flat[key][global_idx] = value

    return {
        key: (
            values.reshape(spatial_shape)
            if key != "h"
            else values.reshape(spatial_shape)
        )
        for key, values in results_flat.items()
    }


def partial_mann_kendall_xarray(
    data: xr.DataArray,
    covariate: xr.DataArray,
    dim: str = "time",
    alpha: float = 0.05,
    method: str = "theilslopes",
) -> xr.Dataset:
    """Apply the partial Mann-Kendall test to xarray data."""
    data_aligned, covariate_aligned = xr.align(
        data, covariate, join="exact", copy=False
    )
    covariate_broadcast = covariate_aligned.broadcast_like(data_aligned)

    time_axis = data_aligned.get_axis_num(dim)
    results = partial_mann_kendall_multidim(
        data_aligned.values,
        covariate_broadcast.values,
        axis=time_axis,
        alpha=alpha,
        method=method,
    )

    output_dims = [name for name in data_aligned.dims if name != dim]
    output_coords = {name: data_aligned.coords[name] for name in output_dims}

    if not output_dims:
        data_vars = {}
        for key, value in results.items():
            scalar_value = value if np.isscalar(value) else np.asarray(value).item()
            data_vars[key] = ([], scalar_value)
    else:
        data_vars = {
            key: (output_dims, np.asarray(value)) for key, value in results.items()
        }

    return xr.Dataset(
        data_vars=data_vars,
        coords=output_coords,
        attrs={
            "title": "Partial Mann-Kendall Trend Analysis",
            "alpha": alpha,
            "method": method,
            "analyzed_dim": dim,
            "covariate_dims": str(covariate_aligned.dims),
        },
    )


partial_test = partial_mann_kendall_test
partial_multidim = partial_mann_kendall_multidim
