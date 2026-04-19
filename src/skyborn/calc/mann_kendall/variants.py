"""Variant-specific batch transforms for Mann-Kendall test families."""

from typing import Optional, Tuple

import numpy as np
import scipy.stats as stats

from .bindings import _score_variance_batch


def _variance_without_ties(sample_size: int) -> float:
    """Return the no-tie MK variance formula for the given sample size."""
    return sample_size * (sample_size - 1) * (2 * sample_size + 5) / 18.0


def _lag1_autocorrelation_batch(data_2d: np.ndarray) -> np.ndarray:
    """Compute lag-1 autocorrelation for each dense series column."""
    centered = np.asarray(data_2d, dtype=np.float64) - np.mean(data_2d, axis=0)
    acov0 = np.sum(centered * centered, axis=0)
    acov1 = np.sum(centered[:-1] * centered[1:], axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return acov1 / acov0


def _pre_whiten_batch(data_2d: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Pre-whiten dense 2D series and flag degenerate zero-variance columns."""
    acf1 = _lag1_autocorrelation_batch(data_2d)
    whitened = np.asarray(data_2d[1:], dtype=np.float64).copy()
    whitened -= np.asarray(data_2d[:-1], dtype=np.float64) * acf1[np.newaxis, :]
    degenerate = ~np.isfinite(acf1)
    return whitened, degenerate


def _pre_whiten_series(data_1d: np.ndarray) -> Tuple[np.ndarray, bool]:
    """Pre-whiten one finite series and flag degenerate zero-variance input."""
    whitened, degenerate = _pre_whiten_batch(
        np.asarray(data_1d, dtype=np.float64).reshape(-1, 1)
    )
    return whitened[:, 0], bool(degenerate[0])


def _pre_whiten_score_variance_batch(
    data_2d: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return S and variance of the whitened series without exporting the intermediate."""
    whitened, degenerate = _pre_whiten_batch(data_2d)
    stat_n = whitened.shape[0]
    s_values = np.zeros(data_2d.shape[1], dtype=np.float64)
    var_values = np.full(
        data_2d.shape[1], _variance_without_ties(stat_n), dtype=np.float64
    )
    nondegenerate_mask = ~degenerate
    if np.any(nondegenerate_mask):
        s_values[nondegenerate_mask], var_values[nondegenerate_mask] = (
            _score_variance_batch(whitened[:, nondegenerate_mask], modified=False)
        )
    return s_values, var_values


def _trend_free_pre_whiten_batch(
    data_2d: np.ndarray, slopes: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Apply trend-free pre-whitening and return the retrended dense series batch."""
    n = data_2d.shape[0]
    detrended = np.asarray(data_2d, dtype=np.float64).copy()
    detrended -= (
        np.arange(1, n + 1, dtype=np.float64)[:, np.newaxis] * slopes[np.newaxis, :]
    )
    acf1 = _lag1_autocorrelation_batch(detrended)
    transformed = detrended[1:].copy()
    transformed -= detrended[:-1] * acf1[np.newaxis, :]
    transformed += (
        np.arange(1, n, dtype=np.float64)[:, np.newaxis] * slopes[np.newaxis, :]
    )
    degenerate = ~np.isfinite(acf1)
    return transformed, degenerate


def _trend_free_pre_whiten_score_variance_batch(
    data_2d: np.ndarray, slopes: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Return S and variance after trend-free pre-whitening without exporting intermediates."""
    transformed, degenerate = _trend_free_pre_whiten_batch(data_2d, slopes)
    stat_n = transformed.shape[0]
    s_values = np.zeros(data_2d.shape[1], dtype=np.float64)
    var_values = np.full(
        data_2d.shape[1], _variance_without_ties(stat_n), dtype=np.float64
    )
    nondegenerate_mask = ~degenerate
    if np.any(nondegenerate_mask):
        s_values[nondegenerate_mask], var_values[nondegenerate_mask] = (
            _score_variance_batch(transformed[:, nondegenerate_mask], modified=False)
        )
    return s_values, var_values


def _hamed_rao_variance_batch(
    data_2d: np.ndarray,
    base_var_values: np.ndarray,
    slopes: np.ndarray,
    alpha: float,
    lag: Optional[int] = None,
) -> np.ndarray:
    """Apply the Hamed-Rao variance correction to a dense 2D series batch."""
    n = data_2d.shape[0]
    max_lag = n - 1 if lag is None else max(0, min(int(lag), n - 1))
    if n < 3:
        return np.asarray(base_var_values, dtype=np.float64)

    detrended = np.asarray(data_2d, dtype=np.float64).copy()
    detrended -= (
        np.arange(1, n + 1, dtype=np.float64)[:, np.newaxis] * slopes[np.newaxis, :]
    )
    ranked = np.asarray(stats.rankdata(detrended, axis=0), dtype=np.float64)
    centered = ranked - np.mean(ranked, axis=0)
    denom = np.sum(centered * centered, axis=0)

    interval = stats.norm.ppf(1 - alpha / 2) / np.sqrt(n)
    sni = np.zeros(data_2d.shape[1], dtype=np.float64)

    with np.errstate(divide="ignore", invalid="ignore"):
        for lag_index in range(1, max_lag + 1):
            acf = (
                np.sum(
                    centered[:-lag_index] * centered[lag_index:],
                    axis=0,
                    dtype=np.float64,
                )
                / denom
            )
            inside_interval = (acf <= interval) & (acf >= -interval)
            weight = (n - lag_index) * (n - lag_index - 1) * (n - lag_index - 2)
            sni += np.where(inside_interval, 0.0, weight * acf)

    n_ns = 1.0 + (2.0 * sni) / (n * (n - 1) * (n - 2))
    return np.asarray(base_var_values, dtype=np.float64) * n_ns


def _reshape_seasonal_1d(values_1d: np.ndarray, period: int) -> np.ndarray:
    """Pad a 1D seasonal series with NaN and reshape it to (cycles, period)."""
    values = np.asarray(values_1d, dtype=np.float64)
    if values.ndim != 1:
        raise ValueError("Seasonal helpers expect a 1D array.")
    if period <= 0:
        raise ValueError("period must be a positive integer.")
    remainder = values.size % period
    if remainder != 0:
        values = np.pad(
            values,
            (0, period - remainder),
            mode="constant",
            constant_values=np.nan,
        )
    return values.reshape(values.size // period, period)


def _seasonal_score_variance_scalar(
    values_1d: np.ndarray, period: int
) -> Tuple[float, float, float]:
    """Return seasonal S, variance, and tau denominator for one 1D series."""
    matrix = _reshape_seasonal_1d(values_1d, period)
    s_value = 0.0
    var_s = 0.0
    denom = 0.0

    for season in range(matrix.shape[1]):
        season_values = matrix[:, season]
        finite = season_values[np.isfinite(season_values)]
        n = finite.size
        if n == 0:
            continue

        season_s, season_var = _score_variance_batch(
            finite.reshape(-1, 1), modified=False
        )
        s_value += float(season_s[0])
        var_s += float(season_var[0])
        denom += 0.5 * n * (n - 1)

    return s_value, var_s, denom


def _seasonal_score_variance_batch(
    data_2d: np.ndarray, period: int
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Return seasonal S and variance for a clean 2D batch with common period."""
    if period <= 0:
        raise ValueError("period must be a positive integer.")

    data = np.asarray(data_2d, dtype=np.float64)
    if data.ndim != 2:
        raise ValueError("Expected a 2D array for seasonal batch helpers.")

    s_values = np.zeros(data.shape[1], dtype=np.float64)
    var_values = np.zeros(data.shape[1], dtype=np.float64)
    denom = 0.0

    for season in range(period):
        season_data = data[season::period, :]
        n = season_data.shape[0]
        if n == 0:
            continue

        season_s, season_var = _score_variance_batch(season_data, modified=False)
        s_values += season_s
        var_values += season_var
        denom += 0.5 * n * (n - 1)

    return s_values, var_values, denom


def _seasonal_sens_slope_scalar(
    values_1d: np.ndarray, period: int
) -> Tuple[float, float]:
    """Return the seasonal Sen slope and intercept for one 1D series."""
    matrix = _reshape_seasonal_1d(values_1d, period)
    slopes = []

    for season in range(matrix.shape[1]):
        season_values = matrix[:, season]
        n = season_values.size
        for i in range(n - 1):
            if not np.isfinite(season_values[i]):
                continue
            for j in range(i + 1, n):
                if not np.isfinite(season_values[j]):
                    continue
                slopes.append((season_values[j] - season_values[i]) / (j - i))

    if slopes:
        slope = float(np.nanmedian(np.asarray(slopes, dtype=np.float64)))
    else:
        slope = np.nan

    values = np.asarray(values_1d, dtype=np.float64)
    finite_index = np.arange(values.size, dtype=np.float64)[np.isfinite(values)]
    intercept = float(np.nanmedian(values) - np.median(finite_index) / period * slope)
    return slope, intercept


def _grouped_time_index(n_steps: int, n_groups: int) -> np.ndarray:
    """Return the flattened grouped time coordinate used by grouped MK variants."""
    return np.arange(n_steps * n_groups, dtype=np.float64) / float(n_groups)


def _multivariate_score_variance_scalar(
    matrix_2d: np.ndarray,
) -> Tuple[float, float, float]:
    """Return multivariate S, variance, and tau denominator for one matrix."""
    matrix = np.asarray(matrix_2d, dtype=np.float64)
    if matrix.ndim == 1:
        matrix = matrix.reshape(-1, 1)
    if matrix.ndim != 2:
        raise ValueError("Expected a 1D or 2D array for multivariate statistics.")

    s_value = 0.0
    var_s = 0.0
    denom = 0.0

    for group in range(matrix.shape[1]):
        group_values = matrix[:, group]
        finite = group_values[np.isfinite(group_values)]
        n = finite.size
        if n == 0:
            continue

        group_s, group_var = _score_variance_batch(
            finite.reshape(-1, 1), modified=False
        )
        s_value += float(group_s[0])
        var_s += float(group_var[0])
        denom += 0.5 * n * (n - 1)

    return s_value, var_s, denom


def _multivariate_score_variance_batch(
    data_3d: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Return multivariate S and variance for a clean grouped batch."""
    data = np.asarray(data_3d, dtype=np.float64)
    if data.ndim != 3:
        raise ValueError("Expected a 3D array with shape (time, group, series).")

    s_values = np.zeros(data.shape[2], dtype=np.float64)
    var_values = np.zeros(data.shape[2], dtype=np.float64)
    denom = 0.0

    for group in range(data.shape[1]):
        group_s, group_var = _score_variance_batch(data[:, group, :], modified=False)
        s_values += group_s
        var_values += group_var
        n = data[:, group, :].shape[0]
        denom += 0.5 * n * (n - 1)

    return s_values, var_values, denom


def _multivariate_sens_slope_scalar(matrix_2d: np.ndarray) -> Tuple[float, float]:
    """Return multivariate Sen slope and intercept for one grouped series matrix."""
    matrix = np.asarray(matrix_2d, dtype=np.float64)
    if matrix.ndim == 1:
        matrix = matrix.reshape(-1, 1)
    if matrix.ndim != 2:
        raise ValueError("Expected a 1D or 2D array for multivariate slope.")

    slopes = []
    for group in range(matrix.shape[1]):
        group_values = matrix[:, group]
        n = group_values.size
        for i in range(n - 1):
            if not np.isfinite(group_values[i]):
                continue
            for j in range(i + 1, n):
                if not np.isfinite(group_values[j]):
                    continue
                slopes.append((group_values[j] - group_values[i]) / (j - i))

    if slopes:
        slope = float(np.nanmedian(np.asarray(slopes, dtype=np.float64)))
    else:
        slope = np.nan

    flattened = matrix.reshape(-1)
    finite_index = np.arange(flattened.size, dtype=np.float64)[np.isfinite(flattened)]
    intercept = float(
        np.nanmedian(flattened) - np.median(finite_index) / matrix.shape[1] * slope
    )
    return slope, intercept


def _correlated_multivariate_stats_scalar(
    matrix_2d: np.ndarray,
) -> Tuple[float, float, float]:
    """Return correlated multivariate S, variance, and tau denominator for one matrix."""
    x = np.asarray(matrix_2d, dtype=np.float64)
    if x.ndim != 2:
        raise ValueError("Expected a 2D matrix for correlated multivariate statistics.")

    valid_rows = ~np.isnan(x).any(axis=1)
    x = x[valid_rows]
    n, c = x.shape if x.ndim == 2 else (len(x), 1)

    s_value = 0.0
    denom = 0.0
    for i in range(c):
        season_s, _ = _score_variance_batch(x[:, i].reshape(-1, 1), modified=False)
        s_value += float(season_s[0])
        denom += 0.5 * n * (n - 1)

    gamma = np.ones((c, c), dtype=np.float64)

    def _k_stat(left: np.ndarray, right: np.ndarray) -> float:
        k_value = 0.0
        for i in range(n - 1):
            j = np.arange(i, n)
            k_value += np.sum(np.sign((left[j] - left[i]) * (right[j] - right[i])))
        return float(k_value)

    def _r_stat(values: np.ndarray) -> np.ndarray:
        result = np.empty(n, dtype=np.float64)
        indices = np.arange(n)
        for j in range(n):
            s = np.sum(np.sign(values[j] - values[indices]))
            result[j] = (n + 1 + s) / 2.0
        return result

    for i in range(1, c):
        for j in range(i):
            k_value = _k_stat(x[:, i], x[:, j])
            ri = _r_stat(x[:, i])
            rj = _r_stat(x[:, j])
            gamma[i, j] = (k_value + 4.0 * np.sum(ri * rj) - n * (n + 1) ** 2) / 3.0
            gamma[j, i] = gamma[i, j]

    for i in range(c):
        k_value = _k_stat(x[:, i], x[:, i])
        ri = _r_stat(x[:, i])
        gamma[i, i] = (k_value + 4.0 * np.sum(ri * ri) - n * (n + 1) ** 2) / 3.0

    return s_value, float(np.sum(gamma)), denom


def _correlated_seasonal_stats_scalar(
    values_1d: np.ndarray, period: int
) -> Tuple[float, float, float]:
    """Return correlated seasonal S, variance, and tau denominator for one 1D series."""
    matrix = _reshape_seasonal_1d(values_1d, period)
    return _correlated_multivariate_stats_scalar(matrix)
