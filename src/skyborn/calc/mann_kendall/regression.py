"""Regression-side batch helpers used by the Mann-Kendall API layer."""

import numpy as np


def _linregress_slope_batch(data_2d: np.ndarray, x: np.ndarray) -> np.ndarray:
    """
    Vectorized linear regression slope calculation.

    Much faster than calling stats.linregress for each series individually.
    """
    x_mean = np.mean(x)
    x_centered = x - x_mean
    numerator = x_centered @ data_2d
    denominator = np.sum(x_centered**2)
    return numerator / denominator


def _theil_std_error_batch(
    data_2d: np.ndarray, x: np.ndarray, slopes: np.ndarray
) -> np.ndarray:
    """Vectorized standard error calculation for Theil-Sen slopes."""
    time_steps = data_2d.shape[0]
    workspace = np.empty_like(data_2d, dtype=np.float64)
    np.multiply(x[:, np.newaxis], slopes[np.newaxis, :], out=workspace)
    workspace *= -1.0
    workspace += data_2d
    intercepts = np.median(workspace, axis=0)
    workspace -= intercepts
    return np.std(workspace, axis=0) / np.sqrt(time_steps)


def _linregress_std_error_batch(
    data_2d: np.ndarray, x: np.ndarray, slopes: np.ndarray
) -> np.ndarray:
    """
    Vectorized residual standard error calculation for linear regression.

    This mirrors the scalar `mann_kendall_test` implementation, which returns
    the standard error of the detrended residuals rather than the regression
    slope standard error.
    """
    time_steps = data_2d.shape[0]
    x_mean = np.mean(x)
    y_mean = np.mean(data_2d, axis=0)
    intercepts = y_mean - slopes * x_mean

    residuals = np.empty_like(data_2d, dtype=np.float64)
    np.multiply(x[:, np.newaxis], slopes[np.newaxis, :], out=residuals)
    residuals += intercepts
    residuals *= -1.0
    residuals += data_2d
    return np.std(residuals, axis=0) / np.sqrt(time_steps)
