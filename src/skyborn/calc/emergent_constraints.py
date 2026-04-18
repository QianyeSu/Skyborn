"""
Emergent Constraint Methods for Climate Data Analysis.

This module implements the emergent constraint method based on Cox et al. (2013).
It provides functions for calculating probability density functions and performing
constrained projections on climate model data.

References
----------
Cox, P. M., et al. (2013). Sensitivity of tropical carbon to climate change
constrained by carbon dioxide variability. Nature, 494(7437), 341-344.

Implementation adapted from:
https://github.com/blackcata/Emergent_Constraints/tree/master

Author: KM.Noh
Date: 2023.03.15
Modified for Skyborn package with type annotations and improved naming
"""

import warnings
from typing import Tuple, Union

import numpy as np
import xarray as xr

__all__ = [
    "gaussian_pdf",
    "emergent_constraint_posterior",
    "emergent_constraint_prior",
    "calc_GAUSSIAN_PDF",
    "calc_PDF_EC",
    "find_std_from_PDF",
    "calc_PDF_EC_PRIOR",
]


def _preferred_float_dtype(*values) -> np.dtype:
    """Pick a floating work dtype from the provided inputs."""
    preferred = None

    for value in values:
        raw = getattr(value, "values", value)
        array = np.asarray(raw)
        dtype = np.dtype(array.dtype)

        if not np.issubdtype(dtype, np.floating):
            continue

        if dtype.itemsize < np.dtype(np.float32).itemsize:
            dtype = np.dtype(np.float32)

        if preferred is None or dtype.itemsize > preferred.itemsize:
            preferred = dtype

    return preferred if preferred is not None else np.dtype(np.float32)


def gaussian_pdf(
    mu: float, sigma: float, x: Union[np.ndarray, float]
) -> Union[np.ndarray, float]:
    """
    Calculate Gaussian probability density function.

    Parameters
    ----------
    mu : float
        Mean of the Gaussian distribution.
    sigma : float
        Standard deviation of the Gaussian distribution.
    x : Union[np.ndarray, float]
        Input values at which to evaluate the PDF.

    Returns
    -------
    Union[np.ndarray, float]
        Probability density function values.

    References
    ----------
    Adapted from: https://github.com/blackcata/Emergent_Constraints/tree/master
    """
    work_dtype = _preferred_float_dtype(x)
    x_values = np.asarray(x, dtype=work_dtype)
    mu_value = work_dtype.type(mu)
    sigma_value = work_dtype.type(sigma)

    if sigma_value < 0:
        warnings.warn(
            "Standard deviation must be non-negative for gaussian_pdf.",
            RuntimeWarning,
            stacklevel=2,
        )
        pdf = np.full_like(x_values, np.nan, dtype=work_dtype)
    elif sigma_value == 0:
        warnings.warn(
            "Standard deviation is zero in gaussian_pdf; returning a degenerate limit.",
            RuntimeWarning,
            stacklevel=2,
        )
        pdf = np.where(
            x_values == mu_value,
            work_dtype.type(np.inf),
            work_dtype.type(0.0),
        )
    else:
        scale = work_dtype.type(1.0) / (
            np.sqrt(work_dtype.type(2.0 * np.pi)) * sigma_value
        )
        exponent = work_dtype.type(-0.5) * ((x_values - mu_value) / sigma_value) ** 2
        pdf = scale * np.exp(exponent)

    if np.ndim(x_values) == 0:
        return float(pdf)
    return pdf


def _fit_linear_relationship(
    x_models: np.ndarray, y_models: np.ndarray, work_dtype: np.dtype
) -> Tuple[np.ndarray, np.ndarray, float, float, float, float]:
    """Fit the inter-model linear relationship used by the EC formulas."""
    x_models = np.asarray(x_models, dtype=work_dtype)
    y_models = np.asarray(y_models, dtype=work_dtype)

    valid_mask = np.isfinite(x_models) & np.isfinite(y_models)
    x_models = x_models[valid_mask]
    y_models = y_models[valid_mask]

    if x_models.size != y_models.size or x_models.size == 0:
        raise ValueError(
            "constraint_data and target_data must contain valid paired data."
        )

    x_mean = np.mean(x_models)
    y_mean = np.mean(y_models)
    x_centered = x_models - x_mean
    y_centered = y_models - y_mean
    ss_x = np.sum(x_centered**2)

    if ss_x > work_dtype.type(0.0):
        slope = np.sum(x_centered * y_centered) / ss_x
    else:
        slope = work_dtype.type(0.0)

    intercept = y_mean - slope * x_mean
    residuals = y_models - (intercept + slope * x_models)

    if x_models.size > 2:
        prediction_error = np.sqrt(
            np.sum(residuals**2) / work_dtype.type(x_models.size - 2)
        )
    else:
        prediction_error = np.sqrt(np.mean(residuals**2))

    if not np.isfinite(prediction_error):
        prediction_error = work_dtype.type(0.0)

    return x_models, y_models, slope, intercept, prediction_error, ss_x


def _prediction_std(
    prediction_error: float,
    sample_size: int,
    x_values: np.ndarray,
    x_mean: float,
    ss_x: float,
) -> np.ndarray:
    """Return predictive standard deviation for regression-conditioned values."""
    work_dtype = np.dtype(np.asarray(x_values).dtype)

    if prediction_error == 0.0:
        return np.zeros_like(x_values, dtype=work_dtype)

    base = np.ones_like(x_values, dtype=work_dtype) + work_dtype.type(1.0 / sample_size)
    if ss_x > work_dtype.type(0.0):
        base = base + ((x_values - x_mean) ** 2) / ss_x

    return prediction_error * np.sqrt(base)


def emergent_constraint_posterior(
    constraint_data: xr.DataArray,
    target_data: xr.DataArray,
    constraint_grid: np.ndarray,
    target_grid: np.ndarray,
    obs_pdf: np.ndarray,
) -> Tuple[np.ndarray, float, float]:
    """
    Calculate posterior PDF using emergent constraint method.

    This function applies the emergent constraint method to reduce uncertainty
    in climate projections by utilizing observational constraints.

    Parameters
    ----------
    constraint_data : xr.DataArray
        Inter-model spread data for the constraint variable (e.g., model sensitivity).
    target_data : xr.DataArray
        Inter-model spread data for the target variable (e.g., future projection).
    constraint_grid : np.ndarray
        Grid values for the constraint variable.
    target_grid : np.ndarray
        Grid values for the target variable.
    obs_pdf : np.ndarray
        Observational PDF of the constraint variable.

    Returns
    -------
    Tuple[np.ndarray, float, float]
        A tuple containing:
        - posterior_pdf : np.ndarray - Posterior PDF of the target variable
        - posterior_std : float - Standard deviation of the target variable
        - posterior_mean : float - Mean of the target variable

    References
    ----------
    Adapted from: https://github.com/blackcata/Emergent_Constraints/tree/master
    Cox, P. M., et al. (2013). Nature, 494(7437), 341-344.
    """
    work_dtype = _preferred_float_dtype(constraint_data, target_data)
    constraint_grid = np.asarray(constraint_grid, dtype=work_dtype)
    target_grid = np.asarray(target_grid, dtype=work_dtype)
    obs_pdf = np.asarray(obs_pdf, dtype=work_dtype)
    dx = constraint_grid[1] - constraint_grid[0]

    x_models, y_models, slope, intercept, prediction_error, ss_x = (
        _fit_linear_relationship(constraint_data.values, target_data.values, work_dtype)
    )
    n_models = len(x_models)
    regression_line = intercept + slope * constraint_grid
    sigma_prediction = _prediction_std(
        prediction_error,
        n_models,
        constraint_grid,
        np.mean(x_models),
        ss_x,
    )

    posterior_pdf = np.zeros(len(target_grid), dtype=work_dtype)
    nonzero_sigma = sigma_prediction > 0.0

    if np.any(nonzero_sigma):
        likelihood = np.zeros(
            (len(target_grid), len(constraint_grid)), dtype=work_dtype
        )
        sigma_used = sigma_prediction[nonzero_sigma]
        centered = (
            target_grid[:, np.newaxis] - regression_line[nonzero_sigma][np.newaxis, :]
        )
        likelihood[:, nonzero_sigma] = np.exp(
            work_dtype.type(-0.5) * (centered / sigma_used[np.newaxis, :]) ** 2
        ) / (np.sqrt(work_dtype.type(2.0 * np.pi)) * sigma_used[np.newaxis, :])
        posterior_pdf = likelihood @ (obs_pdf * dx)
    else:
        nearest = np.argmin(
            np.abs(target_grid[:, np.newaxis] - regression_line[np.newaxis, :]), axis=0
        )
        np.add.at(posterior_pdf, nearest, obs_pdf * dx)

    # Calculate statistics
    threshold = 0.341  # For 1-sigma equivalent
    posterior_std = _calculate_std_from_pdf(threshold, target_grid, posterior_pdf)
    posterior_mean = target_grid[posterior_pdf.argmax()]

    return posterior_pdf, posterior_std, posterior_mean


def _calculate_std_from_pdf(
    threshold: float, values: np.ndarray, pdf: np.ndarray
) -> float:
    """
    Calculate standard deviation from probability density function.

    Parameters
    ----------
    threshold : float
        Threshold value for probability integration (e.g., 0.341 for 1-sigma).
    values : np.ndarray
        Grid values corresponding to the PDF.
    pdf : np.ndarray
        Probability density function values.

    Returns
    -------
    float
        Standard deviation of the distribution.

    References
    ----------
    Adapted from: https://github.com/blackcata/Emergent_Constraints/tree/master
    """
    pdf_sum = np.sum(pdf)
    if pdf_sum <= 0 or not np.isfinite(pdf_sum):
        return np.nan

    max_index = int(pdf.argmax())

    for i in range(max_index, -1, -1):
        pdf_integral = pdf[i : max_index + 1].sum() / pdf_sum
        if pdf_integral >= threshold:
            return float(abs(values[max_index] - values[i]))

    std_dev = np.sqrt(np.average((values - values[max_index]) ** 2, weights=pdf))

    return float(std_dev)


def emergent_constraint_prior(
    constraint_data: xr.DataArray,
    target_data: xr.DataArray,
    constraint_grid: np.ndarray,
    target_grid: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate prior probability distribution for emergent constraints.

    Parameters
    ----------
    constraint_data : xr.DataArray
        Inter-model spread data for the constraint variable.
    target_data : xr.DataArray
        Inter-model spread data for the target variable.
    constraint_grid : np.ndarray
        Grid values for the constraint variable.
    target_grid : np.ndarray
        Grid values for the target variable.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        A tuple containing:
        - prior_pdf : np.ndarray - Prior PDF
        - prediction_error : np.ndarray - Prediction error array
        - regression_line : np.ndarray - Linear regression values

    References
    ----------
    Adapted from: https://github.com/blackcata/Emergent_Constraints/tree/master
    """
    work_dtype = _preferred_float_dtype(constraint_data, target_data)
    constraint_grid = np.asarray(constraint_grid, dtype=work_dtype)
    target_grid = np.asarray(target_grid, dtype=work_dtype)
    x_models, y_models, slope, intercept, prediction_error_base, ss_x = (
        _fit_linear_relationship(constraint_data.values, target_data.values, work_dtype)
    )
    n_models = len(x_models)

    regression_line = intercept + slope * constraint_grid

    prediction_error = _prediction_std(
        prediction_error_base,
        n_models,
        constraint_grid,
        np.mean(x_models),
        ss_x,
    )

    prior_pdf = np.zeros((len(target_grid), len(constraint_grid)), dtype=work_dtype)
    nonzero_sigma = prediction_error > 0.0
    if np.any(nonzero_sigma):
        centered = (
            target_grid[:, np.newaxis] - regression_line[nonzero_sigma][np.newaxis, :]
        )
        sigma_used = prediction_error[nonzero_sigma]
        prior_pdf[:, nonzero_sigma] = np.exp(
            work_dtype.type(-0.5) * (centered / sigma_used[np.newaxis, :]) ** 2
        ) / (np.sqrt(work_dtype.type(2.0 * np.pi)) * sigma_used[np.newaxis, :])

    if np.any(~nonzero_sigma):
        nearest = np.argmin(
            np.abs(
                target_grid[:, np.newaxis]
                - regression_line[~nonzero_sigma][np.newaxis, :]
            ),
            axis=0,
        )
        prior_pdf[nearest, np.where(~nonzero_sigma)[0]] = 1.0

    return prior_pdf, prediction_error, regression_line


# Legacy function names for backward compatibility
def calc_GAUSSIAN_PDF(
    mu: float, sigma: float, x: Union[np.ndarray, float]
) -> Union[np.ndarray, float]:
    """Legacy function name. Use gaussian_pdf() instead."""
    return gaussian_pdf(mu, sigma, x)


def calc_PDF_EC(
    tmp_x: xr.DataArray,
    tmp_y: xr.DataArray,
    x: np.ndarray,
    y: np.ndarray,
    PDF_x: np.ndarray,
) -> Tuple[np.ndarray, float, float]:
    """Legacy function name. Use emergent_constraint_posterior() instead."""
    return emergent_constraint_posterior(tmp_x, tmp_y, x, y, PDF_x)


def find_std_from_PDF(thres: float, y: np.ndarray, PDF: np.ndarray) -> float:
    """Legacy function name. Use _calculate_std_from_pdf() instead."""
    return _calculate_std_from_pdf(thres, y, PDF)


def calc_PDF_EC_PRIOR(
    tmp_x: xr.DataArray, tmp_y: xr.DataArray, x: np.ndarray, y: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Legacy function name. Use emergent_constraint_prior() instead."""
    return emergent_constraint_prior(tmp_x, tmp_y, x, y)
