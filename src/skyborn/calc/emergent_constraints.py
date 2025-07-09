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

import numpy as np
import xarray as xr
from typing import Tuple, Union


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
    pdf = 1 / np.sqrt(2 * np.pi * sigma**2) * np.exp(-1 / 2 * ((x - mu) / sigma) ** 2)
    return pdf


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
    dx = constraint_grid[1] - constraint_grid[0]

    # Extract model data
    x_models = constraint_data.values
    y_models = target_data.values

    # Linear regression for inter-model relationship
    n_models = len(x_models)
    sigma_x = np.sqrt(np.var(x_models))
    sigma_xy = np.sqrt(np.cov(x_models, y_models))[0, 1]

    slope = (sigma_xy / sigma_x) ** 2
    intercept = -1 / n_models * (slope * x_models - y_models).sum()

    regression_line = intercept + slope * x_models

    # Prediction error
    prediction_error = np.sqrt(
        1 / (n_models - 2) * ((y_models - regression_line) ** 2).sum()
    )
    slope_uncertainty = (prediction_error / sigma_x) * np.sqrt(n_models)

    # Calculate posterior PDF
    posterior_pdf = np.zeros(len(target_grid))
    for i_target in range(len(target_grid)):
        for i_constraint in range(len(constraint_grid)):
            sigma_prediction = prediction_error * np.sqrt(
                1
                + 1 / n_models
                + (constraint_grid[i_constraint] - x_models.mean()) ** 2
                / (n_models * sigma_x**2)
            )

            likelihood = (
                1
                / np.sqrt(2 * np.pi * sigma_prediction**2)
                * np.exp(
                    -(
                        (
                            target_grid[i_target]
                            - (intercept + slope * constraint_grid[i_constraint])
                        )
                        ** 2
                    )
                    / (2 * sigma_prediction**2)
                )
            )

            posterior_pdf[i_target] += likelihood * obs_pdf[i_constraint] * dx

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
    max_index = pdf.argmax()

    for i in range(len(values)):
        pdf_integral = pdf[i : max_index + 1].sum() / pdf.sum()
        if pdf_integral < threshold:
            std_dev = values[max_index] - values[i]
            break
    else:
        # If no break occurred, use a default calculation
        std_dev = np.sqrt(np.average((values - values[max_index]) ** 2, weights=pdf))

    return std_dev


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
    x_models = constraint_data.values
    y_models = target_data.values

    n_models = len(x_models)
    sigma_x = np.sqrt(np.var(x_models))
    sigma_xy = np.sqrt(np.cov(x_models, y_models))[0, 1]

    slope = (sigma_xy / sigma_x) ** 2
    intercept = -1 / n_models * (slope * x_models - y_models).sum()

    regression_line = intercept + slope * constraint_grid

    # Prediction error
    prediction_error_base = np.sqrt(
        1 / (n_models - 2) * ((y_models - (intercept + slope * x_models)) ** 2).sum()
    )
    slope_uncertainty = (prediction_error_base / sigma_x) * np.sqrt(n_models)

    prediction_error = prediction_error_base * np.sqrt(
        1
        + 1 / n_models
        + (constraint_grid - x_models.mean()) ** 2 / (n_models * sigma_x**2)
    )

    prior_pdf = (
        1
        / np.sqrt(2 * np.pi * prediction_error**2)
        * np.exp(
            -((target_grid[:, np.newaxis] - regression_line) ** 2)
            / (2 * prediction_error**2)
        )
    )

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
