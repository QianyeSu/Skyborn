import numpy as np
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
from typing import Tuple, Union
from sklearn.feature_selection import f_regression


def linear_regression(
    data: Union[np.ndarray, xr.DataArray], predictor: Union[np.ndarray, xr.DataArray]
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform linear regression between a 3D data array and a predictor sequence.
    Handles both numpy arrays and xarray DataArrays.

    Parameters
    ----------
    data : np.ndarray or xr.DataArray
        A 3D array of shape (n_samples, dim1, dim2) containing dependent variables.
    predictor : np.ndarray or xr.DataArray
        A 1D array of shape (n_samples,) containing the independent variable.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        A tuple containing:
        - regression_coefficients: The slope of the regression line with shape (dim1, dim2)
        - p_values: The p-values of the regression with shape (dim1, dim2)

    Raises
    ------
    ValueError
        If the number of samples in data doesn't match the length of the predictor.
    """
    # Extract numpy arrays regardless of input type
    data = getattr(data, "values", data)
    predictor = getattr(predictor, "values", predictor)

    if len(data) != len(predictor):
        raise ValueError(
            f"Number of samples in data ({data.shape[0]}) must match "
            f"length of predictor ({len(predictor)})"
        )

    # Create design matrix with predictor and intercept
    design_matrix = np.column_stack((predictor, np.ones(predictor.shape[0])))

    # Get original dimensions and reshape for regression
    n_samples, dim1, dim2 = data.shape
    data_flat = data.reshape((n_samples, dim1 * dim2))

    # Perform linear regression
    regression_coef = np.linalg.lstsq(design_matrix, data_flat, rcond=None)[0][0]

    # Reshape results back to original dimensions
    regression_coef = regression_coef.reshape((dim1, dim2))

    # Calculate p-values
    p_values = f_regression(data_flat, predictor)[1].reshape(dim1, dim2)

    return regression_coef, p_values


def convert_longitude_range(
    data: Union[xr.DataArray, xr.Dataset], lon: str = "lon", center_on_180: bool = True
) -> Union[xr.DataArray, xr.Dataset]:
    """
    Wrap longitude coordinates of DataArray or Dataset to either -180..179 or 0..359.

    Parameters
    ----------
    data : xr.DataArray or xr.Dataset
        An xarray DataArray or Dataset object containing longitude coordinates.
    lon : str, optional
        The name of the longitude coordinate, default is 'lon'.
    center_on_180 : bool, optional
        If True, wrap longitude from 0..359 to -180..179;
        If False, wrap longitude from -180..179 to 0..359.

    Returns
    -------
    xr.DataArray or xr.Dataset
        The DataArray or Dataset with wrapped longitude coordinates.
    """
    # Wrap -180..179 to 0..359
    if center_on_180:
        data = data.assign_coords(**{lon: (lambda x: (x[lon] % 360))})
    # Wrap 0..359 to -180..179
    else:
        data = data.assign_coords(**{lon: (lambda x: ((x[lon] + 180) % 360) - 180)})
    return data.sortby(lon, ascending=True)


def pearson_correlation(
    x: Union[np.ndarray, xr.DataArray], y: Union[np.ndarray, xr.DataArray]
) -> float:
    """
    Calculate Pearson correlation coefficient.

    Parameters
    ----------
    x, y : array-like
        Input data arrays.

    Returns
    -------
    float
        Pearson correlation coefficient.
    """
    x = getattr(x, "values", x)
    y = getattr(y, "values", y)
    return np.corrcoef(x.flatten(), y.flatten())[0, 1]


def spearman_correlation(
    x: Union[np.ndarray, xr.DataArray], y: Union[np.ndarray, xr.DataArray]
) -> float:
    """
    Calculate Spearman rank correlation coefficient.

    Parameters
    ----------
    x, y : array-like
        Input data arrays.

    Returns
    -------
    float
        Spearman correlation coefficient.
    """
    from scipy.stats import spearmanr

    x = getattr(x, "values", x)
    y = getattr(y, "values", y)
    correlation, _ = spearmanr(x.flatten(), y.flatten())
    return correlation


def kendall_correlation(
    x: Union[np.ndarray, xr.DataArray], y: Union[np.ndarray, xr.DataArray]
) -> float:
    """
    Calculate Kendall's tau correlation coefficient.

    Parameters
    ----------
    x, y : array-like
        Input data arrays.

    Returns
    -------
    float
        Kendall's tau correlation coefficient.
    """
    from scipy.stats import kendalltau

    x = getattr(x, "values", x)
    y = getattr(y, "values", y)
    correlation, _ = kendalltau(x.flatten(), y.flatten())
    return correlation


def calculate_potential_temperature(
    temperature: Union[np.ndarray, xr.DataArray],
    pressure: Union[np.ndarray, xr.DataArray],
    reference_pressure: float = 1000.0,
) -> Union[np.ndarray, xr.DataArray]:
    """
    Calculate potential temperature.

    Parameters
    ----------
    temperature : array-like
        Temperature values in Kelvin.
    pressure : array-like
        Pressure values in hPa.
    reference_pressure : float, optional
        Reference pressure in hPa, default is 1000.0.

    Returns
    -------
    array-like
        Potential temperature values.
    """
    R_over_cp = 0.286  # R/cp for dry air
    potential_temp = temperature * (reference_pressure / pressure) ** R_over_cp

    if hasattr(temperature, "attrs"):
        if isinstance(potential_temp, np.ndarray):
            return xr.DataArray(
                potential_temp,
                attrs={"units": "K", "long_name": "Potential Temperature"},
            )
        else:
            potential_temp.attrs = {"units": "K", "long_name": "Potential Temperature"}

    return potential_temp
