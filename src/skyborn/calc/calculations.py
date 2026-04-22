from typing import Tuple, Union

import numpy as np
import xarray as xr

# Lazy imports to avoid loading heavy dependencies at startup


def _get_metpy_calc():
    """Lazy import of metpy.calc to avoid startup overhead"""
    import metpy.calc as mpcalc

    return mpcalc


def _get_metpy_units():
    """Lazy import of metpy.units to avoid startup overhead"""
    from metpy.units import units

    return units


def _get_f_regression():
    """Lazy import of sklearn.feature_selection.f_regression"""
    from sklearn.feature_selection import f_regression

    return f_regression


def _get_pearsonr():
    """Lazy import of scipy.stats.pearsonr for p-value calculation"""
    from scipy.stats import pearsonr

    return pearsonr


__all__ = [
    "linear_regression",
    "spatial_correlation",
    "convert_longitude_range",
    "pearson_correlation",
    "spearman_correlation",
    "kendall_correlation",
    "calculate_potential_temperature",
    "calculate_dcape",
]


def linear_regression(
    data: Union[np.ndarray, xr.DataArray], predictor: Union[np.ndarray, xr.DataArray]
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform linear regression between a 3D data array and a predictor sequence.
    Handles both numpy arrays and xarray DataArrays with NaN handling.

    Parameters
    ----------
    data : np.ndarray or xr.DataArray
        A 3D array of shape (n_samples, dim1, dim2) containing dependent variables.
        Missing values should be represented as NaN.
    predictor : np.ndarray or xr.DataArray
        A 1D array of shape (n_samples,) containing the independent variable.
        Missing values should be represented as NaN.

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

    # Check for NaN values in both data and predictor
    has_nan = np.any(np.isnan(data)) or np.any(np.isnan(predictor))
    # Calculate p-values using lazy import
    f_regression = _get_f_regression()
    if has_nan:
        # Optimize np.nan access - 33% faster than repeated np.nan lookups
        undef = np.nan
        # Handle NaN case: record locations and replace with 0 in-place
        nan_mask_data = np.isnan(data)
        nan_mask_predictor = np.isnan(predictor)
        data[nan_mask_data] = 0  # Replace NaN with 0 in original array
        predictor_work = predictor.copy()
        # Replace NaN with 0 in predictor
        predictor_work[nan_mask_predictor] = 0

        # Create design matrix with predictor and intercept
        design_matrix = np.column_stack(
            (predictor_work, np.ones(predictor_work.shape[0]))
        )

        # Get original dimensions and reshape for regression
        n_samples, dim1, dim2 = data.shape
        data_flat = data.reshape((n_samples, dim1 * dim2))

        # Perform linear regression
        regression_coef, intercept = np.linalg.lstsq(
            design_matrix, data_flat, rcond=None
        )[0]
        regression_coef, intercept = regression_coef.reshape(
            (dim1, dim2)
        ), intercept.reshape((dim1, dim2))

        p_values = f_regression(data_flat, predictor_work)[1].reshape(dim1, dim2)

        # Restore original NaN values in data array
        data[nan_mask_data] = undef

        # Set results back to NaN where original data had too many NaN
        # Only consider data NaN, not predictor NaN (predictor NaN affects all gridpoints equally)
        nan_mask_gridpoint = np.any(nan_mask_data, axis=0)
        regression_coef = np.where(nan_mask_gridpoint, undef, regression_coef)
        p_values = np.where(nan_mask_gridpoint, undef, p_values)

    else:
        # No NaN case: use original efficient algorithm
        # Create design matrix with predictor and intercept
        design_matrix = np.column_stack((predictor, np.ones(predictor.shape[0])))

        # Get original dimensions and reshape for regression
        n_samples, dim1, dim2 = data.shape
        data_flat = data.reshape((n_samples, dim1 * dim2))

        # Perform linear regression
        regression_coef, _ = np.linalg.lstsq(design_matrix, data_flat, rcond=None)[0]
        regression_coef = regression_coef.reshape((dim1, dim2))

        p_values = f_regression(data_flat, predictor)[1].reshape(dim1, dim2)

    return regression_coef, p_values


def spatial_correlation(
    data: Union[np.ndarray, xr.DataArray], predictor: Union[np.ndarray, xr.DataArray]
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform fast spatial correlation between a 3D data array and a predictor time series.
    Optimized for vectorized operations to avoid slow loops over lat/lon grid points.

    Parameters
    ----------
    data : np.ndarray or xr.DataArray
        A 3D array of shape (time, lat, lon) containing spatial data.
        Missing values should be represented as NaN.
    predictor : np.ndarray or xr.DataArray
        A 1D array of shape (time,) containing the predictor time series.
        Missing values should be represented as NaN.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        A tuple containing:
        - correlation_coefficients: Pearson correlation coefficients with shape (lat, lon)
        - p_values: The p-values of the correlations with shape (lat, lon)

    Raises
    ------
    ValueError
        If the time dimension of data doesn't match the length of the predictor.
    """
    # Extract numpy arrays regardless of input type
    data = getattr(data, "values", data)
    predictor = getattr(predictor, "values", predictor)

    if len(data) != len(predictor):
        raise ValueError(
            f"Time dimension of data ({data.shape[0]}) must match "
            f"length of predictor ({len(predictor)})"
        )

    # Check for NaN values
    has_nan = np.any(np.isnan(data)) or np.any(np.isnan(predictor))

    if has_nan:
        # Optimize np.nan access - 33% faster than repeated np.nan lookups
        undef = np.nan

        # Handle NaN case: record locations and replace with 0 in-place
        nan_mask_data = np.isnan(data)
        nan_mask_predictor = np.isnan(predictor)

        # Make copies to avoid modifying original data
        data_work = data.copy()
        predictor_work = predictor.copy()

        # Replace NaN with 0 in working arrays
        data_work[nan_mask_data] = 0
        predictor_work[nan_mask_predictor] = 0

        # Get dimensions
        n_time, n_lat, n_lon = data_work.shape

        # Reshape for vectorized operations
        data_flat = data_work.reshape((n_time, n_lat * n_lon))

        # Calculate means (only for valid data points)
        # Create mask for valid data at each grid point
        valid_data_mask = ~nan_mask_data
        valid_predictor_mask = ~nan_mask_predictor

        # Calculate valid counts for each grid point
        valid_counts = np.sum(
            valid_data_mask & valid_predictor_mask[:, np.newaxis, np.newaxis], axis=0
        )

        # Only calculate correlation where we have at least 3 valid pairs
        sufficient_data = valid_counts >= 3

        # For grid points with sufficient data, calculate means
        pred_sum = np.sum(
            predictor_work[:, np.newaxis, np.newaxis]
            * (valid_data_mask & valid_predictor_mask[:, np.newaxis, np.newaxis]),
            axis=0,
        )
        data_sum = np.sum(
            data_work
            * (valid_data_mask & valid_predictor_mask[:, np.newaxis, np.newaxis]),
            axis=0,
        )

        # Avoid division by zero
        pred_means = np.divide(
            pred_sum, valid_counts, out=np.zeros_like(pred_sum), where=valid_counts > 0
        )
        data_means = np.divide(
            data_sum, valid_counts, out=np.zeros_like(data_sum), where=valid_counts > 0
        )

        # Calculate correlation using vectorized operations
        # Center the data (only for valid points)
        valid_mask_3d = (
            valid_data_mask & valid_predictor_mask[:, np.newaxis, np.newaxis]
        )

        pred_centered = (
            predictor_work[:, np.newaxis, np.newaxis] - pred_means
        ) * valid_mask_3d
        data_centered = (data_work - data_means) * valid_mask_3d

        # Vectorized correlation calculation
        numerator = np.sum(pred_centered * data_centered, axis=0)
        pred_ss = np.sum(pred_centered**2, axis=0)
        data_ss = np.sum(data_centered**2, axis=0)

        # Calculate correlation coefficients
        correlation_coef = np.divide(
            numerator,
            np.sqrt(pred_ss * data_ss),
            out=np.full((n_lat, n_lon), undef),
            where=(pred_ss > 0) & (data_ss > 0) & sufficient_data,
        )

        # Calculate p-values for valid correlations
        p_values = np.full((n_lat, n_lon), undef)
        valid_r_mask = (
            sufficient_data
            & (pred_ss > 0)
            & (data_ss > 0)
            & (np.abs(correlation_coef) < 1.0)
        )

        if np.any(valid_r_mask):
            r_valid = correlation_coef[valid_r_mask]
            n_valid = valid_counts[valid_r_mask]
            t_stat = r_valid * np.sqrt((n_valid - 2) / (1 - r_valid**2))
            from scipy.stats import t

            p_valid = 2 * (1 - t.cdf(np.abs(t_stat), n_valid - 2))
            p_values[valid_r_mask] = p_valid

        # Set p-value to 0 for perfect correlations
        perfect_r_mask = (
            sufficient_data
            & (pred_ss > 0)
            & (data_ss > 0)
            & (np.abs(correlation_coef) >= 1.0)
        )
        p_values[perfect_r_mask] = 0.0

        # Set results back to NaN where we don't have sufficient valid data
        insufficient_mask = ~sufficient_data
        correlation_coef[insufficient_mask] = undef
        p_values[insufficient_mask] = undef

    else:
        # No NaN case: use highly optimized vectorized algorithm
        n_time, n_lat, n_lon = data.shape

        # Reshape for vectorized operations
        data_flat = data.reshape((n_time, n_lat * n_lon))

        # Calculate means
        pred_mean = np.mean(predictor)
        data_means = np.mean(data_flat, axis=0)

        # Center the data
        pred_centered = predictor - pred_mean
        data_centered = data_flat - data_means[np.newaxis, :]

        # Vectorized correlation calculation across all grid points
        numerator = np.sum(pred_centered[:, np.newaxis] * data_centered, axis=0)
        pred_ss = np.sum(pred_centered**2)
        data_ss = np.sum(data_centered**2, axis=0)

        # Avoid division by zero
        valid_variance = (pred_ss > 0) & (data_ss > 0)
        correlation_coef = np.full(n_lat * n_lon, np.nan)
        correlation_coef[valid_variance] = numerator[valid_variance] / np.sqrt(
            pred_ss * data_ss[valid_variance]
        )

        # Calculate p-values vectorized
        p_values = np.full(n_lat * n_lon, np.nan)
        r_valid = correlation_coef[valid_variance]

        # Only calculate p-values where correlation is valid and not perfect
        calc_p = valid_variance & (np.abs(correlation_coef) < 1.0)
        if np.any(calc_p):
            r_calc = correlation_coef[calc_p]
            t_stat = r_calc * np.sqrt((n_time - 2) / (1 - r_calc**2))
            from scipy.stats import t

            p_values[calc_p] = 2 * (1 - t.cdf(np.abs(t_stat), n_time - 2))

        # Set p-value to 0 for perfect correlations
        p_values[valid_variance & (np.abs(correlation_coef) >= 1.0)] = 0.0

        # Reshape back to 2D
        correlation_coef = correlation_coef.reshape((n_lat, n_lon))
        p_values = p_values.reshape((n_lat, n_lon))

    return correlation_coef, p_values


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
    return data.assign_coords(
        **{
            lon: (
                lambda x: (
                    ((x[lon] + 180) % 360 - 180)
                    if not center_on_180
                    else (x[lon] % 360)
                )
            )
        }
    ).sortby(lon, ascending=True)


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
    Calculate potential temperature using fast numpy operations.

    This implementation uses lazy imports and avoids heavy metpy dependencies
    for simple potential temperature calculations.

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
        Potential temperature values in Kelvin.

    Notes
    -----
    Uses the standard formula: theta = T * (P0/P)^(R/cp)
    where R/cp = 0.286 for dry air
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


def calculate_theta_se(
    temperature: Union[np.ndarray, xr.DataArray],
    pressure: Union[np.ndarray, xr.DataArray],
    mixing_ratio: Union[np.ndarray, xr.DataArray],
    dewpoint: Union[np.ndarray, xr.DataArray],
) -> Union[np.ndarray, xr.DataArray]:
    """
    Calculate Pseudo-equivalent potential temperature (theta-se)

    the fomula is θ_se = T * (1000 / (p - e)) ** (Ra / cpd) * exp(L * r / (cpd * Tc))
    from http://stream1.cmatc.cn/cmatcvod/12/tqx/second_content.html
    Parameters
    ----------
    temperature : array-like
        Temperature in Kelvin.
    pressure : array-like
        Pressure in hPa.
    mixing_ratio : array-like
        Water vapor mixing ratio in kg/kg.
    dewpoint : array-like
        Dewpoint temperature in Celsius.

    Returns
    -------
    array-like
        Pseudo-Equivalent potential temperature in Kelvin.

    Notes
    -----
    This function uses MetPy's `vapor_pressure` and `lcl` functions internally.
    """
    # Lazy imports
    mpcalc = _get_metpy_calc()
    units = _get_metpy_units()

    # Extract values if xarray
    is_xarray = hasattr(temperature, "attrs")
    if is_xarray:
        temp_values = temperature.values
        pres_values = pressure.values
        mixr_values = mixing_ratio.values
        dewp_values = dewpoint.values
    else:
        temp_values = temperature
        pres_values = pressure
        mixr_values = mixing_ratio
        dewp_values = dewpoint

    # Convert units
    p = pres_values * units.hPa
    T = temp_values * units.kelvin
    r = mixr_values * units("kg/kg")
    td = dewp_values * units.degC

    # Convert mixing ratio to g/kg for vapor_pressure
    r_gkg = r.to("g/kg")

    # Calculate vapor pressure
    e = mpcalc.vapor_pressure(p, r_gkg)

    # Calculate LCL temperature
    _, T_lcl = mpcalc.lcl(p, T, td)

    # Convert LCL temperature to Kelvin
    T_lcl_K = T_lcl.to("kelvin").magnitude

    # Constants
    Rd = 287.0  # J/kg/K
    cp_d = 1004.0  # J/kg/K
    L = 2.5e6  # J/kg

    # Dry air pressure
    p_dry = p - e

    # Theta_e calculation
    theta_e_part = temp_values * (1000.0 / p_dry) ** (Rd / cp_d)
    latent_part = np.exp((L * mixr_values) / (cp_d * T_lcl_K))
    theta_se = theta_e_part * latent_part

    # Extract magnitude if result is Quantity
    if hasattr(theta_se, "magnitude"):
        theta_se = theta_se.magnitude

    # Preserve xarray structure if input is xarray
    if is_xarray:
        return xr.DataArray(
            theta_se,
            dims=temperature.dims,
            coords=temperature.coords,
            attrs={
                "units": "K",
                "long_name": "Pseudo-Equivalent Potential Temperature",
            },
        )
    return theta_se


# =============================================================================
# DCAPE calculation (Numba-accelerated)
# =============================================================================

try:
    from numba import njit, prange

    NUMBA_AVAILABLE = True
except Exception:  # pragma: no cover
    NUMBA_AVAILABLE = False

    def njit(*args, **kwargs):  # type: ignore
        def wrap(func):
            return func

        return wrap

    def prange(*args):  # type: ignore
        return range(*args)


_EPS, _RD, _CP, _LV = 0.622, 287.05, 1004.0, 2.5e6
_KAPPA = _RD / _CP


@njit(cache=True)
def _sat_vapor_pressure_hpa(tc):
    return 6.112 * np.exp(17.67 * tc / (tc + 243.5))


@njit(cache=True)
def _sat_vapor_pressure_pa_from_tk(tk):
    return _sat_vapor_pressure_hpa(tk - 273.15) * 100.0


@njit(cache=True)
def _mixing_ratio_from_e_pa(p_pa, e_pa):
    e = max(e_pa, 1e-8)
    e = min(e, 0.99 * p_pa)
    return _EPS * e / (p_pa - e)


@njit(cache=True)
def _virtual_temp_k(tk, mixing_ratio):
    return tk * (mixing_ratio + _EPS) / (_EPS * (1.0 + mixing_ratio))


@njit(cache=True)
def _thetae_bolton_like_metpy(p_hpa, tc, tdc):
    t = tc + 273.15
    td = max(tdc + 273.15, 173.15)

    p_pa = p_hpa * 100.0
    e_pa = _sat_vapor_pressure_hpa(tdc) * 100.0
    r = _mixing_ratio_from_e_pa(p_pa, e_pa)

    t_l = 56.0 + 1.0 / (1.0 / (td - 56.0) + np.log(t / td) / 800.0)
    th_l = t * (100000.0 / (p_pa - e_pa)) ** _KAPPA * (t / t_l) ** (0.28 * r)
    return th_l * np.exp(r * (1.0 + 0.448 * r) * (3036.0 / t_l - 1.78))


@njit(cache=True)
def _dt_dp_moist(p_pa, t_k):
    es = _sat_vapor_pressure_pa_from_tk(t_k)
    rs = _mixing_ratio_from_e_pa(p_pa, es)
    num = _RD * t_k + _LV * rs
    den = _CP + (_LV * _LV * rs * _EPS) / (_RD * t_k * t_k)
    return (num / den) / p_pa


@njit(cache=True)
def _rk4_integrate_moist_t(p_start_pa, t_start_k, p_target_pa):
    dp_total = p_target_pa - p_start_pa
    if np.abs(dp_total) < 1e-9:
        return t_start_k

    nstep = int(np.abs(dp_total) / 200.0) + 1
    h = dp_total / nstep
    p, t = p_start_pa, t_start_k

    for _ in range(nstep):
        k1 = _dt_dp_moist(p, t)
        k2 = _dt_dp_moist(p + 0.5 * h, t + 0.5 * h * k1)
        k3 = _dt_dp_moist(p + 0.5 * h, t + 0.5 * h * k2)
        k4 = _dt_dp_moist(p + h, t + h * k3)
        t += (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        p += h
    return t


@njit(cache=True)
def _lcl_bolton(p_hpa, tc, tdc):
    t = tc + 273.15
    td = max(tdc + 273.15, 173.15)
    tlcl = 56.0 + 1.0 / (1.0 / (td - 56.0) + np.log(t / td) / 800.0)
    plcl = p_hpa * (tlcl / t) ** (1.0 / _KAPPA)
    return plcl, tlcl


@njit(cache=True)
def _interp_linear_p(p0, v0, p1, v1, p_target):
    if np.abs(p1 - p0) < 1e-12:
        return v0
    if p0 > 0.0 and p1 > 0.0 and p_target > 0.0:
        lp0, lp1, lpt = np.log(p0), np.log(p1), np.log(p_target)
        if np.abs(lp1 - lp0) < 1e-12:
            return v0
        w = (lpt - lp0) / (lp1 - lp0)
    else:
        w = (p_target - p0) / (p1 - p0)
    return v0 + w * (v1 - v0)


@njit(cache=True)
def _thetae_min_in_700_500(p, t, td, m):
    min_te = 1.0e99
    start_p = np.nan
    start_t = np.nan
    start_td = np.nan
    found = False

    for k in range(m):
        pk = p[k]
        if 500.0 <= pk <= 700.0:
            te = _thetae_bolton_like_metpy(pk, t[k], td[k])
            if te < min_te:
                min_te, start_p, start_t, start_td, found = te, pk, t[k], td[k], True

    for bound in (700.0, 500.0):
        for k in range(m - 1):
            p0, p1 = p[k], p[k + 1]
            if (p0 >= bound >= p1) or (p0 <= bound <= p1):
                tb = _interp_linear_p(p0, t[k], p1, t[k + 1], bound)
                tdb = _interp_linear_p(p0, td[k], p1, td[k + 1], bound)
                te = _thetae_bolton_like_metpy(bound, tb, tdb)
                if te < min_te:
                    min_te, start_p, start_t, start_td, found = te, bound, tb, tdb, True
                break

    return found, start_p, start_t, start_td


@njit(cache=True)
def _compute_dcape_profile_numba(pressure_hpa, temperature_c, dewpoint_c):
    n = pressure_hpa.size
    p = np.empty(n, dtype=np.float64)
    t = np.empty(n, dtype=np.float64)
    td = np.empty(n, dtype=np.float64)

    m = 0
    for k in range(n):
        pk, tk, tdk = pressure_hpa[k], temperature_c[k], dewpoint_c[k]
        if np.isfinite(pk) and np.isfinite(tk) and np.isfinite(tdk):
            p[m], t[m], td[m] = pk, tk, tdk
            m += 1
    if m < 4:
        return np.nan

    if p[0] < p[m - 1]:
        for k in range(m // 2):
            k2 = m - 1 - k
            tmp = p[k]
            p[k] = p[k2]
            p[k2] = tmp
            tmp = t[k]
            t[k] = t[k2]
            t[k2] = tmp
            tmp = td[k]
            td[k] = td[k2]
            td[k2] = tmp

    found, start_p_hpa, start_t_c, start_td_c = _thetae_min_in_700_500(p, t, td, m)
    if not found:
        return np.nan

    start_idx = -1
    for k in range(m):
        if p[k] >= start_p_hpa:
            start_idx = k
        else:
            break
    if start_idx <= 0:
        return np.nan

    plcl_hpa, tlcl_k = _lcl_bolton(start_p_hpa, start_t_c, start_td_c)
    tw_start_k = _rk4_integrate_moist_t(plcl_hpa * 100.0, tlcl_k, start_p_hpa * 100.0)

    parcel_tk = np.empty(start_idx + 1, dtype=np.float64)
    parcel_tk[start_idx] = _rk4_integrate_moist_t(
        start_p_hpa * 100.0, tw_start_k, p[start_idx] * 100.0
    )
    for k in range(start_idx, 0, -1):
        parcel_tk[k - 1] = _rk4_integrate_moist_t(
            p[k] * 100.0, parcel_tk[k], p[k - 1] * 100.0
        )

    integral = 0.0
    prev_diff = 0.0
    prev_lnp = 0.0
    first = True

    for k in range(start_idx + 1):
        p_pa = p[k] * 100.0
        e_env = _sat_vapor_pressure_hpa(td[k]) * 100.0
        r_env = _mixing_ratio_from_e_pa(p_pa, e_env)
        tv_env = _virtual_temp_k(t[k] + 273.15, r_env)

        e_par = _sat_vapor_pressure_pa_from_tk(parcel_tk[k])
        r_par = _mixing_ratio_from_e_pa(p_pa, e_par)
        tv_par = _virtual_temp_k(parcel_tk[k], r_par)

        diff = tv_env - tv_par
        lnp = np.log(p[k])
        if first:
            prev_diff, prev_lnp, first = diff, lnp, False
        else:
            integral += 0.5 * (prev_diff + diff) * (lnp - prev_lnp)
            prev_diff, prev_lnp = diff, lnp

    dcape = -_RD * integral
    if not np.isfinite(dcape):
        return np.nan
    return 0.0 if dcape < 0.0 else dcape


@njit(cache=True, parallel=True)
def _compute_dcape_grid_numba(pressure_3d, t_3d, td_3d, j_map, i_map):
    nj, ni = j_map.size, i_map.size
    out = np.empty((nj, ni), dtype=np.float32)

    for sj in prange(nj):
        j = j_map[sj]
        for si in range(ni):
            i = i_map[si]
            out[sj, si] = _compute_dcape_profile_numba(
                pressure_3d[:, j, i], t_3d[:, j, i], td_3d[:, j, i]
            )

    return out


def calculate_dcape(
    pressure: Union[np.ndarray, xr.DataArray],
    temperature: Union[np.ndarray, xr.DataArray],
    dewpoint: Union[np.ndarray, xr.DataArray],
) -> Union[np.ndarray, xr.DataArray, float]:
    """
    Calculate Downdraft Convective Available Potential Energy (DCAPE).

    This implementation uses a Numba-accelerated Bolton-like moist-adiabatic
    integration with RK4. It finds the minimum theta-e layer between 700 and
    500 hPa, computes the wet-bulb temperature along a pseudoadiabat down to
    the surface, and integrates the negative buoyancy.

    Parameters
    ----------
    pressure : np.ndarray or xr.DataArray
        Pressure in hPa. Can be 1D (level,) for a single profile, or 3D
        (level, lat, lon) for a spatial grid.
    temperature : np.ndarray or xr.DataArray
        Temperature in Celsius. Same shape as ``pressure``.
    dewpoint : np.ndarray or xr.DataArray
        Dewpoint temperature in Celsius. Same shape as ``pressure``.

    Returns
    -------
    np.ndarray or xr.DataArray or float
        DCAPE in J kg^-1. For 3D input returns a 2D array (lat, lon).
        For 1D input returns a scalar float.

    Notes
    -----
    - If numba is not available the calculation still works but will be much
      slower for grid data.
    - Profiles with fewer than 4 valid levels, or without a detectable
      theta-e minimum between 500--700 hPa, return NaN.
    """
    # Extract numpy arrays
    p = getattr(pressure, "values", pressure)
    t = getattr(temperature, "values", temperature)
    td = getattr(dewpoint, "values", dewpoint)

    is_xarray = hasattr(pressure, "attrs")

    # 1D profile case
    if p.ndim == 1:
        result = _compute_dcape_profile_numba(
            p.astype(np.float64),
            t.astype(np.float64),
            td.astype(np.float64),
        )
        return float(result)

    # 3D grid case: (level, lat, lon)
    if p.ndim == 3:
        nz, ny, nx = p.shape
        j_idx = np.arange(ny, dtype=np.int64)
        i_idx = np.arange(nx, dtype=np.int64)

        out = _compute_dcape_grid_numba(
            p.astype(np.float64),
            t.astype(np.float64),
            td.astype(np.float64),
            j_idx,
            i_idx,
        )

        if is_xarray:
            # Attempt to preserve spatial coordinates from the input
            coords = {}
            dims = []
            if hasattr(pressure, "coords"):
                coord_names = list(pressure.coords)
                # Assume last two dims are the spatial ones
                spatial_dims = pressure.dims[-2:]
                for d in spatial_dims:
                    if d in pressure.coords:
                        coords[d] = pressure.coords[d]
                dims = list(spatial_dims)
            return xr.DataArray(
                out,
                dims=dims or ["y", "x"],
                coords=coords,
                attrs={"units": "J kg-1", "long_name": "DCAPE"},
            )
        return out

    raise ValueError(
        f"pressure must be 1D or 3D, got shape {p.shape}"
    )

