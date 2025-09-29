"""
Core SPI calculation functions.

This module provides the mathematical implementation of the Standardized
Precipitation Index using vectorized operations for efficient computation
on multi-dimensional datasets.
"""

import numpy as np
from scipy import stats
from typing import Union, Optional, Tuple
import warnings


def _gamma_fit_vectorized(data: np.ndarray, axis: int = 0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit gamma distribution parameters to precipitation data along specified axis.
    
    Parameters
    ----------
    data : np.ndarray
        Precipitation data with shape (..., time, ...)
    axis : int, default 0
        Axis along which to fit the distribution (typically time axis)
        
    Returns
    -------
    alpha : np.ndarray
        Shape parameter of gamma distribution
    beta : np.ndarray  
        Scale parameter of gamma distribution
    """
    # Move time axis to the end for easier processing
    data_moved = np.moveaxis(data, axis, -1)
    original_shape = data_moved.shape
    
    # Reshape to (spatial_points, time)
    data_reshaped = data_moved.reshape(-1, original_shape[-1])
    
    # Initialize output arrays
    n_points = data_reshaped.shape[0]
    alpha = np.full(n_points, np.nan)
    beta = np.full(n_points, np.nan)
    
    # Fit gamma distribution for each spatial point
    for i in range(n_points):
        series = data_reshaped[i, :]
        
        # Remove zeros and NaNs for gamma fitting
        valid_data = series[~np.isnan(series) & (series > 0)]
        
        if len(valid_data) < 10:  # Need sufficient data for fitting
            continue
            
        try:
            # Fit gamma distribution using method of moments as initial guess
            mean_val = np.mean(valid_data)
            var_val = np.var(valid_data)
            
            if var_val > 0 and mean_val > 0:
                # Method of moments estimates
                alpha_est = mean_val**2 / var_val
                beta_est = var_val / mean_val
                
                # Use scipy's gamma fit with good initial guess
                alpha_fit, _, beta_fit = stats.gamma.fit(valid_data, fa=alpha_est, scale=beta_est)
                
                alpha[i] = alpha_fit
                beta[i] = beta_fit
                
        except (RuntimeError, ValueError):
            # If fitting fails, use method of moments
            try:
                mean_val = np.mean(valid_data) 
                var_val = np.var(valid_data)
                if var_val > 0 and mean_val > 0:
                    alpha[i] = mean_val**2 / var_val
                    beta[i] = var_val / mean_val
            except:
                continue
    
    # Reshape back to original spatial dimensions
    spatial_shape = original_shape[:-1]
    alpha = alpha.reshape(spatial_shape)
    beta = beta.reshape(spatial_shape)
    
    return alpha, beta


def _calculate_spi_values(precip: np.ndarray, alpha: np.ndarray, beta: np.ndarray, 
                         axis: int = 0) -> np.ndarray:
    """
    Calculate SPI values using fitted gamma parameters.
    
    Parameters
    ----------
    precip : np.ndarray
        Precipitation data
    alpha : np.ndarray
        Shape parameter of gamma distribution  
    beta : np.ndarray
        Scale parameter of gamma distribution
    axis : int, default 0
        Time axis
        
    Returns
    -------
    spi : np.ndarray
        Standardized Precipitation Index values
    """
    # Move time axis to the end
    precip_moved = np.moveaxis(precip, axis, -1)
    original_shape = precip_moved.shape
    
    # Reshape for vectorized operations
    precip_reshaped = precip_moved.reshape(-1, original_shape[-1])
    alpha_flat = alpha.flatten()
    beta_flat = beta.flatten()
    
    # Initialize output
    spi_reshaped = np.full_like(precip_reshaped, np.nan)
    
    for i in range(precip_reshaped.shape[0]):
        if np.isnan(alpha_flat[i]) or np.isnan(beta_flat[i]):
            continue
            
        series = precip_reshaped[i, :]
        
        # Calculate cumulative probability using gamma distribution
        # Handle zeros separately
        prob = np.full_like(series, np.nan)
        
        # For zero precipitation values
        zero_mask = (series == 0) & ~np.isnan(series)
        nonzero_mask = (series > 0) & ~np.isnan(series)
        
        if np.any(nonzero_mask):
            # Probability for non-zero values
            prob[nonzero_mask] = stats.gamma.cdf(series[nonzero_mask], 
                                                a=alpha_flat[i], scale=beta_flat[i])
        
        if np.any(zero_mask):
            # For zeros, use the probability of zero precipitation
            # This is often estimated as the fraction of zero values in the dataset
            zero_count = np.sum(series == 0)
            total_count = np.sum(~np.isnan(series))
            if total_count > 0:
                prob_zero = zero_count / total_count
                prob[zero_mask] = prob_zero
        
        # Convert to standard normal distribution
        # Ensure probabilities are in valid range (0, 1)
        prob = np.clip(prob, 1e-6, 1-1e-6)
        spi_reshaped[i, :] = stats.norm.ppf(prob)
    
    # Reshape back to original shape
    spi = spi_reshaped.reshape(original_shape)
    
    # Move time axis back to original position
    spi = np.moveaxis(spi, -1, axis)
    
    return spi


def _rolling_sum(data: np.ndarray, window: int, axis: int = 0) -> np.ndarray:
    """
    Calculate rolling sum along specified axis.
    
    Parameters
    ----------
    data : np.ndarray
        Input data
    window : int
        Rolling window size
    axis : int, default 0
        Axis along which to calculate rolling sum
        
    Returns
    -------
    np.ndarray
        Rolling sum values (same shape as input)
    """
    if window == 1:
        return data.copy()
    
    # Move target axis to the end
    data_moved = np.moveaxis(data, axis, -1)
    original_shape = data_moved.shape
    
    # Reshape to (spatial_points, time)
    reshaped = data_moved.reshape(-1, original_shape[-1])
    
    # Initialize result with NaN
    result = np.full(reshaped.shape, np.nan, dtype=np.float64)
    
    for i in range(reshaped.shape[0]):
        series = reshaped[i, :].astype(np.float64)  # Ensure float type
        
        # Calculate rolling sum for each position
        for j in range(len(series)):
            start_idx = max(0, j - window + 1)
            end_idx = j + 1
            window_data = series[start_idx:end_idx]
            
            # Only calculate sum if we have complete window and all values are valid
            if len(window_data) == window and np.all(np.isfinite(window_data)):
                result[i, j] = np.sum(window_data)
    
    # Reshape back to original shape and move axis back
    result = result.reshape(original_shape)
    result = np.moveaxis(result, -1, axis)
    
    return result


def standardized_precipitation_index(
    precipitation: np.ndarray,
    time_scale: int = 3,
    axis: int = 0,
    distribution: str = 'gamma'
) -> np.ndarray:
    """
    Calculate the Standardized Precipitation Index (SPI).
    
    The SPI is a widely used index to characterize meteorological drought
    on a range of time scales. This implementation provides efficient
    calculation for multi-dimensional datasets.
    
    Parameters
    ----------
    precipitation : np.ndarray
        Precipitation data. Can be multi-dimensional with time along any axis.
    time_scale : int, default 3
        Time scale for SPI calculation in months (1, 3, 6, 12, etc.)
    axis : int, default 0
        Axis along which time varies
    distribution : str, default 'gamma'
        Distribution to fit to precipitation data. Currently only 'gamma' is supported.
        
    Returns
    -------
    np.ndarray
        Standardized Precipitation Index values with same shape as input
        
    Notes
    -----
    The SPI calculation involves:
    1. Aggregating precipitation over the specified time scale
    2. Fitting a probability distribution (gamma) to the aggregated data
    3. Transforming to standard normal distribution
    
    SPI values interpretation:
    - SPI ≥ 2.0: Extremely wet
    - 1.5 ≤ SPI < 2.0: Very wet  
    - 1.0 ≤ SPI < 1.5: Moderately wet
    - -1.0 < SPI < 1.0: Near normal
    - -1.5 < SPI ≤ -1.0: Moderately dry
    - -2.0 < SPI ≤ -1.5: Severely dry
    - SPI ≤ -2.0: Extremely dry
    
    Examples
    --------
    >>> import numpy as np
    >>> from skyborn.calc.spi import standardized_precipitation_index
    
    # Generate sample precipitation data (time, lat, lon)
    >>> precip = np.random.gamma(2, 2, size=(120, 10, 15))  # 10 years monthly data
    >>> spi_3m = standardized_precipitation_index(precip, time_scale=3, axis=0)
    >>> print(spi_3m.shape)
    (120, 10, 15)
    """
    if distribution != 'gamma':
        raise ValueError("Currently only 'gamma' distribution is supported")
    
    if time_scale < 1:
        raise ValueError("Time scale must be >= 1")
    
    # Step 1: Calculate rolling sum for the specified time scale
    if time_scale > 1:
        precip_aggregated = _rolling_sum(precipitation, time_scale, axis=axis)
    else:
        precip_aggregated = precipitation.copy()
    
    # Step 2: Fit gamma distribution parameters
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        alpha, beta = _gamma_fit_vectorized(precip_aggregated, axis=axis)
    
    # Step 3: Calculate SPI values
    spi_values = _calculate_spi_values(precip_aggregated, alpha, beta, axis=axis)
    
    return spi_values


# Convenience alias
spi = standardized_precipitation_index