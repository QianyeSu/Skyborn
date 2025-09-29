"""
Xarray interface for Standardized Precipitation Index calculation.

This module provides xarray-compatible functions for calculating SPI
with automatic dimension handling and metadata preservation.
"""

import numpy as np
import xarray as xr
from typing import Union, Optional, Dict, Any
import warnings

from .core import standardized_precipitation_index


def spi_xarray(
    precipitation: xr.DataArray,
    time_scale: int = 3,
    time_dim: Optional[str] = None,
    distribution: str = 'gamma',
    **kwargs
) -> xr.DataArray:
    """
    Calculate Standardized Precipitation Index for xarray DataArrays.
    
    This function provides a convenient interface for calculating SPI on
    xarray DataArrays with automatic dimension detection and metadata
    preservation.
    
    Parameters
    ----------
    precipitation : xr.DataArray
        Precipitation data as xarray DataArray. Should have a time dimension.
    time_scale : int, default 3
        Time scale for SPI calculation in time units (1, 3, 6, 12, etc.)
    time_dim : str, optional
        Name of the time dimension. If None, will attempt to detect automatically.
    distribution : str, default 'gamma'
        Distribution to fit to precipitation data. Currently only 'gamma' is supported.
    **kwargs
        Additional keyword arguments passed to core SPI function
        
    Returns
    -------
    xr.DataArray
        SPI values with same coordinates and dimensions as input, 
        with updated attributes describing the calculation
        
    Raises
    ------
    ValueError
        If time dimension cannot be found or identified
        
    Examples
    --------
    >>> import xarray as xr
    >>> import numpy as np
    >>> from skyborn.calc.spi import spi_xarray
    
    # Create sample precipitation data
    >>> time = pd.date_range('2000-01-01', periods=120, freq='M')
    >>> lat = np.linspace(-30, 30, 10) 
    >>> lon = np.linspace(0, 350, 15)
    >>> precip_data = np.random.gamma(2, 2, size=(120, 10, 15))
    >>> precip = xr.DataArray(
    ...     precip_data,
    ...     coords={'time': time, 'lat': lat, 'lon': lon},
    ...     dims=['time', 'lat', 'lon'],
    ...     attrs={'units': 'mm', 'long_name': 'precipitation'}
    ... )
    
    # Calculate 3-month SPI
    >>> spi_3m = spi_xarray(precip, time_scale=3)
    >>> print(spi_3m.attrs['long_name'])
    '3-month Standardized Precipitation Index'
    
    # Calculate 12-month SPI 
    >>> spi_12m = spi_xarray(precip, time_scale=12)
    """
    
    # Validate input
    if not isinstance(precipitation, xr.DataArray):
        raise TypeError("precipitation must be an xarray DataArray")
    
    # Detect time dimension
    if time_dim is None:
        time_dim = _detect_time_dimension(precipitation)
    
    if time_dim not in precipitation.dims:
        raise ValueError(f"Time dimension '{time_dim}' not found in DataArray dimensions: {precipitation.dims}")
    
    # Get the axis number for the time dimension
    time_axis = precipitation.get_axis_num(time_dim)
    
    # Extract data and calculate SPI
    precip_data = precipitation.values
    
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        spi_data = standardized_precipitation_index(
            precip_data, 
            time_scale=time_scale,
            axis=time_axis,
            distribution=distribution,
            **kwargs
        )
    
    # Create result DataArray with same coordinates and dimensions
    spi_result = xr.DataArray(
        spi_data,
        coords=precipitation.coords,
        dims=precipitation.dims,
        name=f'spi_{time_scale}m'
    )
    
    # Update attributes
    attrs = _create_spi_attributes(precipitation.attrs, time_scale, distribution)
    spi_result.attrs.update(attrs)
    
    return spi_result


def spi_dataset(
    dataset: xr.Dataset, 
    precip_var: str,
    time_scales: Union[int, list] = [1, 3, 6, 12],
    time_dim: Optional[str] = None,
    distribution: str = 'gamma',
    **kwargs
) -> xr.Dataset:
    """
    Calculate SPI for multiple time scales and return as Dataset.
    
    Parameters
    ----------
    dataset : xr.Dataset
        Dataset containing precipitation data
    precip_var : str
        Name of precipitation variable in dataset
    time_scales : int or list of int, default [1, 3, 6, 12]
        Time scale(s) for SPI calculation
    time_dim : str, optional
        Name of time dimension. If None, will attempt to detect automatically.
    distribution : str, default 'gamma'
        Distribution to fit to precipitation data
    **kwargs
        Additional keyword arguments passed to SPI calculation
        
    Returns
    -------
    xr.Dataset
        Dataset containing SPI variables for each time scale
        
    Examples
    --------
    >>> # Calculate multiple SPI time scales
    >>> spi_ds = spi_dataset(dataset, 'precipitation', time_scales=[3, 6, 12])
    >>> print(list(spi_ds.data_vars))
    ['spi_3m', 'spi_6m', 'spi_12m']
    """
    
    if precip_var not in dataset:
        raise ValueError(f"Precipitation variable '{precip_var}' not found in dataset")
    
    precipitation = dataset[precip_var]
    
    # Ensure time_scales is a list
    if isinstance(time_scales, int):
        time_scales = [time_scales]
    
    # Calculate SPI for each time scale
    spi_vars = {}
    for ts in time_scales:
        spi_result = spi_xarray(
            precipitation, 
            time_scale=ts,
            time_dim=time_dim,
            distribution=distribution,
            **kwargs
        )
        spi_vars[f'spi_{ts}m'] = spi_result
    
    # Create new dataset with SPI variables
    result_ds = xr.Dataset(spi_vars)
    
    # Copy coordinates from original dataset
    for coord_name, coord_data in dataset.coords.items():
        if coord_name not in result_ds.coords:
            result_ds.coords[coord_name] = coord_data
    
    # Update global attributes
    result_ds.attrs.update(dataset.attrs)
    result_ds.attrs['spi_calculation'] = f'SPI calculated for time scales: {time_scales}'
    result_ds.attrs['spi_distribution'] = distribution
    
    return result_ds


def _detect_time_dimension(da: xr.DataArray) -> str:
    """
    Attempt to automatically detect time dimension in DataArray.
    
    Parameters
    ----------
    da : xr.DataArray
        Input DataArray
        
    Returns
    -------
    str
        Name of detected time dimension
        
    Raises
    ------
    ValueError
        If time dimension cannot be detected
    """
    
    # Common time dimension names
    time_names = ['time', 'Time', 'TIME', 't', 'date', 'Date']
    
    # Check for exact matches first
    for name in time_names:
        if name in da.dims:
            return name
    
    # Check coordinate types
    for dim_name in da.dims:
        if dim_name in da.coords:
            coord = da.coords[dim_name]
            
            # Check if coordinate has datetime-like dtype
            if np.issubdtype(coord.dtype, np.datetime64):
                return dim_name
            
            # Check for time-related attributes
            if hasattr(coord, 'attrs'):
                long_name = coord.attrs.get('long_name', '').lower()
                standard_name = coord.attrs.get('standard_name', '').lower()
                
                if any(keyword in long_name or keyword in standard_name 
                       for keyword in ['time', 'date']):
                    return dim_name
    
    # If no clear time dimension found, use first dimension as fallback
    if da.dims:
        warnings.warn(
            f"Could not detect time dimension. Using first dimension '{da.dims[0]}' as time axis. "
            f"Specify time_dim parameter explicitly if this is incorrect.",
            UserWarning
        )
        return da.dims[0]
    
    raise ValueError("Cannot detect time dimension in DataArray")


def _create_spi_attributes(
    original_attrs: Dict[str, Any], 
    time_scale: int, 
    distribution: str
) -> Dict[str, Any]:
    """
    Create appropriate attributes for SPI DataArray.
    
    Parameters
    ----------
    original_attrs : dict
        Original attributes from precipitation data
    time_scale : int
        Time scale used for SPI calculation
    distribution : str
        Distribution used for fitting
        
    Returns
    -------
    dict
        Attributes for SPI DataArray
    """
    
    attrs = {
        'long_name': f'{time_scale}-month Standardized Precipitation Index',
        'standard_name': 'standardized_precipitation_index',
        'units': '1',  # Dimensionless
        'description': (
            f'Standardized Precipitation Index calculated over {time_scale}-month time scale. '
            f'Based on {distribution} distribution fitting to precipitation data.'
        ),
        'spi_time_scale': time_scale,
        'spi_distribution': distribution,
        'interpretation': (
            'SPI >= 2.0: Extremely wet; 1.5 <= SPI < 2.0: Very wet; '
            '1.0 <= SPI < 1.5: Moderately wet; -1.0 < SPI < 1.0: Near normal; '
            '-1.5 < SPI <= -1.0: Moderately dry; -2.0 < SPI <= -1.5: Severely dry; '
            'SPI <= -2.0: Extremely dry'
        ),
        'references': (
            'McKee, T.B., Doesken, N.J. and Kleist, J., 1993. The relationship of drought '
            'frequency and duration to time scales. In Proceedings of the 8th Conference on '
            'Applied Climatology (Vol. 17, No. 22, pp. 179-183).'
        )
    }
    
    # Preserve some original attributes if relevant
    if 'source' in original_attrs:
        attrs['source'] = original_attrs['source']
    if 'history' in original_attrs:
        attrs['history'] = original_attrs['history'] + f'; SPI calculated with time_scale={time_scale}'
    else:
        attrs['history'] = f'SPI calculated with time_scale={time_scale}'
    
    return attrs