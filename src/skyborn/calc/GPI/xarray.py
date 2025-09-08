"""
XArray interface for Genesis Potential Index (GPI) and Tropical Cyclone Potential Intensity calculations.

This module provides xarray-native interfaces for tropical cyclone potential intensity
calculations, with automatic coordinate handling and metadata preservation.
"""

import warnings
from typing import Optional, Tuple, Union

import numpy as np
import xarray as xr

from .interface import (
    calculate_potential_intensity_3d,
    calculate_potential_intensity_4d,
    calculate_potential_intensity_profile,
)


def _detect_atmospheric_dimensions(
    data_array: xr.DataArray,
) -> Tuple[Optional[int], Optional[int], int, Optional[int]]:
    """
    Auto-detect dimension indices for atmospheric data in xarray DataArray.
    Supports 1D, 3D, and 4D data.

    Parameters
    ----------
    data_array : xr.DataArray
        Atmospheric data to analyze

    Returns
    -------
    xdim : int, optional
        Longitude dimension index (None for 1D data without longitude)
    ydim : int, optional
        Latitude dimension index (None for 1D data without latitude)
    levdim : int
        Vertical level dimension index
    timedim : int, optional
        Time dimension index (None if not found)

    Raises
    ------
    ValueError
        If level dimension cannot be identified
    """
    dims = data_array.dims

    # Common dimension name patterns
    lon_names = {
        "lon",
        "longitude",
        "x",
        "X",
        "LON",
        "XLON",
        "LONS",
        "LONG",
        "LONGITUDE",
        "mlon",
    }
    lat_names = {
        "lat",
        "latitude",
        "y",
        "Y",
        "LAT",
        "YLAT",
        "LATS",
        "LATI",
        "LATITUDE",
        "nlat",
    }
    lev_names = {
        "level",
        "lev",
        "plev",
        "plevel",
        "p",
        "plevs",
        "pressure",
        "pressure_level",
        "z",
        "Z",
        "LEV",
        "PRES",
        "LEVEL",
        "PLEVEL",
        "PRES_LEVEL",
        "HEIGHT",
        "height",
        "altitude",
        "ALTITUDE",
        "depth",
        "DEPTH",
        "isobaric",
        "ISOBARIC",
        "model_level",
        "PRESSURE_LEVEL",
        "PRESSURE",
    }
    time_names = {"time", "t", "T", "year", "month", "yr", "mn", "season"}

    dims_found = {"xdim": None, "ydim": None, "levdim": None, "timedim": None}
    dim_types = [
        (lon_names, "xdim"),
        (lat_names, "ydim"),
        (lev_names, "levdim"),
        (time_names, "timedim"),
    ]

    for i, dim_name in enumerate(dims):
        dim_lower = dim_name.lower()
        for names, var_name in dim_types:
            if any(name in dim_lower for name in names):
                dims_found[var_name] = i
                break

    xdim, ydim, levdim, timedim = (
        dims_found["xdim"],
        dims_found["ydim"],
        dims_found["levdim"],
        dims_found["timedim"],
    )

    # Level dimension is required
    if levdim is None:
        raise ValueError(
            f"Could not auto-detect level dimension. Found dims: {dims}. "
            f"Please specify vertical_dim explicitly."
        )

    # Check dimension requirements based on data dimensionality
    ndim = data_array.ndim
    if ndim == 1:
        # 1D profile data - only need level dimension
        pass
    elif ndim == 3:
        # 3D data - need both lat and lon dimensions
        if xdim is None or ydim is None:
            raise ValueError(
                f"For 3D data, need both lat and lon dimensions. Found dims: {dims}. "
                f"Please specify spatial dimensions explicitly."
            )
    elif ndim == 4:
        # 4D data - need time, level, lat, and lon dimensions
        if xdim is None or ydim is None or timedim is None:
            raise ValueError(
                f"For 4D data, need time, lat, and lon dimensions. Found dims: {dims}. "
                f"Please specify dimensions explicitly."
            )

    return xdim, ydim, levdim, timedim


def _check_units(
    data_array: xr.DataArray, variable_name: str, expected_units: list
) -> tuple[xr.DataArray, bool]:
    """
    Check and validate units for atmospheric variables with comprehensive unit recognition.

    Parameters
    ----------
    data_array : xr.DataArray
        Input data array to check
    variable_name : str
        Name of the variable for error messages
    expected_units : list
        List of acceptable units

    Returns
    -------
    data_array : xr.DataArray
        Potentially converted data array
    converted : bool
        Whether unit conversion was applied
    """
    if not hasattr(data_array, "attrs") or "units" not in data_array.attrs:
        warnings.warn(
            f"{variable_name} has no 'units' attribute. Assuming correct units: {expected_units[0]}",
            UserWarning,
        )
        return data_array, False

    # Get units and clean up common formatting issues
    units = data_array.attrs["units"].strip()
    units = units.replace("**", "^")  # Handle ** for exponents
    units_lower = units.lower()

    if units in expected_units:
        return data_array, False

    # Handle common unit conversions
    converted_data = data_array.copy()

    # Temperature conversions
    temp_var_names = [
        "temperature",
        "temp",
        "t",
        "air",
        "sst",
        "sea_surface_temperature",
        "tair",
        "air_temperature",
        "t2m",
        "tas",
        "ta",
    ]
    if variable_name.lower() in temp_var_names and expected_units[0] == "K":
        # Celsius variations
        celsius_units = [
            "C",
            "°C",
            "celsius",
            "Celsius",
            "degrees C",
            "degC",
            "deg C",
            "degrees_C",
            "degrees_celsius",
            "degree_C",
            "degree_Celsius",
            "Degrees_Celsius",
            "CELSIUS",
            "celcius",
            "Celcius",  # common misspellings
            "℃",
            "ºC",
            "oC",
            "deg_C",
            "degree celsius",
            "degrees celsius",
        ]
        if units in celsius_units or units_lower in [u.lower() for u in celsius_units]:
            converted_data = data_array + 273.15
            converted_data.attrs["units"] = "K"
            warnings.warn(f"Converted {variable_name} from {units} to K", UserWarning)
            return converted_data, True

        # Fahrenheit variations
        fahrenheit_units = [
            "F",
            "°F",
            "fahrenheit",
            "Fahrenheit",
            "degrees F",
            "degF",
            "deg F",
            "degrees_F",
            "degrees_fahrenheit",
            "degree_F",
            "degree_Fahrenheit",
            "FAHRENHEIT",
            "℉",
            "ºF",
            "oF",
            "deg_F",
        ]
        if units in fahrenheit_units or units_lower in [
            u.lower() for u in fahrenheit_units
        ]:
            converted_data = (data_array - 32) * 5 / 9 + 273.15
            converted_data.attrs["units"] = "K"
            warnings.warn(f"Converted {variable_name} from {units} to K", UserWarning)
            return converted_data, True

        # Kelvin variations (already correct but may need standardization)
        kelvin_units = [
            "K",
            "kelvin",
            "Kelvin",
            "KELVIN",
            "degrees_kelvin",
            "degrees K",
            "degK",
            "deg K",
            "degree_K",
            "degree_Kelvin",
        ]
        if units in kelvin_units or units_lower in [u.lower() for u in kelvin_units]:
            if units != "K":
                converted_data.attrs["units"] = "K"
                warnings.warn(
                    f"Standardized {variable_name} units from {units} to K", UserWarning
                )
            return converted_data, False

    # Pressure conversions
    pressure_var_names = [
        "pressure",
        "pres",
        "p",
        "psl",
        "msl",
        "pressure_levels",
        "plev",
        "level",
        "slp",
        "sea_level_pressure",
        "pmsl",
        "air_pressure",
        "atmospheric_pressure",
    ]
    if variable_name.lower() in pressure_var_names and expected_units[0] in [
        "Pa",
        "hPa",
        "mb",
    ]:
        # Pascal variations
        pascal_units = [
            "Pa",
            "pa",
            "pascal",
            "Pascal",
            "PASCAL",
            "Pascals",
            "pascals",
            "N/m^2",
            "N/m2",
            "N m^-2",
            "N m-2",
            "N.m-2",
            "N·m^-2",
            "newton/m^2",
        ]
        if units in pascal_units or units_lower in [u.lower() for u in pascal_units]:
            if expected_units[0] in ["hPa", "mb"]:
                converted_data = data_array / 100.0
                converted_data.attrs["units"] = "hPa"
                warnings.warn(
                    f"Converted {variable_name} from {units} to hPa", UserWarning
                )
                return converted_data, True
            elif units != "Pa":
                converted_data.attrs["units"] = "Pa"
                warnings.warn(
                    f"Standardized {variable_name} units from {units} to Pa",
                    UserWarning,
                )
                return converted_data, False

        # hPa/mb variations
        hpa_mb_units = [
            "hPa",
            "mb",
            "mbar",
            "millibar",
            "millibars",
            "hpa",
            "hectopascal",
            "hectopascals",
            "hecto-pascal",
            "hecto pascal",
            "hPascal",
            "hpascal",
            "MILLIBAR",
            "MBAR",
            "HPA",
            "HECTOPASCAL",
            "mbars",
        ]
        if units in hpa_mb_units or units_lower in [u.lower() for u in hpa_mb_units]:
            if expected_units[0] == "Pa":
                converted_data = data_array * 100.0
                converted_data.attrs["units"] = "Pa"
                warnings.warn(
                    f"Converted {variable_name} from {units} to Pa", UserWarning
                )
                return converted_data, True
            elif units not in ["hPa", "mb"]:
                converted_data.attrs["units"] = expected_units[0]
                warnings.warn(
                    f"Standardized {variable_name} units from {units} to {expected_units[0]}",
                    UserWarning,
                )
                return converted_data, False

        # kPa to Pa/hPa conversion
        kpa_units = [
            "kPa",
            "kpa",
            "kilopascal",
            "kilopascals",
            "kilo-pascal",
            "kilo pascal",
            "kPascal",
            "kpascal",
            "KILOPASCAL",
            "KPA",
        ]
        if units in kpa_units or units_lower in [u.lower() for u in kpa_units]:
            if expected_units[0] == "Pa":
                converted_data = data_array * 1000.0
                converted_data.attrs["units"] = "Pa"
            else:  # hPa or mb
                converted_data = data_array * 10.0
                converted_data.attrs["units"] = "hPa"
            warnings.warn(
                f"Converted {variable_name} from {units} to {converted_data.attrs['units']}",
                UserWarning,
            )
            return converted_data, True

        # atm to Pa/hPa conversion
        atm_units = [
            "atm",
            "atmosphere",
            "atmospheres",
            "standard atmosphere",
            "ATM",
            "ATMOSPHERE",
        ]
        if units in atm_units or units_lower in [u.lower() for u in atm_units]:
            if expected_units[0] == "Pa":
                converted_data = data_array * 101325.0
                converted_data.attrs["units"] = "Pa"
            else:  # hPa or mb
                converted_data = data_array * 1013.25
                converted_data.attrs["units"] = "hPa"
            warnings.warn(
                f"Converted {variable_name} from {units} to {converted_data.attrs['units']}",
                UserWarning,
            )
            return converted_data, True

    # Mixing ratio/humidity conversions
    humidity_var_names = [
        "mixing_ratio",
        "specific_humidity",
        "q",
        "qv",
        "hus",
        "hur",
        "water_vapor",
        "water_vapor_mixing_ratio",
        "humidity",
        "shum",
        "spec_hum",
        "spec_humidity",
        "qvapor",
        "qvap",
    ]
    if variable_name.lower() in humidity_var_names and expected_units[0] == "kg/kg":
        # g/kg variations
        g_per_kg_units = [
            "g/kg",
            "g kg-1",
            "g kg^-1",
            "g.kg-1",
            "g.kg^-1",
            "g·kg^-1",
            "grams/kg",
            "grams per kilogram",
            "g/kilogram",
            "gram/kg",
            "g kg^(-1)",
            "g*kg^-1",
            "G/KG",
            "g.kg^(-1)",
            "g kg^(-1)",
        ]
        if units in g_per_kg_units or units_lower in [
            u.lower() for u in g_per_kg_units
        ]:
            converted_data = data_array / 1000.0  # g/kg to kg/kg
            converted_data.attrs["units"] = "kg/kg"
            warnings.warn(
                f"Converted {variable_name} from {units} to kg/kg", UserWarning
            )
            return converted_data, True

        # mg/kg variations (less common but possible)
        mg_per_kg_units = ["mg/kg", "mg kg-1", "mg kg^-1", "mg.kg-1", "milligrams/kg"]
        if units in mg_per_kg_units or units_lower in [
            u.lower() for u in mg_per_kg_units
        ]:
            converted_data = data_array / 1000000.0  # mg/kg to kg/kg
            converted_data.attrs["units"] = "kg/kg"
            warnings.warn(
                f"Converted {variable_name} from {units} to kg/kg", UserWarning
            )
            return converted_data, True

        # kg/kg variations (already correct but may need standardization)
        kg_per_kg_units = [
            "kg/kg",
            "kg kg-1",
            "kg kg^-1",
            "kg.kg-1",
            "kg.kg^-1",
            "kg·kg^-1",
            "kilograms/kilogram",
            "kg*kg^-1",
            "KG/KG",
            "kg.kg^(-1)",
            "kg kg^(-1)",
        ]
        dimensionless_units = [
            "dimensionless",
            "1",
            "none",
            "unitless",
            "-",
            "",
            "fraction",
            "ratio",
            "nondimensional",
            "no unit",
            "DIMENSIONLESS",
        ]
        all_correct_units = kg_per_kg_units + dimensionless_units
        if units in all_correct_units or units_lower in [
            u.lower() for u in all_correct_units
        ]:
            if units != "kg/kg":
                converted_data.attrs["units"] = "kg/kg"
                warnings.warn(
                    f"Standardized {variable_name} units from {units} to kg/kg",
                    UserWarning,
                )
            return converted_data, False

        # Percentage to kg/kg (rare but possible for humidity)
        percentage_units = ["%", "percent", "percentage", "pct", "PERCENT", "％"]
        if units in percentage_units or units_lower in [
            u.lower() for u in percentage_units
        ]:
            converted_data = data_array / 100.0  # percent to fraction
            converted_data.attrs["units"] = "kg/kg"
            warnings.warn(
                f"Converted {variable_name} from {units} to kg/kg (assuming fraction representation)",
                UserWarning,
            )
            return converted_data, True

        # ppmv to kg/kg (parts per million by volume - used in some datasets)
        ppmv_units = ["ppmv", "ppm", "parts per million", "PPM", "PPMV", "ppm(v)"]
        if units in ppmv_units or units_lower in [u.lower() for u in ppmv_units]:
            # ppmv to mass mixing ratio approximation for water vapor
            # This is a rough conversion; actual conversion depends on molecular weights
            converted_data = data_array * 0.622e-6  # approximate for water vapor
            converted_data.attrs["units"] = "kg/kg"
            warnings.warn(
                f"Converted {variable_name} from {units} to kg/kg (approximate conversion)",
                UserWarning,
            )
            return converted_data, True

    # If no conversion available, raise error
    raise ValueError(
        f"{variable_name} has units '{units}' but expected one of {expected_units}. "
        f"Please convert to appropriate units before calling this function."
    )


def _create_output_dataset(
    min_pressure,
    pi,
    error_flag,
    input_coords=None,
    sst_val=None,
    psl_val=None,
    vertical_levels=None,
    vertical_dim="level",
    data_type="profile",
):
    """
    Helper function to create output dataset with proper metadata.

    Parameters
    ----------
    min_pressure, pi, error_flag : array-like
        Calculation results
    input_coords : dict, optional
        Coordinate dictionary for spatial/temporal dimensions
    sst_val, psl_val : float, optional
        Input SST and PSL values for metadata
    vertical_levels : int, optional
        Number of vertical levels for metadata
    vertical_dim : str
        Name of vertical dimension
    data_type : str
        Type of calculation (profile, 3D, 4D)
    """
    # Determine output dimensions based on input data
    dims = list(input_coords.keys()) if input_coords else []
    coords = input_coords or {}

    # Create data variables
    data_vars = {
        "min_pressure": (
            dims,
            min_pressure,
            {
                "long_name": "Minimum central pressure",
                "units": "mb",
                "description": "Tropical cyclone minimum central pressure",
            },
        ),
        "pi": (
            dims,
            pi,
            {
                "long_name": "Potential intensity",
                "units": "m/s",
                "description": "Tropical cyclone maximum potential intensity (maximum sustained wind speed)",
            },
        ),
        "error_flag": (
            [],  # error_flag is always a scalar, regardless of input dimensions
            int(error_flag) if np.isscalar(error_flag) else int(error_flag),
            {
                "long_name": "Error flag",
                "units": "dimensionless",
                "description": "Error status (1 = success, other values = error)",
            },
        ),
    }

    # Create global attributes
    attrs = {
        "title": "Tropical Cyclone Potential Intensity Analysis",
        "description": f"{data_type} potential intensity calculation",
        "vertical_dim": vertical_dim,
        "method": "Emanuel potential intensity calculation",
        "reference": "Emanuel (1995), Bister & Emanuel (2002)",
    }

    if sst_val is not None:
        attrs.update(
            {
                "sst_input": float(sst_val) if np.isscalar(sst_val) else "variable",
                "sst_units": "K",
            }
        )
    if psl_val is not None:
        attrs.update(
            {
                "psl_input": float(psl_val) if np.isscalar(psl_val) else "variable",
                "psl_units": "Pa",
            }
        )
    if vertical_levels is not None:
        attrs["vertical_levels"] = vertical_levels

    return xr.Dataset(data_vars=data_vars, coords=coords, attrs=attrs)


def _potential_intensity_profile(
    sst: Union[float, xr.DataArray],
    psl: Union[float, xr.DataArray],
    pressure_levels: xr.DataArray,
    temperature: xr.DataArray,
    mixing_ratio: xr.DataArray,
) -> xr.Dataset:
    """
    Calculate tropical cyclone potential intensity for single atmospheric profile using xarray.

    This function provides an xarray-native interface for single column potential intensity
    calculations, handling coordinate preservation and metadata.

    Parameters
    ----------
    sst : float or xr.DataArray
        Sea surface temperature [K]. If DataArray, should be 0-dimensional (scalar).
    psl : float or xr.DataArray
        Sea level pressure [Pa]. If DataArray, should be 0-dimensional (scalar).
    pressure_levels : xr.DataArray
        Atmospheric pressure levels [mb]. Should be 1D with dimension specified by vertical_dim.
    temperature : xr.DataArray
        Temperature profile [K]. Should be 1D with vertical dimension.
    mixing_ratio : xr.DataArray
        Water vapor mixing ratio profile [kg/kg]. Should be 1D with same dimension as temperature.

    Returns
    -------
    result : xr.Dataset
        Dataset containing potential intensity results as scalar variables.
    """

    # Auto-detect vertical dimension
    _, _, levdim, _ = _detect_atmospheric_dimensions(temperature)
    vertical_dim = list(temperature.dims)[levdim]

    # Sort data by pressure levels in descending order (from surface to top)
    pressure_levels = pressure_levels.sortby(vertical_dim, ascending=False)
    temperature = temperature.sortby(vertical_dim, ascending=False)
    mixing_ratio = mixing_ratio.sortby(vertical_dim, ascending=False)

    # Ensure 1D profiles
    if temperature.ndim != 1 or mixing_ratio.ndim != 1 or pressure_levels.ndim != 1:
        raise ValueError(
            "Only 1D profiles are supported. All arrays should have only the vertical dimension."
        )

    # Extract scalar values for SST and PSL
    def _extract_scalar(value, name):
        if isinstance(value, xr.DataArray):
            if value.ndim != 0:
                raise ValueError(f"{name} DataArray must be 0-dimensional (scalar)")
            return float(value.values)
        return float(value)

    sst_val = _extract_scalar(sst, "SST")
    psl_val = _extract_scalar(psl, "PSL")

    # Extract numpy arrays for calculation
    p_levels = pressure_levels.values
    temp_vals = temperature.values
    mixr_vals = mixing_ratio.values

    # Perform calculation using existing interface
    min_pressure, pi, error_flag = calculate_potential_intensity_profile(
        sst_val, psl_val, p_levels, temp_vals, mixr_vals
    )

    # Create output dataset
    return _create_output_dataset(
        min_pressure,
        pi,
        error_flag,
        sst_val=sst_val,
        psl_val=psl_val,
        vertical_levels=len(p_levels),
        vertical_dim=vertical_dim,
        data_type="Single column profile",
    )


def _potential_intensity_3d(
    sst: xr.DataArray,
    psl: xr.DataArray,
    pressure_levels: xr.DataArray,
    temperature: xr.DataArray,
    mixing_ratio: xr.DataArray,
) -> xr.Dataset:
    """
    Calculate tropical cyclone potential intensity for 3D gridded data using xarray.

    Parameters
    ----------
    sst : xr.DataArray, shape (nlat, nlon) or (y, x)
        Sea surface temperature [K]
    psl : xr.DataArray, shape (nlat, nlon) or (y, x)
        Sea level pressure [Pa]
    pressure_levels : xr.DataArray, shape (num_levels,)
        Atmospheric pressure levels [mb]
    temperature : xr.DataArray, shape (num_levels, nlat, nlon) or (num_levels, y, x)
        Temperature profiles [K]
    mixing_ratio : xr.DataArray, shape (num_levels, nlat, nlon) or (num_levels, y, x)
        Water vapor mixing ratio profiles [kg/kg]

    Returns
    -------
    result : xr.Dataset
        Dataset containing potential intensity results with spatial dimensions.
    """

    # Auto-detect dimensions
    _, _, levdim, _ = _detect_atmospheric_dimensions(temperature)
    vertical_dim = list(temperature.dims)[levdim]

    # Sort data by pressure levels in descending order (from surface to top)
    pressure_levels = pressure_levels.sortby(vertical_dim, ascending=False)
    temperature = temperature.sortby(vertical_dim, ascending=False)
    mixing_ratio = mixing_ratio.sortby(vertical_dim, ascending=False)

    # Ensure temperature and mixing_ratio are 3D (level + 2 spatial dims)
    if temperature.ndim != 3 or mixing_ratio.ndim != 3:
        raise ValueError(
            "temperature and mixing_ratio must be 3D arrays (vertical + 2 spatial dimensions)"
        )

    # Ensure SST and PSL are 2D
    if sst.ndim != 2 or psl.ndim != 2:
        raise ValueError("sst and psl must be 2D arrays (2 spatial dimensions)")

    # Get spatial dimensions (all dims except vertical_dim from temperature)
    spatial_dims = [d for d in temperature.dims if d != vertical_dim]
    if len(spatial_dims) != 2:
        raise ValueError("Expected exactly 2 spatial dimensions")

    # Reorder dimensions to match expected Fortran order: (vertical, spatial1, spatial2)
    temp_reordered = temperature.transpose(vertical_dim, *spatial_dims)
    mixr_reordered = mixing_ratio.transpose(vertical_dim, *spatial_dims)

    # Extract numpy arrays
    sst_vals = sst.values
    psl_vals = psl.values
    p_levels = pressure_levels.values
    temp_vals = temp_reordered.values
    mixr_vals = mixr_reordered.values

    # Perform calculation
    min_pressure, pi, error_flag = calculate_potential_intensity_3d(
        sst_vals, psl_vals, p_levels, temp_vals, mixr_vals
    )

    # Create output coordinates (spatial dimensions only)
    output_coords = {dim: temperature.coords[dim] for dim in spatial_dims}

    return _create_output_dataset(
        min_pressure,
        pi,
        error_flag,
        input_coords=output_coords,
        vertical_levels=len(p_levels),
        vertical_dim=vertical_dim,
        data_type="3D gridded",
    )


def _potential_intensity_4d(
    sst: xr.DataArray,
    psl: xr.DataArray,
    pressure_levels: xr.DataArray,
    temperature: xr.DataArray,
    mixing_ratio: xr.DataArray,
) -> xr.Dataset:
    """
    Calculate tropical cyclone potential intensity for 4D time series data using xarray.

    Parameters
    ----------
    sst : xr.DataArray, shape (ntimes, nlat, nlon) or (time, y, x)
        Sea surface temperature [K]
    psl : xr.DataArray, shape (ntimes, nlat, nlon) or (time, y, x)
        Sea level pressure [Pa]
    pressure_levels : xr.DataArray, shape (num_levels,)
        Atmospheric pressure levels [mb]
    temperature : xr.DataArray, shape (ntimes, num_levels, nlat, nlon) or (time, level, y, x)
        Temperature profiles [K]
    mixing_ratio : xr.DataArray, shape (ntimes, num_levels, nlat, nlon) or (time, level, y, x)
        Water vapor mixing ratio profiles [kg/kg]

    Returns
    -------
    result : xr.Dataset
        Dataset containing potential intensity results with time and spatial dimensions.
    """

    # Auto-detect dimensions
    _, _, levdim, timedim = _detect_atmospheric_dimensions(temperature)
    vertical_dim = list(temperature.dims)[levdim]
    time_dim = list(temperature.dims)[timedim]

    # Sort data by pressure levels in descending order (from surface to top)
    pressure_levels = pressure_levels.sortby(vertical_dim, ascending=False)
    temperature = temperature.sortby(vertical_dim, ascending=False)
    mixing_ratio = mixing_ratio.sortby(vertical_dim, ascending=False)

    # Ensure correct number of dimensions
    if temperature.ndim != 4 or mixing_ratio.ndim != 4:
        raise ValueError(
            "temperature and mixing_ratio must be 4D arrays (time + vertical + 2 spatial dimensions)"
        )

    if sst.ndim != 3 or psl.ndim != 3:
        raise ValueError("sst and psl must be 3D arrays (time + 2 spatial dimensions)")

    # Get output dimensions (all dims except vertical_dim from temperature)
    output_dims = [d for d in temperature.dims if d != vertical_dim]
    if len(output_dims) != 3:  # time + 2 spatial
        raise ValueError("Expected time dimension plus exactly 2 spatial dimensions")

    # Reorder dimensions to match expected Fortran order: (time, vertical, spatial1, spatial2)
    expected_dims = [time_dim, vertical_dim] + [d for d in output_dims if d != time_dim]
    temp_reordered = temperature.transpose(*expected_dims)
    mixr_reordered = mixing_ratio.transpose(*expected_dims)

    # Extract numpy arrays
    sst_vals = sst.values
    psl_vals = psl.values
    p_levels = pressure_levels.values
    temp_vals = temp_reordered.values
    mixr_vals = mixr_reordered.values

    # Perform calculation
    min_pressure, pi, error_flag = calculate_potential_intensity_4d(
        sst_vals, psl_vals, p_levels, temp_vals, mixr_vals
    )

    # Create output coordinates (time + spatial dimensions)
    output_coords = {dim: temperature.coords[dim] for dim in output_dims}

    return _create_output_dataset(
        min_pressure,
        pi,
        error_flag,
        input_coords=output_coords,
        vertical_levels=len(p_levels),
        vertical_dim=vertical_dim,
        data_type="4D time series",
    )


def potential_intensity(
    sst: xr.DataArray,
    psl: xr.DataArray,
    pressure_levels: xr.DataArray,
    temperature: xr.DataArray,
    mixing_ratio: xr.DataArray,
) -> xr.Dataset:
    """
    Calculate tropical cyclone potential intensity with automatic dimension detection.

    This function automatically detects the input data dimensions and calls the
    appropriate specific function (profile, 3D, or 4D).

    Parameters
    ----------
    sst : xr.DataArray
        Sea surface temperature [K]
    psl : xr.DataArray
        Sea level pressure [Pa]
    pressure_levels : xr.DataArray
        Atmospheric pressure levels [mb]
    temperature : xr.DataArray
        Temperature data [K]
    mixing_ratio : xr.DataArray
        Water vapor mixing ratio [kg/kg]

    Returns
    -------
    result : xr.Dataset
        Dataset containing potential intensity results with appropriate dimensions.

    Notes
    -----
    The function determines the calculation type based on temperature array dimensions:
    - 1D (vertical only): Single profile calculation
    - 3D (vertical + 2 spatial): Gridded data calculation
    - 4D (time + vertical + 2 spatial): Time series calculation
    """

    # Check and convert units if needed
    sst, _ = _check_units(sst, "SST", ["K"])
    psl, _ = _check_units(psl, "PSL", ["Pa"])
    pressure_levels, _ = _check_units(pressure_levels, "pressure_levels", ["mb", "hPa"])
    temperature, _ = _check_units(temperature, "temperature", ["K"])
    mixing_ratio, _ = _check_units(mixing_ratio, "mixing_ratio", ["kg/kg"])
    # Auto-detect calculation type based on temperature dimensions
    ndim = temperature.ndim

    # Function mapping for different dimensions
    func_map = {
        1: _potential_intensity_profile,  # Single profile
        3: _potential_intensity_3d,  # 3D gridded
        4: _potential_intensity_4d,  # 4D time series
    }

    if ndim not in func_map:
        raise ValueError(
            f"Unsupported number of dimensions: {ndim}. Expected 1, 3, or 4."
        )

    return func_map[ndim](sst, psl, pressure_levels, temperature, mixing_ratio)
