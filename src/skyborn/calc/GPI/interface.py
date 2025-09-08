"""
High-level Python interface for Tropical Cyclone Genesis Potential Index (GPI)
and Potential Intensity (PI) calculations.

This module provides user-friendly interfaces for calculating tropical cyclone
potential intensity from atmospheric and oceanic data, with support for
multi-dimensional data arrays and proper handling of missing values.

The interface handles automatic data validation, unit conversions, and
integration with the optimized Fortran backend.
"""

import warnings
from typing import Any, Dict, Optional, Tuple, Union

import numpy as np

from . import tropical_cyclone_potential_intensity as _gpi_module

# Fortran UNDEF value constant
UNDEF = -9.99e33


def _postprocess_results(min_pressure, max_wind):
    """Convert UNDEF values to NaN in results."""
    min_pressure = np.where(min_pressure == UNDEF, np.nan, min_pressure)
    max_wind = np.where(max_wind == UNDEF, np.nan, max_wind)
    return min_pressure, max_wind


def _validate_input_arrays(*arrays, names=None):
    """Validate input arrays for NaN/inf values, convert to float32, and detect missing values."""
    if names is None:
        names = [f"array_{i}" for i in range(len(arrays))]

    validated = []
    has_missing = False
    missing_sources = []

    for arr, name in zip(arrays, names):
        arr = np.asarray(arr, dtype=np.float32)

        # Check for various types of missing values
        has_undef = np.any(arr == UNDEF)
        has_nan = np.any(np.isnan(arr))
        has_inf = np.any(np.isinf(arr))

        if has_undef or has_nan:
            has_missing = True
            if has_undef:
                missing_sources.append(f"{name}(UNDEF)")
            if has_nan:
                missing_sources.append(f"{name}(NaN)")

        if has_inf:
            warnings.warn(
                f"Infinite values detected in {name}. Results may be unreliable."
            )

        # Convert NaN to UNDEF for consistent Fortran handling
        if has_nan:
            arr = np.where(np.isnan(arr), UNDEF, arr)

        validated.append(arr)

    if has_missing:
        print(
            f"Missing values detected in: {', '.join(missing_sources)}. Using missing value handling version."
        )

    result = validated if len(validated) > 1 else validated[0]
    return result, has_missing


def _ensure_pressure_ordering(pressure_levels, temperature, mixing_ratio):
    """
    Ensure pressure levels are ordered from surface to top (high to low pressure).

    The Fortran code expects pressure levels with:
    - Index 1 = highest pressure (surface/ground level)
    - Index N = lowest pressure (top of atmosphere)

    Parameters
    ----------
    pressure_levels : ndarray
        Pressure levels [mb]
    temperature : ndarray
        Temperature data with pressure as first or second dimension
    mixing_ratio : ndarray
        Mixing ratio data with same shape as temperature

    Returns
    -------
    pressure_levels_ordered : ndarray
        Pressure levels ordered surface to top (high to low)
    temperature_ordered : ndarray
        Temperature data reordered to match pressure ordering
    mixing_ratio_ordered : ndarray
        Mixing ratio data reordered to match pressure ordering
    """
    pressure_levels = np.asarray(pressure_levels)

    # Check if pressure levels need reordering (only reverse if necessary)
    if len(pressure_levels) > 1 and pressure_levels[0] < pressure_levels[-1]:
        warnings.warn(
            "Pressure levels appear to be ordered from top to surface (low to high pressure). "
            "Reordering to surface to top (high to low pressure) as required by the calculation."
        )

        # Reverse arrays - views are memory-efficient, only copy when necessary
        pressure_axis = 0 if temperature.ndim in [1, 3] else 1
        return (
            pressure_levels[::-1],
            np.flip(temperature, axis=pressure_axis),
            np.flip(mixing_ratio, axis=pressure_axis),
        )
    else:
        # Already correctly ordered or single level - return original arrays
        return pressure_levels, temperature, mixing_ratio


def _validate_dimensions(sst, psl, pressure_levels, temp, mixing_ratio, data_type="3D"):
    """Validate that input arrays have compatible dimensions based on temperature array."""

    # Temperature and mixing_ratio must have same shape
    if temp.shape != mixing_ratio.shape:
        raise ValueError(
            f"Temperature shape {temp.shape} doesn't match mixing ratio shape {mixing_ratio.shape}"
        )

    expected_levels = len(pressure_levels)

    if data_type == "profile":
        # Profile: temp.shape = (num_levels,), SST/PSL scalars
        if temp.shape != (expected_levels,):
            raise ValueError(
                f"Temperature shape {temp.shape} doesn't match expected profile shape ({expected_levels},)"
            )
        if not (np.isscalar(sst) or sst.ndim == 0) or not (
            np.isscalar(psl) or psl.ndim == 0
        ):
            raise ValueError("SST and PSL must be scalars for profile data")

    elif data_type == "3D":
        # 3D: temp.shape = (num_levels, nlat, nlon), SST/PSL.shape = (nlat, nlon)
        if temp.ndim != 3 or temp.shape[0] != expected_levels:
            raise ValueError(
                f"Temperature shape {temp.shape} doesn't match expected 3D shape ({expected_levels}, nlat, nlon)"
            )
        expected_sst_shape = temp.shape[1:]  # (nlat, nlon)
        if sst.shape != expected_sst_shape or psl.shape != expected_sst_shape:
            raise ValueError(
                f"SST/PSL shape mismatch - expected {expected_sst_shape}, got SST:{sst.shape}, PSL:{psl.shape}"
            )

    elif data_type == "4D":
        # 4D: temp.shape = (ntimes, num_levels, nlat, nlon), SST/PSL.shape = (ntimes, nlat, nlon)
        if temp.ndim != 4 or temp.shape[1] != expected_levels:
            raise ValueError(
                f"Temperature shape {temp.shape} doesn't match expected 4D shape (ntimes, {expected_levels}, nlat, nlon)"
            )
        expected_sst_shape = (
            temp.shape[0],
            temp.shape[2],
            temp.shape[3],
        )  # (ntimes, nlat, nlon)
        if sst.shape != expected_sst_shape or psl.shape != expected_sst_shape:
            raise ValueError(
                f"SST/PSL shape mismatch - expected {expected_sst_shape}, got SST:{sst.shape}, PSL:{psl.shape}"
            )

    else:
        raise ValueError(f"Unsupported data type: {data_type}")


def calculate_potential_intensity_3d(
    sst: np.ndarray,
    psl: np.ndarray,
    pressure_levels: np.ndarray,
    temperature: np.ndarray,
    mixing_ratio: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Calculate tropical cyclone potential intensity for 3D gridded data (spatial grid with vertical levels).

    Parameters
    ----------
    sst : ndarray, shape (nlat, nlon)
        Sea surface temperature [K]
    psl : ndarray, shape (nlat, nlon)
        Sea level pressure [Pa]
    pressure_levels : ndarray, shape (num_levels,)
        Atmospheric pressure levels [mb]. Can be in any order - will be automatically
        reordered from surface to top (high to low pressure) as required by the calculation
    temperature : ndarray, shape (num_levels, nlat, nlon)
        Temperature profiles [K]
    mixing_ratio : ndarray, shape (num_levels, nlat, nlon)
        Water vapor mixing ratio [kg/kg]

    Returns
    -------
    min_pressure : ndarray, shape (nlat, nlon)
        Minimum central pressure [mb]
    max_wind : ndarray, shape (nlat, nlon)
        Maximum sustained wind speed [m/s]
    error_flag : int
        Error status (0 = success, non-zero = error)

    Examples
    --------
    >>> # 3D calculation with automatic missing value detection
    >>> min_p, max_w, err = calculate_potential_intensity_3d(
    ...     sst, psl, p_levels, temp, mixr)
    """
    # Validate and convert inputs, detect missing values
    (sst, psl, pressure_levels, temperature, mixing_ratio), has_missing = (
        _validate_input_arrays(
            sst,
            psl,
            pressure_levels,
            temperature,
            mixing_ratio,
            names=["SST", "PSL", "pressure_levels", "temperature", "mixing_ratio"],
        )
    )

    # Ensure correct pressure level ordering (surface to top)
    pressure_levels, temperature, mixing_ratio = _ensure_pressure_ordering(
        pressure_levels, temperature, mixing_ratio
    )

    # Validate dimensions
    _validate_dimensions(sst, psl, pressure_levels, temperature, mixing_ratio, "3D")

    # Choose appropriate Fortran function based on missing value detection
    func = (
        _gpi_module.calculate_pi_gridded_with_missing
        if has_missing
        else _gpi_module.calculate_pi_gridded_data
    )

    # Call Fortran function (dimensions are automatically inferred)
    min_pressure, max_wind, error_flag = func(
        sst, psl, pressure_levels, temperature, mixing_ratio
    )

    # Convert UNDEF to NaN in results
    min_pressure, max_wind = _postprocess_results(min_pressure, max_wind)

    return min_pressure, max_wind, error_flag


def calculate_potential_intensity_4d(
    sst: np.ndarray,
    psl: np.ndarray,
    pressure_levels: np.ndarray,
    temperature: np.ndarray,
    mixing_ratio: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Calculate tropical cyclone potential intensity for 4D time series data.

    Parameters
    ----------
    sst : ndarray, shape (ntimes, nlat, nlon)
        Sea surface temperature [K]
    psl : ndarray, shape (ntimes, nlat, nlon)
        Sea level pressure [Pa]
    pressure_levels : ndarray, shape (num_levels,)
        Atmospheric pressure levels [mb]. Can be in any order - will be automatically
        reordered from surface to top (high to low pressure) as required by the calculation
    temperature : ndarray, shape (ntimes, num_levels, nlat, nlon)
        Temperature profiles [K]
    mixing_ratio : ndarray, shape (ntimes, num_levels, nlat, nlon)
        Water vapor mixing ratio [kg/kg]

    Returns
    -------
    min_pressure : ndarray, shape (ntimes, nlat, nlon)
        Minimum central pressure [mb]
    max_wind : ndarray, shape (ntimes, nlat, nlon)
        Maximum sustained wind speed [m/s]
    error_flag : int
        Error status (0 = success, non-zero = error)

    Examples
    --------
    >>> # 4D time series calculation
    >>> min_p, max_w, err = calculate_potential_intensity_4d(
    ...     sst_4d, psl_4d, p_levels, temp_4d, mixr_4d)
    """
    # Validate and convert inputs, detect missing values
    (sst, psl, pressure_levels, temperature, mixing_ratio), has_missing = (
        _validate_input_arrays(
            sst,
            psl,
            pressure_levels,
            temperature,
            mixing_ratio,
            names=["SST", "PSL", "pressure_levels", "temperature", "mixing_ratio"],
        )
    )

    # Ensure correct pressure level ordering (surface to top)
    pressure_levels, temperature, mixing_ratio = _ensure_pressure_ordering(
        pressure_levels, temperature, mixing_ratio
    )

    # Validate dimensions
    _validate_dimensions(sst, psl, pressure_levels, temperature, mixing_ratio, "4D")

    # Choose appropriate Fortran function based on missing value detection
    func = (
        _gpi_module.calculate_pi_4d_with_missing
        if has_missing
        else _gpi_module.calculate_pi_4d_data
    )

    # Call Fortran function (dimensions are automatically inferred)
    min_pressure, max_wind, error_flag = func(
        sst, psl, pressure_levels, temperature, mixing_ratio
    )

    # Convert UNDEF to NaN in results
    min_pressure, max_wind = _postprocess_results(min_pressure, max_wind)

    return min_pressure, max_wind, error_flag


def calculate_potential_intensity_profile(
    sst: float,
    psl: float,
    pressure_levels: np.ndarray,
    temperature: np.ndarray,
    mixing_ratio: np.ndarray,
    actual_levels: Optional[int] = None,
) -> Tuple[float, float, int]:
    """
    Calculate tropical cyclone potential intensity for a single atmospheric profile.

    Parameters
    ----------
    sst : float
        Sea surface temperature [K]
    psl : float
        Sea level pressure [Pa]
    pressure_levels : ndarray, shape (num_levels,)
        Atmospheric pressure levels [mb]. Can be in any order - will be automatically
        reordered from surface to top (high to low pressure) as required by the calculation
    temperature : ndarray, shape (num_levels,)
        Temperature profile [K]
    mixing_ratio : ndarray, shape (num_levels,)
        Water vapor mixing ratio profile [kg/kg]
    actual_levels : int, optional
        Number of actual levels to use (default: len(pressure_levels))

    Returns
    -------
    min_pressure : float
        Minimum central pressure [mb]
    max_wind : float
        Maximum sustained wind speed [m/s]
    error_flag : int
        Error status (0 = success, non-zero = error)

    Examples
    --------
    >>> # Single profile calculation
    >>> min_p, max_w, err = calculate_potential_intensity_profile(
    ...     28.5, 1013.25, p_levels, temp_profile, mixr_profile)
    """
    # Validate and convert inputs (profile data typically doesn't have missing values)
    (pressure_levels, temperature, mixing_ratio), has_missing = _validate_input_arrays(
        pressure_levels,
        temperature,
        mixing_ratio,
        names=["pressure_levels", "temperature", "mixing_ratio"],
    )

    # Ensure correct pressure level ordering (surface to top)
    pressure_levels, temperature, mixing_ratio = _ensure_pressure_ordering(
        pressure_levels, temperature, mixing_ratio
    )

    if actual_levels is None:
        actual_levels = len(pressure_levels)

    # Validate profile dimensions
    expected_len = len(pressure_levels)
    if len(temperature) != expected_len or len(mixing_ratio) != expected_len:
        raise ValueError(
            f"Profile lengths mismatch - pressure: {expected_len}, temperature: {len(temperature)}, mixing_ratio: {len(mixing_ratio)}"
        )

    # Call Fortran function
    min_pressure, max_wind, error_flag = _gpi_module.calculate_pi_single_profile(
        float(sst),
        float(psl),
        pressure_levels,
        temperature,
        mixing_ratio,
        actual_levels,
    )

    # Convert UNDEF to NaN in results (though less common for profile data)
    min_pressure = np.nan if min_pressure == UNDEF else min_pressure
    max_wind = np.nan if max_wind == UNDEF else max_wind

    return min_pressure, max_wind, error_flag


class PotentialIntensityCalculator:
    """
    Class-based interface for tropical cyclone potential intensity calculations.

    This class provides a high-level interface for calculating potential intensity
    with automatic dimension handling and result caching.

    Parameters
    ----------
    sst : ndarray
        Sea surface temperature data [K]
    psl : ndarray
        Sea level pressure data [Pa]
    pressure_levels : ndarray
        Atmospheric pressure levels [mb]
    temperature : ndarray
        Temperature profile data [K]
    mixing_ratio : ndarray
        Water vapor mixing ratio data [kg/kg]

    Examples
    --------
    >>> # Create calculator instance
    >>> pi_calc = PotentialIntensityCalculator(sst, psl, p_levels, temp, mixr)
    >>>
    >>> # Calculate potential intensity
    >>> min_p, max_w, err = pi_calc.calculate()
    >>>
    >>> # Access results
    >>> results = pi_calc.results
    """

    def __init__(
        self,
        sst: np.ndarray,
        psl: np.ndarray,
        pressure_levels: np.ndarray,
        temperature: np.ndarray,
        mixing_ratio: np.ndarray,
    ):
        self.sst = np.asarray(sst)
        self.psl = np.asarray(psl)
        self.pressure_levels = np.asarray(pressure_levels)
        self.temperature = np.asarray(temperature)
        self.mixing_ratio = np.asarray(mixing_ratio)
        self._results = None

        # Auto-detect data type
        self._data_type = self._detect_data_type()

    def _detect_data_type(self) -> str:
        """Detect whether data is 3D, 4D, or single profile based on temperature/mixing_ratio dimensions."""
        # Use temperature dimensions to detect data type since it contains all dimension info
        temp_ndim = self.temperature.ndim
        # (num_levels,) | (num_levels, nlat, nlon) | (ntimes, num_levels, nlat, nlon)
        data_types = {1: "profile", 3: "3D", 4: "4D"}

        if temp_ndim not in data_types:
            raise ValueError(
                f"Unsupported temperature dimensions: {temp_ndim}. Expected 1, 3, or 4 dimensions."
            )
        return data_types[temp_ndim]

    def calculate(self) -> Tuple[np.ndarray, np.ndarray, int]:
        """
        Calculate potential intensity based on input data dimensions.

        Returns
        -------
        min_pressure : ndarray or float
            Minimum central pressure [mb]
        max_wind : ndarray or float
            Maximum sustained wind speed [m/s]
        error_flag : int
            Error status (0 = success, non-zero = error)
        """
        # Direct function calls based on data type
        if self._data_type == "3D":
            result = calculate_potential_intensity_3d(
                self.sst,
                self.psl,
                self.pressure_levels,
                self.temperature,
                self.mixing_ratio,
            )
        elif self._data_type == "4D":
            result = calculate_potential_intensity_4d(
                self.sst,
                self.psl,
                self.pressure_levels,
                self.temperature,
                self.mixing_ratio,
            )
        elif self._data_type == "profile":
            result = calculate_potential_intensity_profile(
                float(self.sst),
                float(self.psl),
                self.pressure_levels,
                self.temperature,
                self.mixing_ratio,
            )
        else:
            raise ValueError(f"Unsupported data type: {self._data_type}")

        self._results = {
            "min_pressure": result[0],
            "max_wind": result[1],
            "error_flag": result[2],
            "data_type": self._data_type,
        }

        return result

    @property
    def results(self) -> Optional[Dict[str, Any]]:
        """Return calculation results if available."""
        return self._results

    @property
    def data_type(self) -> str:
        """Return detected data type."""
        return self._data_type


# Convenience functions for direct calculation


def potential_intensity(
    sst: np.ndarray,
    psl: np.ndarray,
    pressure_levels: np.ndarray,
    temperature: np.ndarray,
    mixing_ratio: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Calculate tropical cyclone potential intensity with automatic dimension detection.

    This is a convenience function that automatically detects input data dimensions
    and calls the appropriate calculation function.

    Parameters
    ----------
    sst : ndarray
        Sea surface temperature [K]
    psl : ndarray
        Sea level pressure [Pa]
    pressure_levels : ndarray
        Atmospheric pressure levels [mb]
    temperature : ndarray
        Temperature data [K]
    mixing_ratio : ndarray
        Water vapor mixing ratio [kg/kg]

    Returns
    -------
    min_pressure : ndarray or float
        Minimum central pressure [mb]
    max_wind : ndarray or float
        Maximum sustained wind speed [m/s]
    error_flag : int
        Error status (0 = success, non-zero = error)
    """
    calc = PotentialIntensityCalculator(
        sst, psl, pressure_levels, temperature, mixing_ratio
    )
    return calc.calculate()
