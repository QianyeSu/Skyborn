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
from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple, Union

import numpy as np

from . import tropical_cyclone_potential_intensity as _gpi_module

# Fortran UNDEF value constant
UNDEF = -9.99e33
_OUTFLOW_SOURCE_FLAGS = {"cape_star": 0, "cape_env": 1}
_TEMPERATURE_KIND_BY_NDIM = {1: "profile", 3: "3D", 4: "4D"}
_SSTUnits = str

_GRIDDED_BACKENDS = {
    ("3D", False, False): _gpi_module.calculate_pi_gridded_data,
    ("3D", False, True): _gpi_module.calculate_pi_gridded_with_missing,
    ("3D", True, False): _gpi_module.calculate_pi_gridded_diagnostics,
    ("3D", True, True): _gpi_module.calculate_pi_gridded_diagnostics_with_missing,
    ("4D", False, False): _gpi_module.calculate_pi_4d_data,
    ("4D", False, True): _gpi_module.calculate_pi_4d_with_missing,
    ("4D", True, False): _gpi_module.calculate_pi_4d_diagnostics,
    ("4D", True, True): _gpi_module.calculate_pi_4d_diagnostics_with_missing,
}
_PROFILE_BACKENDS = {
    False: _gpi_module.calculate_pi_single_profile,
    True: _gpi_module.calculate_pi_profile_diagnostics,
}


@dataclass(frozen=True)
class _PreparedInputs:
    """Normalized inputs ready for a compiled backend call."""

    kind: str
    sst: Any
    psl: Any
    pressure_levels: np.ndarray
    temperature: np.ndarray
    mixing_ratio: np.ndarray
    has_missing: bool = False
    actual_levels: Optional[int] = None


def _postprocess_results(min_pressure, max_wind):
    """Convert UNDEF values to NaN in results."""
    min_pressure = np.where(np.isclose(min_pressure, UNDEF), np.nan, min_pressure)
    max_wind = np.where(np.isclose(max_wind, UNDEF), np.nan, max_wind)
    return min_pressure, max_wind


def _postprocess_scalar(value: Union[float, np.floating]) -> float:
    """Convert scalar UNDEF output to NaN."""
    value = float(value)
    return np.nan if np.isclose(value, UNDEF) else value


def _detect_temperature_kind(temperature: Any) -> str:
    """Map thermodynamic input dimensionality to the supported API kind."""
    temp_ndim = np.ndim(temperature)
    try:
        return _TEMPERATURE_KIND_BY_NDIM[temp_ndim]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported temperature dimensions: {temp_ndim}. Expected 1, 3, or 4 dimensions."
        ) from exc


def _normalize_outflow_source(outflow_source: str) -> int:
    """Map the public outflow-source string to the Fortran flag."""
    try:
        return _OUTFLOW_SOURCE_FLAGS[outflow_source]
    except KeyError as exc:
        raise ValueError(
            f"Invalid outflow_source={outflow_source!r}; expected 'cape_star' or 'cape_env'."
        ) from exc


def _sst_to_kelvin(sst: Any, sst_units: _SSTUnits) -> Any:
    """Convert sea-surface temperature to Kelvin."""
    units = sst_units.upper()
    if units == "K":
        return sst
    if units == "C":
        return np.asarray(sst) + 273.15
    raise ValueError(f"Unsupported sst_units={sst_units!r}; expected 'C' or 'K'.")


def _compiled_float32_array(values: Any, *, order: Optional[str] = None) -> np.ndarray:
    """Materialize a float32 array with an explicit compiled-boundary layout."""
    arr = np.asarray(values, dtype=np.float32)

    if order == "F":
        if arr.ndim <= 1:
            return np.ascontiguousarray(arr)
        return arr if arr.flags.f_contiguous else np.asfortranarray(arr)

    if order == "C":
        return arr if arr.flags.c_contiguous else np.ascontiguousarray(arr)

    if arr.flags.c_contiguous or arr.flags.f_contiguous:
        return arr
    return np.ascontiguousarray(arr)


def _diagnostic_result(
    min_pressure,
    max_wind,
    error_flag,
    outflow_temp,
    outflow_level,
    lnpi,
    lneff,
    lndiseq,
    lnCKCD,
) -> Dict[str, Any]:
    """Normalize compiled diagnostic outputs to Python scalars/arrays."""
    if np.asarray(min_pressure).ndim == 0:
        return {
            "max_wind": _postprocess_scalar(max_wind),
            "min_pressure": _postprocess_scalar(min_pressure),
            "error_flag": int(error_flag),
            "t0": _postprocess_scalar(outflow_temp),
            "otl": _postprocess_scalar(outflow_level),
            "lnpi": _postprocess_scalar(lnpi),
            "lneff": _postprocess_scalar(lneff),
            "lndiseq": _postprocess_scalar(lndiseq),
            "lnCKCD": _postprocess_scalar(lnCKCD),
        }

    min_pressure, max_wind = _postprocess_results(min_pressure, max_wind)
    outflow_temp = np.where(np.isclose(outflow_temp, UNDEF), np.nan, outflow_temp)
    outflow_level = np.where(np.isclose(outflow_level, UNDEF), np.nan, outflow_level)
    lnpi = np.where(np.isclose(lnpi, UNDEF), np.nan, lnpi)
    lneff = np.where(np.isclose(lneff, UNDEF), np.nan, lneff)
    lndiseq = np.where(np.isclose(lndiseq, UNDEF), np.nan, lndiseq)

    return {
        "max_wind": max_wind,
        "min_pressure": min_pressure,
        "error_flag": int(error_flag),
        "t0": outflow_temp,
        "otl": outflow_level,
        "lnpi": lnpi,
        "lneff": lneff,
        "lndiseq": lndiseq,
        "lnCKCD": _postprocess_scalar(lnCKCD),
    }


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

    # Already correctly ordered or single level - return original arrays
    return pressure_levels, temperature, mixing_ratio


def _validate_dimensions(sst, psl, pressure_levels, temp, mixing_ratio, data_type="3D"):
    """Validate that input arrays have compatible dimensions based on temperature array."""
    expected_levels = len(pressure_levels)

    if data_type == "profile":
        if temp.shape != (expected_levels,):
            raise ValueError(
                f"Temperature shape {temp.shape} doesn't match expected profile shape ({expected_levels},)"
            )
        if not (np.isscalar(sst) or sst.ndim == 0) or not (
            np.isscalar(psl) or psl.ndim == 0
        ):
            raise ValueError("SST and PSL must be scalars for profile data")

    elif data_type == "3D":
        if temp.ndim != 3 or temp.shape[0] != expected_levels:
            raise ValueError(
                f"Temperature shape {temp.shape} doesn't match expected 3D shape ({expected_levels}, nlat, nlon)"
            )
        expected_sst_shape = temp.shape[1:]
        if sst.shape != expected_sst_shape or psl.shape != expected_sst_shape:
            raise ValueError(
                f"SST/PSL shape mismatch - expected {expected_sst_shape}, got SST:{sst.shape}, PSL:{psl.shape}"
            )

    elif data_type == "4D":
        if temp.ndim != 4 or temp.shape[1] != expected_levels:
            raise ValueError(
                f"Temperature shape {temp.shape} doesn't match expected 4D shape (ntimes, {expected_levels}, nlat, nlon)"
            )
        expected_sst_shape = (temp.shape[0], temp.shape[2], temp.shape[3])
        if sst.shape != expected_sst_shape or psl.shape != expected_sst_shape:
            raise ValueError(
                f"SST/PSL shape mismatch - expected {expected_sst_shape}, got SST:{sst.shape}, PSL:{psl.shape}"
            )

    else:
        raise ValueError(f"Unsupported data type: {data_type}")

    if temp.shape != mixing_ratio.shape:
        raise ValueError(
            f"Temperature shape {temp.shape} doesn't match mixing ratio shape {mixing_ratio.shape}"
        )


def _prepare_gridded_inputs(
    kind: str,
    sst: np.ndarray,
    psl: np.ndarray,
    pressure_levels: np.ndarray,
    temperature: np.ndarray,
    mixing_ratio: np.ndarray,
) -> _PreparedInputs:
    """Validate, reorder, and materialize 3D/4D inputs for the Fortran backend."""
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

    pressure_levels, temperature, mixing_ratio = _ensure_pressure_ordering(
        pressure_levels, temperature, mixing_ratio
    )
    _validate_dimensions(sst, psl, pressure_levels, temperature, mixing_ratio, kind)

    # Multi-dimensional inputs from xarray/netCDF are usually C contiguous.
    # f2py would need Fortran-layout copies for these signatures anyway, so we
    # materialize the layout explicitly once at the Python boundary.
    return _PreparedInputs(
        kind=kind,
        sst=_compiled_float32_array(sst, order="F"),
        psl=_compiled_float32_array(psl, order="F"),
        pressure_levels=_compiled_float32_array(pressure_levels),
        temperature=_compiled_float32_array(temperature, order="F"),
        mixing_ratio=_compiled_float32_array(mixing_ratio, order="F"),
        has_missing=has_missing,
    )


def _prepare_profile_inputs(
    sst: float,
    psl: float,
    pressure_levels: np.ndarray,
    temperature: np.ndarray,
    mixing_ratio: np.ndarray,
    actual_levels: Optional[int] = None,
) -> _PreparedInputs:
    """Validate, reorder, and materialize single-profile inputs."""
    expected_len = len(np.asarray(pressure_levels))
    if len(temperature) != expected_len or len(mixing_ratio) != expected_len:
        raise ValueError(
            f"Profile lengths mismatch - pressure: {expected_len}, temperature: {len(temperature)}, mixing_ratio: {len(mixing_ratio)}"
        )

    (pressure_levels, temperature, mixing_ratio), _has_missing = _validate_input_arrays(
        pressure_levels,
        temperature,
        mixing_ratio,
        names=["pressure_levels", "temperature", "mixing_ratio"],
    )

    pressure_levels, temperature, mixing_ratio = _ensure_pressure_ordering(
        pressure_levels, temperature, mixing_ratio
    )
    _validate_dimensions(
        np.asarray(sst),
        np.asarray(psl),
        pressure_levels,
        temperature,
        mixing_ratio,
        "profile",
    )

    if actual_levels is None:
        actual_levels = len(pressure_levels)
    if actual_levels < 1 or actual_levels > len(pressure_levels):
        raise ValueError(
            f"actual_levels must be between 1 and {len(pressure_levels)}, got {actual_levels}"
        )

    return _PreparedInputs(
        kind="profile",
        sst=float(sst),
        psl=float(psl),
        pressure_levels=_compiled_float32_array(pressure_levels),
        temperature=_compiled_float32_array(temperature),
        mixing_ratio=_compiled_float32_array(mixing_ratio),
        actual_levels=int(actual_levels),
    )


def _prepare_inputs(
    kind: str,
    sst: Any,
    psl: Any,
    pressure_levels: np.ndarray,
    temperature: np.ndarray,
    mixing_ratio: np.ndarray,
    actual_levels: Optional[int] = None,
) -> _PreparedInputs:
    """Dispatch to the appropriate input-normalization path."""
    if kind == "profile":
        return _prepare_profile_inputs(
            float(sst),
            float(psl),
            pressure_levels,
            temperature,
            mixing_ratio,
            actual_levels=actual_levels,
        )

    return _prepare_gridded_inputs(
        kind,
        sst,
        psl,
        pressure_levels,
        temperature,
        mixing_ratio,
    )


def _run_backend(
    prepared: _PreparedInputs,
    *,
    diagnostics: bool = False,
    outflow_source: str = "cape_star",
    CKCD: float = 0.9,
):
    """Execute the requested compiled backend after normalization."""
    if prepared.kind == "profile":
        func = _PROFILE_BACKENDS[diagnostics]
        if diagnostics:
            outflow_flag = _normalize_outflow_source(outflow_source)
            return _diagnostic_result(
                *func(
                    float(prepared.sst),
                    float(prepared.psl),
                    prepared.pressure_levels,
                    prepared.temperature,
                    prepared.mixing_ratio,
                    prepared.actual_levels,
                    outflow_source_flag=outflow_flag,
                    ckcd_in=float(CKCD),
                )
            )

        min_pressure, max_wind, error_flag = func(
            float(prepared.sst),
            float(prepared.psl),
            prepared.pressure_levels,
            prepared.temperature,
            prepared.mixing_ratio,
            prepared.actual_levels,
        )
        return (
            _postprocess_scalar(min_pressure),
            _postprocess_scalar(max_wind),
            int(error_flag),
        )

    func = _GRIDDED_BACKENDS[(prepared.kind, diagnostics, prepared.has_missing)]
    if diagnostics:
        outflow_flag = _normalize_outflow_source(outflow_source)
        return _diagnostic_result(
            *func(
                prepared.sst,
                prepared.psl,
                prepared.pressure_levels,
                prepared.temperature,
                prepared.mixing_ratio,
                outflow_source_flag=outflow_flag,
                ckcd_in=float(CKCD),
            )
        )

    min_pressure, max_wind, error_flag = func(
        prepared.sst,
        prepared.psl,
        prepared.pressure_levels,
        prepared.temperature,
        prepared.mixing_ratio,
    )
    min_pressure, max_wind = _postprocess_results(min_pressure, max_wind)
    return min_pressure, max_wind, int(error_flag)


def log_decompose_pi(
    pi: Any,
    sst: Any,
    t0: Any,
    CKCD: float = 0.9,
    *,
    sst_units: _SSTUnits = "K",
) -> Tuple[Any, Any, Any, float]:
    """Log-decompose potential intensity into efficiency, disequilibrium, and Ck/Cd.

    Parameters
    ----------
    pi
        Potential intensity wind speed [m/s].
    sst
        Sea surface temperature in units given by ``sst_units``.
    t0
        Outflow temperature [K].
    CKCD : float, default: 0.9
        Ratio of exchange coefficients.
    sst_units : {"K", "C"}, default: "K"
        Units of ``sst``.

    Returns
    -------
    tuple
        ``(lnpi, lneff, lndiseq, lnCKCD)`` where ``lnpi = ln(V^2)``.
    """
    pi_arr = np.asarray(pi, dtype=float)
    t0_arr = np.asarray(t0, dtype=float)
    sst_k = np.asarray(_sst_to_kelvin(sst, sst_units), dtype=float)

    pi_arr, sst_k, t0_arr = np.broadcast_arrays(pi_arr, sst_k, t0_arr)
    lnCKCD = float(np.log(CKCD))

    efficiency = (sst_k - t0_arr) / t0_arr
    valid_eff = efficiency > 0
    valid_pi = pi_arr > 0

    lneff = np.full(efficiency.shape, np.nan, dtype=float)
    lneff[valid_eff] = np.log(efficiency[valid_eff])

    lnpi = np.full(pi_arr.shape, np.nan, dtype=float)
    lnpi[valid_pi] = 2.0 * np.log(pi_arr[valid_pi])

    lndiseq = np.full(pi_arr.shape, np.nan, dtype=float)
    valid = valid_eff & valid_pi
    lndiseq[valid] = lnpi[valid] - lneff[valid] - lnCKCD

    return (
        float(lnpi) if np.ndim(lnpi) == 0 else lnpi,
        float(lneff) if np.ndim(lneff) == 0 else lneff,
        float(lndiseq) if np.ndim(lndiseq) == 0 else lndiseq,
        lnCKCD,
    )


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
        Error status (`1` = success, other values indicate non-convergence or invalid input)

    Examples
    --------
    >>> # 3D calculation with automatic missing value detection
    >>> min_p, max_w, err = calculate_potential_intensity_3d(
    ...     sst, psl, p_levels, temp, mixr)
    """
    prepared = _prepare_inputs(
        "3D", sst, psl, pressure_levels, temperature, mixing_ratio
    )
    return _run_backend(prepared)


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
        Error status (`1` = success, other values indicate non-convergence or invalid input)

    Examples
    --------
    >>> # 4D time series calculation
    >>> min_p, max_w, err = calculate_potential_intensity_4d(
    ...     sst_4d, psl_4d, p_levels, temp_4d, mixr_4d)
    """
    prepared = _prepare_inputs(
        "4D", sst, psl, pressure_levels, temperature, mixing_ratio
    )
    return _run_backend(prepared)


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
        Number of actual levels to use (default: ``len(pressure_levels)``).

    Returns
    -------
    tuple
        ``(min_pressure, max_wind, error_flag)``.
    """
    prepared = _prepare_inputs(
        "profile",
        sst,
        psl,
        pressure_levels,
        temperature,
        mixing_ratio,
        actual_levels=actual_levels,
    )
    return _run_backend(prepared)


def pi_log_decomposition(
    sst: Any,
    psl: Any,
    pressure_levels: np.ndarray,
    temperature: np.ndarray,
    mixing_ratio: np.ndarray,
    CKCD: float = 0.9,
    *,
    outflow_source: str = "cape_star",
) -> Dict[str, Any]:
    """Calculate PI plus the Wing et al. (2015) logarithmic decomposition.

    Parameters
    ----------
    sst
        Sea surface temperature [K].
    psl
        Sea level pressure [Pa].
    pressure_levels : ndarray
        Atmospheric pressure levels [mb].
    temperature : ndarray
        Thermodynamic input field. Supported shapes are ``(level,)``,
        ``(level, y, x)``, and ``(time, level, y, x)``.
    mixing_ratio : ndarray
        Water vapor mixing ratio [kg/kg] with the same shape as
        ``temperature``.
    CKCD : float, default: 0.9
        Ratio of exchange coefficients.
    outflow_source : {"cape_star", "cape_env"}, default: "cape_star"
        Outflow branch used for the backend diagnostics.

    Returns
    -------
    dict
        ``max_wind``, ``min_pressure``, ``error_flag``, ``t0``, ``otl``,
        ``lnpi``, ``lneff``, ``lndiseq``, and ``lnCKCD``.
    """
    kind = _detect_temperature_kind(temperature)
    prepared = _prepare_inputs(
        kind, sst, psl, pressure_levels, temperature, mixing_ratio
    )
    return _run_backend(
        prepared,
        diagnostics=True,
        outflow_source=outflow_source,
        CKCD=CKCD,
    )


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
        Error status (`1` = success, other values indicate non-convergence or invalid input)
    """
    kind = _detect_temperature_kind(temperature)
    prepared = _prepare_inputs(
        kind, sst, psl, pressure_levels, temperature, mixing_ratio
    )
    return _run_backend(prepared)
