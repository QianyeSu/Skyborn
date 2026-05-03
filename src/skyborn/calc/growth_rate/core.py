"""Chemke-style baroclinic and barotropic growth-rate diagnostics."""

from __future__ import annotations

import numpy as np

from skyborn.interp.interpolation import interp_pressure_1d

try:
    from ..troposphere import trop_wmo_profile as _trop_wmo_profile
except ImportError:
    _trop_wmo_profile = None

try:
    from .growth_rate_kernels import dbaroc_growth_rate_1d as _dbaroc_growth_rate_1d
    from .growth_rate_kernels import dbarot_growth_rate_1d as _dbarot_growth_rate_1d
except ImportError:
    _dbaroc_growth_rate_1d = None
    _dbarot_growth_rate_1d = None

GAS_CONSTANT_DRY = 287.04
HEAT_CAPACITY = 1004.7
KAPPA = GAS_CONSTANT_DRY / HEAT_CAPACITY
REFERENCE_PRESSURE_PA = 100000.0
DEFAULT_SOLVER_LEVELS = 45

__all__ = [
    "baroc_growth_rate",
    "barot_growth_rate",
]


def _require_growth_rate_backend(function_name: str, handle) -> None:
    """Raise a clear error when the compiled growth-rate backend is unavailable."""

    if handle is None:
        raise RuntimeError(
            f"{function_name} requires the compiled Fortran backend in "
            "`skyborn.calc.growth_rate.growth_rate_kernels`."
        )


def _as_1d_float64(values, name: str) -> np.ndarray:
    """Convert one-dimensional inputs to eager float64 NumPy arrays."""

    array = np.asarray(getattr(values, "data", values), dtype=np.float64)
    if array.ndim != 1:
        raise ValueError(f"`{name}` must be one-dimensional")
    return np.array(array, dtype=np.float64, copy=True)


def _same_shape_or_raise(*named_arrays) -> None:
    """Validate that all named 1D arrays have identical shapes."""

    names = [name for name, _ in named_arrays]
    arrays = [array for _, array in named_arrays]
    reference_shape = arrays[0].shape
    for name, array in zip(names[1:], arrays[1:]):
        if array.shape != reference_shape:
            raise ValueError(
                "dimension mismatch for growth-rate inputs: "
                + ", ".join(
                    f"{item_name}={item_array.shape}"
                    for item_name, item_array in zip(names, arrays)
                )
            )


def _is_nan_missing_value(value) -> bool:
    """Return True when ``value`` behaves like a NaN missing sentinel."""

    if value is None:
        return True

    try:
        return bool(np.isnan(value))
    except TypeError:
        return False


def _return_missing_value(missing_value):
    """Return a normalized scalar missing value for public APIs."""

    if _is_nan_missing_value(missing_value):
        return float("nan")
    return float(missing_value)


def _contains_missing(values: np.ndarray, missing_value) -> bool:
    """Return True when an input profile contains missing values."""

    if np.any(~np.isfinite(values)):
        return True
    if not _is_nan_missing_value(missing_value):
        return bool(np.any(values == float(missing_value)))
    return False


def _require_monotonic(values: np.ndarray, name: str) -> None:
    """Require a strictly monotonic one-dimensional coordinate."""

    diffs = np.diff(values)
    if not (np.all(diffs > 0.0) or np.all(diffs < 0.0)):
        raise ValueError(f"`{name}` must be strictly monotonic")


def _pressure_to_pa(pressure) -> np.ndarray:
    """Convert a pressure axis from Pa-or-hPa input to Pascals."""

    pressure_values = _as_1d_float64(pressure, "pressure")
    if np.any(pressure_values <= 0.0):
        raise ValueError("pressure values must be strictly positive")

    if np.nanmax(pressure_values) <= 2000.0:
        return pressure_values * 100.0
    return pressure_values


def _target_pressure_to_pa(target_pressure) -> np.ndarray:
    """Convert explicit target pressure levels to Pascals."""

    target_values = _as_1d_float64(target_pressure, "target_pressure")
    if np.any(target_values <= 0.0):
        raise ValueError("target pressure values must be strictly positive")

    if np.nanmax(target_values) <= 2000.0:
        return target_values * 100.0
    return target_values


def _normalize_output_units(output_units: str) -> str:
    """Normalize public output-unit labels."""

    normalized = output_units.strip().lower()
    if normalized in {"s^-1", "s-1", "1/s"}:
        return "s^-1"
    if normalized in {"day^-1", "day-1", "1/day"}:
        return "day^-1"
    raise ValueError("`output_units` must be 's^-1' or 'day^-1'")


def _convert_output_units(value: float, output_units: str) -> float:
    """Convert per-second growth rates into the requested public units."""

    normalized = _normalize_output_units(output_units)
    if normalized == "day^-1":
        return float(value) * 86400.0
    return float(value)


def _normalize_vertical_interp(method: str) -> str:
    """Map public vertical interpolation aliases to interp_pressure_1d modes."""

    normalized = method.strip().lower().replace("_", "")
    if normalized in {"logp", "log"}:
        return "log"
    if normalized in {"linear", "pressure"}:
        return "linear"
    raise ValueError("`vertical_interp` must be 'logp' or 'linear'")


def _infer_tropopause_pressure_pa(
    temperature: np.ndarray,
    pressure_pa: np.ndarray,
) -> float:
    """Diagnose tropopause pressure from the input sounding via WMO criteria."""

    if _trop_wmo_profile is None:
        raise RuntimeError(
            "Automatic tropopause detection requires the compiled "
            "`skyborn.calc.troposphere` backend. Provide "
            "`tropopause_pressure` explicitly if that backend is unavailable."
        )

    pressure_for_tropopause = np.array(pressure_pa, dtype=np.float64, copy=True)
    temperature_for_tropopause = np.array(temperature, dtype=np.float64, copy=True)
    if pressure_for_tropopause[0] > pressure_for_tropopause[-1]:
        pressure_for_tropopause = pressure_for_tropopause[::-1].copy()
        temperature_for_tropopause = temperature_for_tropopause[::-1].copy()

    result = _trop_wmo_profile(
        temperature_for_tropopause,
        pressure_for_tropopause,
        pressure_unit="Pa",
    )
    tropopause_pressure_hpa = float(result["pressure"])
    success = bool(result["success"])
    if (
        not success
        or not np.isfinite(tropopause_pressure_hpa)
        or tropopause_pressure_hpa <= 0.0
    ):
        raise RuntimeError(
            "Automatic WMO tropopause detection did not return a usable "
            "tropopause pressure for this profile; provide "
            "`tropopause_pressure` explicitly."
        )

    return tropopause_pressure_hpa * 100.0


def barot_growth_rate(
    u_barotropic,
    lat,
    *,
    output_units: str = "s-1",
    radius: float = 6_371_000.0,
    omega: float = 7.292e-5,
    missing_value=np.nan,
):
    """Compute maximum barotropic growth rate from ``U(lat)``.

    Notes
    -----
    This first implementation follows the Chemke-style one-dimensional
    normal-mode problem directly. It expects a single latitude profile and does
    not yet provide xarray batch dispatch.
    """

    _require_growth_rate_backend("barot_growth_rate", _dbarot_growth_rate_1d)
    if radius != 6_371_000.0 or omega != 7.292e-5:
        raise NotImplementedError(
            "custom `radius` and `omega` values are not implemented yet for "
            "the compiled barotropic growth-rate kernel"
        )

    u_values = _as_1d_float64(u_barotropic, "u_barotropic")
    lat_values = _as_1d_float64(lat, "lat")
    _same_shape_or_raise(("u_barotropic", u_values), ("lat", lat_values))

    missing = _return_missing_value(missing_value)
    if _contains_missing(u_values, missing_value) or _contains_missing(
        lat_values, missing_value
    ):
        return missing

    _require_monotonic(lat_values, "lat")
    if lat_values[0] > lat_values[-1]:
        lat_values = lat_values[::-1].copy()
        u_values = u_values[::-1].copy()

    max_growth, ier = _dbarot_growth_rate_1d(lat_values, u_values)
    if ier != 0:
        raise RuntimeError(
            f"barot_growth_rate Fortran backend returned ier={ier} for the "
            "requested latitude profile"
        )

    return _convert_output_units(max_growth, output_units)


def baroc_growth_rate(
    u,
    temperature,
    pressure,
    *,
    lat=None,
    tropopause_pressure=None,
    tropopause_height=None,
    output_units: str = "s-1",
    vertical_interp: str = "logp",
    target_pressure=None,
    missing_value=np.nan,
):
    """Compute maximum baroclinic normal-mode growth rate.

    Notes
    -----
    This first implementation expects one atmospheric column. The expensive
    normal-mode solve runs only through the compiled Fortran backend; the Python
    wrapper handles unit normalization and pressure-grid preparation.
    """

    _require_growth_rate_backend("baroc_growth_rate", _dbaroc_growth_rate_1d)

    if lat is None:
        raise ValueError("`lat` is required for baroc_growth_rate")
    if tropopause_height is not None and tropopause_pressure is None:
        raise NotImplementedError(
            "`tropopause_height` cannot yet be converted directly into the "
            "solver pressure bound. Omit it to use automatic WMO tropopause "
            "detection, or provide `tropopause_pressure` explicitly."
        )

    u_values = _as_1d_float64(u, "u")
    temperature_values = _as_1d_float64(temperature, "temperature")
    pressure_pa = _pressure_to_pa(pressure)
    _same_shape_or_raise(
        ("u", u_values),
        ("temperature", temperature_values),
        ("pressure", pressure_pa),
    )

    missing = _return_missing_value(missing_value)
    if (
        _contains_missing(u_values, missing_value)
        or _contains_missing(temperature_values, missing_value)
        or _contains_missing(pressure_pa, missing_value)
    ):
        return missing

    if np.any(temperature_values <= 0.0):
        raise ValueError("temperature values must be strictly positive")

    _require_monotonic(pressure_pa, "pressure")

    if target_pressure is None:
        if tropopause_pressure is None:
            tropopause_pressure_pa = _infer_tropopause_pressure_pa(
                temperature_values,
                pressure_pa,
            )
        else:
            tropopause_pressure_pa = float(
                _target_pressure_to_pa([tropopause_pressure])[0]
            )

        top_pressure = tropopause_pressure_pa
        bottom_pressure = min(REFERENCE_PRESSURE_PA, float(np.max(pressure_pa)))
        if top_pressure >= bottom_pressure:
            raise ValueError(
                "tropopause pressure must be lower than the lower-tropospheric "
                "pressure bound used by the solver grid"
            )
        target_pressure_pa = np.linspace(
            top_pressure,
            bottom_pressure,
            DEFAULT_SOLVER_LEVELS,
            dtype=np.float64,
        )
    else:
        target_pressure_pa = _target_pressure_to_pa(target_pressure)

    _require_monotonic(target_pressure_pa, "target_pressure")
    interp_method = _normalize_vertical_interp(vertical_interp)
    u_solver = np.asarray(
        interp_pressure_1d(
            u_values,
            pressure_pa,
            target_pressure_pa,
            method=interp_method,
            extrapolate=False,
            missing_value=np.nan,
        ),
        dtype=np.float64,
    )
    temperature_solver = np.asarray(
        interp_pressure_1d(
            temperature_values,
            pressure_pa,
            target_pressure_pa,
            method=interp_method,
            extrapolate=False,
            missing_value=np.nan,
        ),
        dtype=np.float64,
    )

    if np.any(~np.isfinite(u_solver)) or np.any(~np.isfinite(temperature_solver)):
        return missing

    if target_pressure_pa[0] > target_pressure_pa[-1]:
        target_pressure_pa = target_pressure_pa[::-1].copy()
        u_solver = u_solver[::-1].copy()
        temperature_solver = temperature_solver[::-1].copy()

    theta_solver = (
        temperature_solver * (REFERENCE_PRESSURE_PA / target_pressure_pa) ** KAPPA
    )
    max_growth, ier = _dbaroc_growth_rate_1d(
        u_solver,
        theta_solver,
        target_pressure_pa,
        temperature_solver,
        float(lat),
    )
    if ier != 0:
        raise RuntimeError(
            f"baroc_growth_rate Fortran backend returned ier={ier} for the "
            "requested atmospheric column"
        )

    return _convert_output_units(max_growth, output_units)
