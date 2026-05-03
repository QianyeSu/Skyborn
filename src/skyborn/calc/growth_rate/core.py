"""Chemke-style baroclinic and barotropic growth-rate diagnostics.

This module prepares one-dimensional atmospheric profiles in Python and
dispatches the expensive normal-mode solves to the compiled Fortran kernels in
``skyborn.calc.growth_rate.growth_rate_kernels``.

References
----------
Chemke, R., and Ming, Y. (2020):
    Large Atmospheric Waves Will Get Stronger, While Small Waves Will Get
    Weaker by the End of the 21st Century. Geophysical Research Letters, 47,
    e2020GL090441.
    https://doi.org/10.1029/2020GL090441

Chemke, R., Zanna, L., Orbe, C., Sentman, L. T., and Polvani, L. M. (2022):
    The Future Intensification of the North Atlantic Winter Storm Track: The
    Key Role of Dynamic Ocean Coupling. Journal of Climate, 35(8), 2407-2421.
    https://doi.org/10.1175/JCLI-D-21-0407.1

Chemke, R. (2022):
    The future poleward shift of Southern Hemisphere summer mid-latitude
    storm tracks stems from ocean coupling. Nature Communications, 13, 2531.
    https://doi.org/10.1038/s41467-022-29392-4

Chemke, R., Ming, Y., and Yuval, J. (2022):
    The intensification of winter mid-latitude storm tracks in the Southern
    Hemisphere. Nature Climate Change, 12, 553-557.
    https://doi.org/10.1038/s41558-022-01368-8
"""

from __future__ import annotations

import numpy as np
import xarray as xr

from skyborn.interp.interpolation import interp_pressure_1d

from ..troposphere import trop_wmo_profile as _trop_wmo_profile
from .growth_rate_kernels import dbaroc_growth_rate_1d as _dbaroc_growth_rate_1d
from .growth_rate_kernels import dbarot_growth_rate_1d as _dbarot_growth_rate_1d

GAS_CONSTANT_DRY = 287.04
HEAT_CAPACITY = 1004.7
KAPPA = GAS_CONSTANT_DRY / HEAT_CAPACITY
REFERENCE_PRESSURE_PA = 100000.0
DEFAULT_SOLVER_LEVELS = 45

ProfileInput = xr.DataArray | np.ndarray
MissingValue = float | int | np.floating | np.integer | None

__all__ = [
    "baroc_growth_rate",
    "barot_growth_rate",
]


def _as_1d_float64(values: ProfileInput, name: str) -> np.ndarray:
    """Convert one-dimensional eager inputs to ``float64`` NumPy arrays."""

    array = np.asarray(getattr(values, "data", values), dtype=np.float64)
    if array.ndim != 1:
        raise ValueError(f"`{name}` must be one-dimensional")
    return np.array(array, dtype=np.float64, copy=True)


def _same_shape_or_raise(*named_arrays: tuple[str, np.ndarray]) -> None:
    """Validate that all named one-dimensional arrays have identical shapes."""

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


def _is_nan_missing_value(value: MissingValue) -> bool:
    """Return True when ``value`` behaves like a NaN missing sentinel."""

    if value is None:
        return True

    if isinstance(value, (float, int, np.floating, np.integer)):
        return bool(np.isnan(value))

    return False


def _return_missing_value(missing_value: MissingValue) -> float:
    """Return a normalized scalar missing value for the public APIs."""

    if _is_nan_missing_value(missing_value):
        return float("nan")
    return float(missing_value)


def _contains_missing(values: np.ndarray, missing_value: MissingValue) -> bool:
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


def _pressure_to_pa(pressure: ProfileInput) -> np.ndarray:
    """Convert a pressure axis from Pa-or-hPa input to Pascals."""

    pressure_values = _as_1d_float64(pressure, "pressure")
    if np.any(pressure_values <= 0.0):
        raise ValueError("pressure values must be strictly positive")

    if np.nanmax(pressure_values) <= 2000.0:
        return pressure_values * 100.0
    return pressure_values


def _target_pressure_to_pa(target_pressure: ProfileInput) -> np.ndarray:
    """Convert explicit target pressure levels from Pa-or-hPa to Pascals."""

    target_values = _as_1d_float64(target_pressure, "target_pressure")
    if np.any(target_values <= 0.0):
        raise ValueError("target pressure values must be strictly positive")

    if np.nanmax(target_values) <= 2000.0:
        return target_values * 100.0
    return target_values


def _normalize_output_units(output_units: str) -> str:
    """Normalize public output-unit aliases to one canonical label."""

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
    """Map public vertical interpolation aliases to ``interp_pressure_1d`` modes."""

    normalized = method.strip().lower().replace("_", "")
    if normalized in {"logp", "log"}:
        return "log"
    if normalized in {"linear", "pressure"}:
        return "linear"
    raise ValueError(
        "`vertical_interp` must be 'log' or 'linear' "
        "(the alias 'logp' is also accepted)"
    )


def _infer_tropopause_pressure_pa(
    temperature: np.ndarray,
    pressure_pa: np.ndarray,
) -> float:
    """Diagnose tropopause pressure from the input sounding via WMO criteria."""

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
            "tropopause pressure for this profile."
        )

    return tropopause_pressure_hpa * 100.0


def _build_solver_pressure_grid_pa(
    pressure_pa: np.ndarray,
    temperature: np.ndarray,
    target_pressure: ProfileInput | None,
    solver_levels: int,
) -> np.ndarray:
    """Build the solver pressure grid in Pascals."""

    if target_pressure is not None:
        target_pressure_pa = _target_pressure_to_pa(target_pressure)
        _require_monotonic(target_pressure_pa, "target_pressure")
        return target_pressure_pa

    tropopause_pressure_pa = _infer_tropopause_pressure_pa(
        temperature,
        pressure_pa,
    )
    bottom_pressure_pa = min(REFERENCE_PRESSURE_PA, float(np.max(pressure_pa)))
    if tropopause_pressure_pa >= bottom_pressure_pa:
        raise ValueError(
            "The diagnosed tropopause pressure must be lower than the lower-"
            "tropospheric pressure bound used by the solver grid."
        )

    return np.linspace(
        tropopause_pressure_pa,
        bottom_pressure_pa,
        solver_levels,
        dtype=np.float64,
    )


def _normalize_solver_levels(solver_levels: int) -> int:
    """Validate the requested automatic solver-grid level count."""

    if isinstance(solver_levels, bool) or not isinstance(
        solver_levels, (int, np.integer)
    ):
        raise TypeError("`solver_levels` must be an integer >= 2")

    normalized = int(solver_levels)
    if normalized < 2:
        raise ValueError("`solver_levels` must be at least 2")
    return normalized


def barot_growth_rate(
    u_barotropic: ProfileInput,
    lat: ProfileInput,
    *,
    output_units: str = "s-1",
    radius: float = 6_371_000.0,
    omega: float = 7.292e-5,
    missing_value: MissingValue = np.nan,
) -> float:
    r"""Compute the maximum barotropic normal-mode growth rate from ``U(lat)``.

    Parameters
    ----------
    u_barotropic : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        One-dimensional vertically averaged zonal-wind profile in ``m s^-1``.

    lat : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        One-dimensional latitude coordinate in degrees north. Values must be
        strictly monotonic.

    output_units : {"s^-1", "day^-1"}, optional
        Output units for the returned maximum growth rate. Defaults to
        ``"s-1"``.

    radius : float, optional
        Earth radius in meters. Only the default value is currently supported
        by the compiled kernel.

    omega : float, optional
        Planetary rotation rate in ``s^-1``. Only the default value is
        currently supported by the compiled kernel.

    missing_value : scalar, optional
        Public missing-value marker. If any input value is missing, the
        function returns this marker normalized to ``float``. Defaults to
        ``np.nan``.

    Returns
    -------
    float
        Maximum barotropic growth rate in the requested units.

    Notes
    -----
    This diagnostic solves the finite-difference generalized eigenproblem

    .. math::

       A(k)\psi = \lambda B(k)\psi

    for each zonal wavenumber :math:`k`, where the continuous barotropic
    operator is

    .. math::

       \left[k\beta - k^3 U + kU\partial_{yy} - kU_{yy}\right]\psi
       = \lambda \left(\partial_{yy} - k^2\right)\psi,

    with

    .. math::

       \beta = \frac{2\Omega\cos(\phi)}{a}.

    The returned value is :math:`\max_k \operatorname{Im}(\lambda_k)`.

    References
    ----------
    Chemke, R., Ming, Y., and Yuval, J. (2022):
        The intensification of winter mid-latitude storm tracks in the
        Southern Hemisphere. Nature Climate Change, 12, 553-557.
        https://doi.org/10.1038/s41558-022-01368-8
    """

    if radius != 6_371_000.0 or omega != 7.292e-5:
        raise NotImplementedError(
            "Custom `radius` and `omega` values are not implemented yet for "
            "the compiled barotropic growth-rate kernel."
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
            "requested latitude profile."
        )

    return _convert_output_units(max_growth, output_units)


def baroc_growth_rate(
    u: ProfileInput,
    temperature: ProfileInput,
    pressure: ProfileInput,
    *,
    lat: float | int | np.floating | np.integer | None = None,
    output_units: str = "s-1",
    vertical_interp: str = "log",
    target_pressure: ProfileInput | None = None,
    solver_levels: int = DEFAULT_SOLVER_LEVELS,
    missing_value: MissingValue = np.nan,
) -> float:
    r"""Compute the maximum baroclinic normal-mode growth rate.

    Parameters
    ----------
    u : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        One-dimensional zonal-wind profile in ``m s^-1``.

    temperature : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        One-dimensional temperature profile in Kelvin.

    pressure : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        One-dimensional pressure coordinate in Pa or hPa. Values must be
        strictly monotonic and strictly positive.

    lat : scalar
        Latitude in degrees north. Required because the compiled solver uses
        the local Coriolis parameter and meridional gradient of planetary
        vorticity.

    output_units : {"s^-1", "day^-1"}, optional
        Output units for the returned maximum growth rate. Defaults to
        ``"s-1"``.

    vertical_interp : {"log", "linear"}, optional
        Vertical interpolation method used to remap the input profiles onto the
        solver pressure grid. ``"log"`` means linear interpolation in
        log-pressure, not logarithmic interpolation of the field itself.
        Defaults to ``"log"``. The legacy alias ``"logp"`` is still accepted.

    target_pressure : :class:`xarray.DataArray`, :class:`numpy.ndarray`, optional
        Optional explicit solver pressure grid in Pa or hPa. If omitted, the
        function diagnoses the WMO tropopause pressure from the input
        ``temperature`` and builds a fixed-pressure solver grid from that
        tropopause to the lower-tropospheric bound.

    solver_levels : int, optional
        Number of pressure levels used when ``target_pressure`` is omitted and
        the solver grid is built automatically. Defaults to
        ``DEFAULT_SOLVER_LEVELS``. Ignored when ``target_pressure`` is given
        explicitly.

    missing_value : scalar, optional
        Public missing-value marker. If any input value is missing, the
        function returns this marker normalized to ``float``. Defaults to
        ``np.nan``.

    Returns
    -------
    float
        Maximum baroclinic growth rate in the requested units.

    Notes
    -----
    When ``target_pressure`` is omitted, this function first diagnoses the
    WMO tropopause pressure and then interpolates ``u`` and ``temperature`` to
    a fixed pressure grid between that tropopause and the lower troposphere.
    The thermodynamic profiles passed to the compiled kernel are

    .. math::

       \theta = T\left(\frac{p_0}{p}\right)^{\kappa},
       \qquad
       \rho = \frac{p}{R_d T},
       \qquad
       N_z = -\frac{1}{\rho\theta}\frac{\partial \theta}{\partial p}.

    The compiled solver then forms the finite-difference generalized
    eigenproblem

    .. math::

       A(k)\psi = \lambda B(k)\psi,

    where the vertical quasi-geostrophic operator is represented in pressure
    coordinates by

    .. math::

       A(k) \sim k\beta - k^3 U
       + kUf^2 \partial_p \left(N_z^{-1}\partial_p\right)
       - kf^2 \left(\partial_p U\right) N_z^{-1}\partial_p,

    and

    .. math::

       B(k) \sim f^2 \partial_p \left(N_z^{-1}\partial_p\right) - k^2,

    and the returned value is :math:`\max_k \operatorname{Im}(\lambda_k)`.

    The default ``DEFAULT_SOLVER_LEVELS = 45`` is a practical resolution
    choice, not a formal optimum. Increasing ``solver_levels`` can improve
    vertical-grid convergence up to a point, but it also increases the
    per-profile eigenvalue cost and does not guarantee a better diagnostic once
    the solution is already converged on the chosen pressure interval.

    References
    ----------
    Chemke, R., and Ming, Y. (2020):
        Large Atmospheric Waves Will Get Stronger, While Small Waves Will Get
        Weaker by the End of the 21st Century. Geophysical Research Letters,
        47, e2020GL090441.
        https://doi.org/10.1029/2020GL090441

    Chemke, R., Zanna, L., Orbe, C., Sentman, L. T., and Polvani, L. M. (2022):
        The Future Intensification of the North Atlantic Winter Storm Track:
        The Key Role of Dynamic Ocean Coupling. Journal of Climate, 35(8),
        2407-2421.
        https://doi.org/10.1175/JCLI-D-21-0407.1

    Chemke, R. (2022):
        The future poleward shift of Southern Hemisphere summer mid-latitude
        storm tracks stems from ocean coupling. Nature Communications, 13,
        2531.
        https://doi.org/10.1038/s41467-022-29392-4

    Chemke, R., Ming, Y., and Yuval, J. (2022):
        The intensification of winter mid-latitude storm tracks in the
        Southern Hemisphere. Nature Climate Change, 12, 553-557.
        https://doi.org/10.1038/s41558-022-01368-8
    """

    if lat is None:
        raise ValueError("`lat` is required for baroc_growth_rate")

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
    solver_levels = _normalize_solver_levels(solver_levels)

    target_pressure_pa = _build_solver_pressure_grid_pa(
        pressure_pa,
        temperature_values,
        target_pressure,
        solver_levels,
    )
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
            "requested atmospheric column."
        )

    return _convert_output_units(max_growth, output_units)
