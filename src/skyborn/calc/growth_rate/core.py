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

import warnings
from typing import Optional, Tuple, Union

import numpy as np
import xarray as xr

from skyborn.interp.interpolation import interp_pressure_1d

from ..troposphere import trop_wmo as _trop_wmo
from ..troposphere import trop_wmo_profile as _trop_wmo_profile
from .growth_rate_kernels import dbaroc_growth_rate_1d as _dbaroc_growth_rate_1d
from .growth_rate_kernels import (
    dbaroc_growth_rate_profiles as _dbaroc_growth_rate_profiles,
)
from .growth_rate_kernels import dbarot_growth_rate_1d as _dbarot_growth_rate_1d

GAS_CONSTANT_DRY = 287.04
HEAT_CAPACITY = 1004.7
KAPPA = GAS_CONSTANT_DRY / HEAT_CAPACITY
REFERENCE_PRESSURE_PA = 100000.0
DEFAULT_SOLVER_LEVELS = 45
DEFAULT_SMOOTH_WINDOW = 1
DEFAULT_WAVENUMBER_MODE = "high"
DEFAULT_HIGH_RES_WAVENUMBER_COUNT = 200
DEFAULT_HIGH_RES_WAVENUMBER_STEP = 1.0e-7
DEFAULT_LOW_RES_LONGITUDE = np.arange(0.0, 361.0, 1.5, dtype=np.float64)

ProfileInput = Union[xr.DataArray, np.ndarray]
LatitudeInput = Union[ProfileInput, float, int, np.floating, np.integer]
PressureInput = Union[ProfileInput, float, int, np.floating, np.integer]
GrowthRateOutput = Union[float, np.ndarray, xr.DataArray]

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


def _coerce_profile_matrix(
    values: ProfileInput,
    name: str,
    pressure_size: int,
    pressure_dim: Optional[str],
) -> np.ndarray:
    """Normalize profile inputs to a ``(nprofile, nlev)`` float64 matrix."""

    array = np.asarray(getattr(values, "data", values), dtype=np.float64)
    if array.ndim not in {1, 2}:
        raise ValueError(f"`{name}` must be one- or two-dimensional")
    array = np.array(array, dtype=np.float64, copy=True)
    if array.ndim == 1:
        if array.shape[0] != pressure_size:
            raise ValueError(f"`{name}` length must match the pressure coordinate")
        return array[np.newaxis, :]

    if isinstance(values, xr.DataArray) and pressure_dim in values.dims:
        vertical_axis = values.dims.index(pressure_dim)  # type: ignore[arg-type]
    else:
        candidate_axes = [
            axis for axis, size in enumerate(array.shape) if size == pressure_size
        ]
        if len(candidate_axes) != 1:
            raise ValueError(
                f"`{name}` must have exactly one axis matching the pressure length"
            )
        vertical_axis = candidate_axes[0]

    if array.shape[vertical_axis] != pressure_size:
        raise ValueError(f"`{name}` vertical axis must match the pressure coordinate")

    return np.array(np.moveaxis(array, vertical_axis, -1), dtype=np.float64, copy=True)


def _broadcast_profile_matrix(
    matrix: np.ndarray,
    nprofile: int,
    name: str,
) -> np.ndarray:
    """Broadcast a ``(1, nlev)`` climatology to ``nprofile`` profiles."""

    if matrix.shape[0] == nprofile:
        return np.array(matrix, dtype=np.float64, copy=True)
    if matrix.shape[0] == 1:
        return np.array(
            np.broadcast_to(matrix, (nprofile, matrix.shape[1])),
            dtype=np.float64,
            copy=True,
        )
    raise ValueError(
        f"`{name}` must either have one profile or match the requested profile count"
    )


def _as_profile_vector(
    values: PressureInput,
    name: str,
    nprofile: int,
) -> np.ndarray:
    """Normalize scalar-or-vector inputs to a length-``nprofile`` array."""

    array = np.asarray(getattr(values, "data", values), dtype=np.float64)
    if array.ndim == 0:
        return np.full(nprofile, float(array), dtype=np.float64)
    if array.ndim != 1:
        raise ValueError(f"`{name}` must be scalar or one-dimensional")
    if array.size == nprofile:
        return np.array(array, dtype=np.float64, copy=True)
    if array.size == 1:
        return np.full(nprofile, float(array[0]), dtype=np.float64)
    raise ValueError(f"`{name}` must have length 1 or match the profile count")


def _missing_profile_mask(values: np.ndarray) -> np.ndarray:
    """Return a per-profile missing mask for ``(nprofile, nlev)`` matrices."""

    return np.any(~np.isfinite(values), axis=1)


def _missing_vector_mask(values: np.ndarray) -> np.ndarray:
    """Return a missing-value mask for profile-wise scalar vectors."""

    return ~np.isfinite(values)


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


def _normalize_method(method: str) -> str:
    """Map public vertical interpolation aliases to ``interp_pressure_1d`` modes."""

    normalized = method.strip().lower().replace("_", "")
    if normalized in {"logp", "log"}:
        return "log"
    if normalized in {"linear", "pressure"}:
        return "linear"
    raise ValueError(
        "`method` must be 'log' or 'linear' " "(the alias 'logp' is also accepted)"
    )


def _normalize_wavenumber_mode(wavenumber_mode: str) -> str:
    """Normalize public wavenumber-resolution aliases."""

    normalized = wavenumber_mode.strip().lower().replace("_", "").replace("-", "")
    if normalized in {"high", "highres", "highresolution"}:
        return "high"
    if normalized in {"low", "lowres", "lowresolution"}:
        return "low"
    raise ValueError("`wavenumber_mode` must be 'high' or 'low'")


def _f_beta_from_lat_bounds(lat_bounds: Tuple[float, float]) -> Tuple[float, float]:
    """Return Chemke-style band-mean Coriolis terms for a latitude range."""

    lat0, lat1 = float(lat_bounds[0]), float(lat_bounds[1])
    sin_avg = (np.sin(np.deg2rad(lat0)) + np.sin(np.deg2rad(lat1))) / 2.0
    cos_avg = (np.cos(np.deg2rad(lat0)) + np.cos(np.deg2rad(lat1))) / 2.0
    f_cor = 2.0 * 7.292e-5 * sin_avg
    beta = 2.0 * 7.292e-5 * cos_avg / 6_371_000.0
    return float(f_cor), float(beta)


def _prepare_wavenumber_inputs(
    wavenumber_mode: str,
    lon: Optional[ProfileInput],
) -> Tuple[int, int, Optional[float]]:
    """Return normalized wavenumber-mode inputs for the compiled backend."""

    normalized_mode = _normalize_wavenumber_mode(wavenumber_mode)
    if normalized_mode == "high":
        return 1, DEFAULT_HIGH_RES_WAVENUMBER_COUNT, None

    lon_values = (
        np.array(DEFAULT_LOW_RES_LONGITUDE, dtype=np.float64, copy=True)
        if lon is None
        else _as_1d_float64(lon, "lon")
    )
    _require_monotonic(lon_values, "lon")
    lon_span_degrees = abs(float(lon_values[0]) - float(lon_values[-1]))
    if lon_span_degrees <= 0.0:
        raise ValueError("`lon` must span a nonzero zonal distance")

    wavenumber_count = int(round(lon_values.size / 3.0))
    if wavenumber_count < 2:
        raise ValueError(
            "`lon` must contain enough points to build at least two low-resolution "
            "zonal wavenumbers"
        )

    return 2, wavenumber_count, lon_span_degrees * np.pi / 180.0


def _reduce_latitude_band_inputs(
    u: ProfileInput,
    temperature: ProfileInput,
    lat: Optional[LatitudeInput],
    lat_bounds: Tuple[float, float],
    pressure_pa: np.ndarray,
    pressure_dim: Optional[str],
) -> Tuple[ProfileInput, ProfileInput, Optional[LatitudeInput]]:
    """Collapse explicit latitude axes using a cosine-weighted band mean."""

    lat_min, lat_max = sorted((float(lat_bounds[0]), float(lat_bounds[1])))
    latitude_dim = None

    for values in (u, temperature):
        if not isinstance(values, xr.DataArray):
            continue

        candidate_dims = [
            dim
            for dim in values.dims
            if dim != pressure_dim and dim.lower() in {"lat", "latitude"}
        ]
        if len(candidate_dims) > 1:
            raise ValueError(
                "Latitude-band growth-rate inputs must have at most one "
                "latitude dimension per DataArray"
            )
        if not candidate_dims:
            continue
        if latitude_dim is None:
            latitude_dim = candidate_dims[0]
            continue
        if candidate_dims[0] != latitude_dim:
            raise ValueError(
                "Latitude-band xarray inputs must share the same latitude " "dimension"
            )

    if latitude_dim is not None:
        if lat is not None:
            raise ValueError(
                "`lat` cannot be provided together with xarray latitude-band "
                "inputs; the latitude coordinate is taken from the DataArray "
                "dimensions"
            )

        reduced_values = []
        for name, values in (("u", u), ("temperature", temperature)):
            if not isinstance(values, xr.DataArray) or latitude_dim not in values.dims:
                reduced_values.append(values)
                continue

            latitude_coord = values[latitude_dim]
            latitude_slice = (
                slice(lat_min, lat_max)
                if float(latitude_coord[0]) <= float(latitude_coord[-1])
                else slice(lat_max, lat_min)
            )
            subset = values.sel({latitude_dim: latitude_slice})
            if subset.sizes[latitude_dim] == 0:
                raise ValueError(
                    f"`lat_bounds` does not select any latitude points for `{name}`"
                )

            subset = subset.where(np.isfinite(subset))

            weights = xr.DataArray(
                np.cos(np.deg2rad(np.asarray(subset[latitude_dim], dtype=np.float64))),
                coords={latitude_dim: subset[latitude_dim]},
                dims=(latitude_dim,),
            )
            reduced_values.append(subset.weighted(weights).mean(latitude_dim))

        return reduced_values[0], reduced_values[1], None

    if (
        max(
            np.asarray(getattr(u, "data", u)).ndim,
            np.asarray(getattr(temperature, "data", temperature)).ndim,
        )
        < 3
    ):
        return u, temperature, lat

    if lat is None:
        raise ValueError(
            "NumPy latitude-band inputs with an explicit latitude axis "
            "require `lat` to provide the one-dimensional latitude "
            "coordinate"
        )

    latitude_values = np.asarray(getattr(lat, "data", lat), dtype=np.float64)
    if latitude_values.ndim != 1:
        raise ValueError(
            "`lat` must be one-dimensional when used as the latitude "
            "coordinate for NumPy latitude-band inputs"
        )

    latitude_mask = (latitude_values >= lat_min) & (latitude_values <= lat_max)
    if not np.any(latitude_mask):
        raise ValueError("`lat_bounds` does not select any latitude points")

    latitude_weights = np.cos(np.deg2rad(latitude_values[latitude_mask])).astype(
        np.float64,
        copy=False,
    )
    reduced_arrays = []
    for name, values in (("u", u), ("temperature", temperature)):
        array = np.asarray(getattr(values, "data", values), dtype=np.float64)
        if array.ndim not in {2, 3}:
            raise ValueError(
                f"`{name}` must be two- or three-dimensional when using "
                "NumPy latitude-band inputs"
            )

        pressure_axes = [
            axis for axis, size in enumerate(array.shape) if size == pressure_pa.size
        ]
        if len(pressure_axes) != 1:
            raise ValueError(
                f"`{name}` must have exactly one axis matching the "
                "pressure coordinate"
            )

        pressure_axis = pressure_axes[0]
        latitude_axes = [
            axis
            for axis, size in enumerate(array.shape)
            if axis != pressure_axis and size == latitude_values.size
        ]
        if len(latitude_axes) != 1:
            raise ValueError(
                f"`{name}` must have exactly one latitude axis matching "
                "`lat` for NumPy latitude-band inputs"
            )

        subset = np.take(array, np.flatnonzero(latitude_mask), axis=latitude_axes[0])
        subset = np.moveaxis(subset, latitude_axes[0], -1)
        subset = np.where(np.isfinite(subset), subset, np.nan)
        weight_view = latitude_weights.reshape(
            (1,) * (subset.ndim - 1) + (latitude_weights.size,)
        )
        masked_weights = np.where(np.isnan(subset), np.nan, weight_view)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            numerator = np.nanmean(subset * weight_view, axis=-1)
            denominator = np.nanmean(masked_weights, axis=-1)
        reduced_arrays.append(
            np.divide(
                numerator,
                denominator,
                out=np.full(numerator.shape, np.nan, dtype=np.float64),
                where=np.isfinite(denominator) & (denominator != 0.0),
            )
        )

    return reduced_arrays[0], reduced_arrays[1], None


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


def _normalize_smooth_window(smooth_window: int) -> int:
    """Validate the zonal-wavenumber smoothing window for growth spectra."""

    if isinstance(smooth_window, bool) or not isinstance(
        smooth_window, (int, np.integer)
    ):
        raise TypeError("`smooth_window` must be a positive odd integer")

    normalized = int(smooth_window)
    if normalized < 1:
        raise ValueError("`smooth_window` must be at least 1")
    if normalized % 2 == 0:
        raise ValueError("`smooth_window` must be an odd integer")
    return normalized


def _wrap_batch_baroc_output(
    values: np.ndarray,
    profile_coord: Optional[xr.DataArray],
) -> Union[np.ndarray, xr.DataArray]:
    """Wrap batched baroclinic growth rates back into NumPy or xarray output."""

    converted = np.array(values, dtype=np.float64, copy=True)
    if profile_coord is None:
        return converted
    return xr.DataArray(
        converted,
        dims=profile_coord.dims,
        coords={profile_coord.dims[0]: profile_coord},
        name="baroc_growth_rate",
        attrs={"units": "s^-1"},
    )


def barot_growth_rate(
    u_barotropic: ProfileInput,
    lat: ProfileInput,
    *,
    radius: float = 6_371_000.0,
    omega: float = 7.292e-5,
) -> float:
    r"""Compute the maximum barotropic normal-mode growth rate from ``U(lat)``.

    Parameters
    ----------
    u_barotropic : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        One-dimensional vertically averaged zonal-wind profile in ``m s^-1``.

    lat : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        One-dimensional latitude coordinate in degrees north. Values must be
        strictly monotonic.

    radius : float, optional
        Earth radius in meters. Only the default value is currently supported
        by the compiled kernel.

    omega : float, optional
        Planetary rotation rate in ``s^-1``. Only the default value is
        currently supported by the compiled kernel.

    Returns
    -------
    float
        Maximum barotropic growth rate in ``s^-1``. Missing or unusable
        profiles return ``NaN``.

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

    if np.any(~np.isfinite(u_values)) or np.any(~np.isfinite(lat_values)):
        return float("nan")

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

    return float(max_growth)


def baroc_growth_rate(
    u: ProfileInput,
    temperature: ProfileInput,
    pressure: ProfileInput,
    *,
    lat: Optional[LatitudeInput] = None,
    lat_bounds: Optional[Tuple[float, float]] = None,
    lon: Optional[ProfileInput] = None,
    wavenumber_mode: str = DEFAULT_WAVENUMBER_MODE,
    method: str = "log",
    tropopause_pressure: Optional[PressureInput] = None,
    solver_levels: int = DEFAULT_SOLVER_LEVELS,
    smooth_window: int = DEFAULT_SMOOTH_WINDOW,
) -> GrowthRateOutput:
    r"""Compute the maximum baroclinic normal-mode growth rate.

    Parameters
    ----------
    u : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        One-dimensional zonal-wind profile in ``m s^-1`` or a two-dimensional
        stack of profiles that shares one vertical axis with ``pressure``.

    temperature : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        One-dimensional temperature profile in Kelvin or a two-dimensional
        stack of profiles that shares one vertical axis with ``pressure``.

    pressure : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        One-dimensional pressure coordinate in Pa or hPa. Values must be
        strictly monotonic and strictly positive.

    lat : scalar, :class:`xarray.DataArray`, :class:`numpy.ndarray`, optional
        Latitude in degrees north. This may be a scalar or, for batched
        profile input, a one-dimensional vector matching the profile
        dimension. The compiled solver uses this latitude to build the local
        Coriolis parameter and meridional gradient of planetary vorticity,
        unless ``lat_bounds`` is provided.

        When ``lat_bounds`` is given and raw NumPy inputs still contain an
        explicit latitude axis, ``lat`` may instead be the one-dimensional
        latitude coordinate used to select and cosine-weight that latitude
        band before the remaining profile solve.

    lat_bounds : tuple of float, optional
        Two-element latitude range in degrees north. When provided, the
        baroclinic solver uses the Chemke-style band-mean definitions

        .. math::

           \sin_{\mathrm{avg}} = \frac{\sin(\phi_1) + \sin(\phi_2)}{2},
           \qquad
           \cos_{\mathrm{avg}} = \frac{\cos(\phi_1) + \cos(\phi_2)}{2},

        together with

        .. math::

           f = 2\Omega \sin_{\mathrm{avg}},
           \qquad
           \beta = \frac{2\Omega\cos_{\mathrm{avg}}}{a},

        to match the latitude-band formulation used in the reference
        Chemke scripts.

        If ``u`` and ``temperature`` still carry a latitude axis, this option
        also reduces that axis with cosine-latitude weighting before the
        vertical-profile solve. For xarray inputs, the latitude coordinate is
        taken from the ``lat`` or ``latitude`` dimension. For NumPy inputs
        with three dimensions, pass the one-dimensional latitude coordinate
        through ``lat=...`` so the requested latitude band can be selected
        explicitly. The recommended NumPy layout is ``(time, level, lat)``.
        More generally, the wrapper requires exactly one axis matching the
        pressure-coordinate length and one non-pressure axis matching the
        latitude-coordinate length.

    lon : :class:`xarray.DataArray`, :class:`numpy.ndarray`, optional
        One-dimensional longitude coordinate in degrees used only when
        ``wavenumber_mode="low"``. The coordinate must be strictly monotonic.
        If omitted in low-resolution mode, the wrapper uses the original
        Chemke Python-share default
        ``np.arange(0.0, 361.0, 1.5)``.

    wavenumber_mode : {"high", "low"}, optional
        Zonal-wavenumber grid used by the compiled instability solver.
        ``"high"`` uses the fixed high-resolution Chemke/MATLAB-style grid
        ``k = 0, 1\times10^{-7}, \dots, 1.99\times10^{-5}\,\mathrm{m}^{-1}``.
        ``"low"`` uses the lower-resolution Chemke Python-share formula

        .. math::

           k = \frac{2\pi w}{\left|\lambda_0 - \lambda_1\right|\pi a \cos\phi / 180},

        where :math:`w = 0, 1, \dots, \mathrm{round}(N_{\lambda}/3)-1`.
        When ``lon`` is omitted, the low-resolution mode uses the original
        Chemke Python-share longitude grid
        ``np.arange(0.0, 361.0, 1.5)``. Defaults to ``"high"``.

    method : {"log", "linear"}, optional
        Vertical interpolation method used to remap the input profiles onto the
        solver pressure grid. ``"log"`` means linear interpolation in
        log-pressure, not logarithmic interpolation of the field itself.
        Defaults to ``"log"``. The legacy alias ``"logp"`` is still accepted.

    tropopause_pressure : scalar, :class:`xarray.DataArray`, :class:`numpy.ndarray`, optional
        Explicit tropopause pressure in Pa or hPa. If provided, the automatic
        WMO tropopause diagnosis is skipped and the solver grid is built from
        this pressure to the lower-tropospheric bound. For batched profile
        input, this may be a scalar climatology or a one-dimensional vector
        matching the profile dimension.

    solver_levels : int, optional
        Number of pressure levels used when the solver grid is built
        automatically from the diagnosed or explicit tropopause pressure.
        Defaults to ``DEFAULT_SOLVER_LEVELS``.

    smooth_window : int, optional
        Optional centered running-mean window applied to the growth-rate
        spectrum over zonal wavenumber inside the compiled Fortran backend
        before the final maximum-growth diagnostic is taken. ``1`` disables
        smoothing. If given, this must be a positive odd integer.

    Returns
    -------
    float, :class:`numpy.ndarray`, :class:`xarray.DataArray`
        Maximum baroclinic growth rate in ``s^-1``. One-dimensional input
        returns a scalar float. Batched profile input returns a one-dimensional
        NumPy array or DataArray over the non-pressure dimension. Missing or
        unusable profiles return ``NaN``.

    Notes
    -----
    This function diagnoses the WMO tropopause pressure by default, or uses
    an explicit ``tropopause_pressure`` when provided, and then interpolates
    ``u`` and ``temperature`` to a fixed pressure grid between that
    tropopause and the lower troposphere. The thermodynamic profiles passed
    to the compiled kernel are

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
    If ``smooth_window`` is greater than 1, the compiled backend first applies
    a centered running mean along the discrete zonal-wavenumber spectrum and
    then evaluates the final Chemke-style turning-point maximum on that
    smoothed spectrum. When ``u`` or ``temperature`` contains multiple
    profiles, one-dimensional climatologies are broadcast across the profile
    dimension in Python, while the per-profile pressure interpolation and
    normal-mode solve loop run in the compiled Fortran backend.

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

    if lat_bounds is not None and len(lat_bounds) != 2:
        raise ValueError("`lat_bounds` must contain exactly two latitude values")

    pressure_pa = _pressure_to_pa(pressure)
    _require_monotonic(pressure_pa, "pressure")
    solver_levels = _normalize_solver_levels(solver_levels)
    smooth_window = _normalize_smooth_window(smooth_window)
    interp_method = _normalize_method(method)
    wavenumber_mode_id, wavenumber_count, lon_span_radians = _prepare_wavenumber_inputs(
        wavenumber_mode, lon
    )

    pressure_dim = pressure.dims[0] if isinstance(pressure, xr.DataArray) else None

    if lat_bounds is not None:
        u, temperature, lat = _reduce_latitude_band_inputs(
            u,
            temperature,
            lat,
            lat_bounds,
            pressure_pa,
            pressure_dim,
        )

    if lat is None and lat_bounds is None:
        raise ValueError("`lat` or `lat_bounds` is required for baroc_growth_rate")
    if lat is not None and lat_bounds is not None:
        raise ValueError("`lat` and `lat_bounds` cannot be provided together")

    u_matrix_raw = _coerce_profile_matrix(u, "u", pressure_pa.size, pressure_dim)
    temperature_matrix_raw = _coerce_profile_matrix(
        temperature,
        "temperature",
        pressure_pa.size,
        pressure_dim,
    )

    if (
        u_matrix_raw.shape[0] == 1
        and temperature_matrix_raw.shape[0] == 1
        and np.asarray(getattr(u, "data", u)).ndim == 1
        and np.asarray(getattr(temperature, "data", temperature)).ndim == 1
    ):
        u_values = u_matrix_raw[0]
        temperature_values = temperature_matrix_raw[0]
        _same_shape_or_raise(
            ("u", u_values),
            ("temperature", temperature_values),
            ("pressure", pressure_pa),
        )

        if (
            np.any(~np.isfinite(u_values))
            or np.any(~np.isfinite(temperature_values))
            or np.any(~np.isfinite(pressure_pa))
        ):
            return float("nan")

        if np.any(temperature_values <= 0.0):
            raise ValueError("temperature values must be strictly positive")

        if tropopause_pressure is None:
            tropopause_pressure_pa = _infer_tropopause_pressure_pa(
                temperature_values,
                pressure_pa,
            )
        else:
            tropopause_values = _as_profile_vector(
                tropopause_pressure,
                "tropopause_pressure",
                1,
            )
            tropopause_pressure_pa = float(tropopause_values[0])
            if np.isfinite(tropopause_pressure_pa) and tropopause_pressure_pa <= 2000.0:
                tropopause_pressure_pa *= 100.0
        bottom_pressure_pa = min(REFERENCE_PRESSURE_PA, float(np.max(pressure_pa)))
        if tropopause_pressure_pa >= bottom_pressure_pa:
            raise ValueError(
                "The diagnosed tropopause pressure must be lower than the lower-"
                "tropospheric pressure bound used by the solver grid."
            )
        solver_pressure_pa = np.linspace(
            tropopause_pressure_pa,
            bottom_pressure_pa,
            solver_levels,
            dtype=np.float64,
        )
        u_solver = np.asarray(
            interp_pressure_1d(
                u_values,
                pressure_pa,
                solver_pressure_pa,
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
                solver_pressure_pa,
                method=interp_method,
                extrapolate=False,
                missing_value=np.nan,
            ),
            dtype=np.float64,
        )

        if np.any(~np.isfinite(u_solver)) or np.any(~np.isfinite(temperature_solver)):
            return float("nan")

        theta_solver = (
            temperature_solver * (REFERENCE_PRESSURE_PA / solver_pressure_pa) ** KAPPA
        )
        if lat_bounds is None:
            latitude_value = float(_as_profile_vector(lat, "lat", 1)[0])
            latitude_radians = np.deg2rad(latitude_value)
            f_cor = float(2.0 * 7.292e-5 * np.sin(latitude_radians))
            beta = float(2.0 * 7.292e-5 * np.cos(latitude_radians) / 6_371_000.0)
            zonal_length = 1.0
            if wavenumber_mode_id == 2:
                zonal_length = float(
                    lon_span_radians * 6_371_000.0 * np.cos(latitude_radians)
                )
        else:
            f_cor, beta = _f_beta_from_lat_bounds(lat_bounds)
            zonal_length = 1.0
            if wavenumber_mode_id == 2:
                cos_avg = (
                    np.cos(np.deg2rad(float(lat_bounds[0])))
                    + np.cos(np.deg2rad(float(lat_bounds[1])))
                ) / 2.0
                zonal_length = float(lon_span_radians * 6_371_000.0 * cos_avg)
        if wavenumber_mode_id == 2 and (
            not np.isfinite(zonal_length) or zonal_length <= 0.0
        ):
            raise ValueError(
                "The low-resolution wavenumber grid requires a positive zonal "
                "length from `lon` and the selected latitude geometry."
            )
        max_growth, ier = _dbaroc_growth_rate_1d(
            u_solver,
            theta_solver,
            solver_pressure_pa,
            temperature_solver,
            f_cor,
            beta,
            smooth_window,
            wavenumber_mode_id,
            wavenumber_count,
            zonal_length,
        )
        if ier != 0:
            raise RuntimeError(
                f"baroc_growth_rate Fortran backend returned ier={ier} for the "
                "requested atmospheric column."
            )

        return float(max_growth)

    nprofile = max(u_matrix_raw.shape[0], temperature_matrix_raw.shape[0])
    profile_coord = None
    for values in (u, temperature):
        if (
            not isinstance(values, xr.DataArray)
            or values.ndim != 2
            or pressure_dim not in values.dims
        ):
            continue
        candidate_dim = (
            values.dims[0] if values.dims[1] == pressure_dim else values.dims[1]
        )
        candidate_coord = values[candidate_dim]
        if profile_coord is None:
            profile_coord = candidate_coord
        elif candidate_dim != profile_coord.dims[0]:
            raise ValueError(
                "Batched xarray inputs must share the same non-pressure dimension"
            )
    u_matrix = _broadcast_profile_matrix(
        u_matrix_raw,
        nprofile,
        "u",
    )
    temperature_matrix = _broadcast_profile_matrix(
        temperature_matrix_raw,
        nprofile,
        "temperature",
    )
    output = np.full(nprofile, np.nan, dtype=np.float64)
    valid_profile_mask = ~(
        _missing_profile_mask(u_matrix) | _missing_profile_mask(temperature_matrix)
    )

    if lat_bounds is None:
        lat_values = _as_profile_vector(lat, "lat", nprofile)
        valid_profile_mask &= ~_missing_vector_mask(lat_values)
        f_values = np.full(nprofile, np.nan, dtype=np.float64)
        beta_values = np.full(nprofile, np.nan, dtype=np.float64)
        valid_lat_mask = valid_profile_mask & np.isfinite(lat_values)
        f_values[valid_lat_mask] = (
            2.0 * 7.292e-5 * np.sin(np.deg2rad(lat_values[valid_lat_mask]))
        )
        beta_values[valid_lat_mask] = (
            2.0
            * 7.292e-5
            * np.cos(np.deg2rad(lat_values[valid_lat_mask]))
            / 6_371_000.0
        )
    else:
        f_cor, beta = _f_beta_from_lat_bounds(lat_bounds)
        f_values = np.full(nprofile, f_cor, dtype=np.float64)
        beta_values = np.full(nprofile, beta, dtype=np.float64)

    tropopause_pressure_pa = None
    if tropopause_pressure is not None:
        tropopause_values = _as_profile_vector(
            tropopause_pressure,
            "tropopause_pressure",
            nprofile,
        )
        tropopause_missing_mask = _missing_vector_mask(tropopause_values)
        tropopause_pressure_pa = np.array(
            tropopause_values, dtype=np.float64, copy=True
        )
        finite_tropopause_values = tropopause_pressure_pa[
            ~tropopause_missing_mask & np.isfinite(tropopause_pressure_pa)
        ]
        if (
            finite_tropopause_values.size > 0
            and np.nanmax(finite_tropopause_values) <= 2000.0
        ):
            tropopause_pressure_pa[~tropopause_missing_mask] *= 100.0
        valid_profile_mask &= ~tropopause_missing_mask
        finite_tropopause_mask = valid_profile_mask & np.isfinite(
            tropopause_pressure_pa
        )
        if np.any(finite_tropopause_mask & (tropopause_pressure_pa <= 0.0)):
            raise ValueError("tropopause pressure values must be strictly positive")
        bottom_pressure_pa = min(REFERENCE_PRESSURE_PA, float(np.max(pressure_pa)))
        if np.any(
            finite_tropopause_mask & (tropopause_pressure_pa >= bottom_pressure_pa)
        ):
            raise ValueError(
                "Each explicit tropopause pressure must be lower than the "
                "lower-tropospheric pressure bound used by the solver grid."
            )

    if np.any(valid_profile_mask & np.any(temperature_matrix <= 0.0, axis=1)):
        raise ValueError("temperature values must be strictly positive")

    if tropopause_pressure_pa is None:
        tropopause_result = _trop_wmo(
            temperature_matrix[valid_profile_mask],
            pressure_pa,
            xdim=-1,
            ydim=0,
            levdim=1,
            timedim=None,
            pressure_unit="Pa",
            missing_value=-999.0,
        )
        tropopause_success = np.asarray(
            tropopause_result["success"],
            dtype=bool,
        )
        tropopause_pressure_pa_valid = (
            np.asarray(tropopause_result["pressure"], dtype=np.float64) * 100.0
        )
        tropopause_pressure_pa = np.full(nprofile, np.nan, dtype=np.float64)
        valid_indices = np.flatnonzero(valid_profile_mask)
        tropopause_pressure_pa[valid_indices[tropopause_success]] = (
            tropopause_pressure_pa_valid[tropopause_success]
        )
        valid_profile_mask[valid_indices[~tropopause_success]] = False
        if not np.any(valid_profile_mask):
            return _wrap_batch_baroc_output(
                output,
                profile_coord,
            )
    bottom_pressure_pa = min(REFERENCE_PRESSURE_PA, float(np.max(pressure_pa)))
    ramp = np.linspace(0.0, 1.0, solver_levels, dtype=np.float64)
    solver_pressure_matrix_pa = (
        tropopause_pressure_pa[:, np.newaxis]
        + (bottom_pressure_pa - tropopause_pressure_pa)[:, np.newaxis]
        * ramp[np.newaxis, :]
    )

    if not np.any(valid_profile_mask):
        return _wrap_batch_baroc_output(
            output,
            profile_coord,
        )

    zonal_length_values = np.full(nprofile, 1.0, dtype=np.float64)
    if wavenumber_mode_id == 2:
        if lat_bounds is None:
            zonal_length_values[valid_profile_mask] = (
                lon_span_radians
                * 6_371_000.0
                * np.cos(np.deg2rad(lat_values[valid_profile_mask]))
            )
        else:
            cos_avg = (
                np.cos(np.deg2rad(float(lat_bounds[0])))
                + np.cos(np.deg2rad(float(lat_bounds[1])))
            ) / 2.0
            zonal_length_values[valid_profile_mask] = (
                lon_span_radians * 6_371_000.0 * cos_avg
            )
        if np.any(
            valid_profile_mask
            & (~np.isfinite(zonal_length_values) | (zonal_length_values <= 0.0))
        ):
            raise ValueError(
                "The low-resolution wavenumber grid requires a positive zonal "
                "length from `lon` and the selected latitude geometry."
            )

    interp_kind = 1 if interp_method == "linear" else 2
    growth_valid, ier_valid = _dbaroc_growth_rate_profiles(
        np.asfortranarray(u_matrix[valid_profile_mask].T),
        np.asfortranarray(temperature_matrix[valid_profile_mask].T),
        pressure_pa,
        np.asfortranarray(solver_pressure_matrix_pa[valid_profile_mask].T),
        np.asfortranarray(f_values[valid_profile_mask]),
        np.asfortranarray(beta_values[valid_profile_mask]),
        interp_kind,
        smooth_window,
        wavenumber_mode_id,
        wavenumber_count,
        np.asfortranarray(zonal_length_values[valid_profile_mask]),
    )
    growth_valid = np.asarray(growth_valid, dtype=np.float64)
    ier_valid = np.asarray(ier_valid, dtype=np.int64)

    valid_indices = np.flatnonzero(valid_profile_mask)
    output[valid_indices[ier_valid == 0]] = growth_valid[ier_valid == 0]
    if np.any((ier_valid != 0) & (ier_valid != 100)):
        failed = valid_indices[(ier_valid != 0) & (ier_valid != 100)]
        failed_codes = ier_valid[(ier_valid != 0) & (ier_valid != 100)]
        raise RuntimeError(
            "baroc_growth_rate Fortran backend returned nonzero ier values for "
            f"profile indices {failed.tolist()} with ier={failed_codes.tolist()}."
        )

    return _wrap_batch_baroc_output(
        output,
        profile_coord,
    )
