"""
Interpolation functions for hybrid-sigma and multidimensional data.

This module provides advanced interpolation capabilities for atmospheric science data,
including hybrid-sigma to pressure level interpolation and multidimensional spatial
interpolation with optional extrapolation.

References
----------
General interpolation formulations adapted for atmospheric datasets.
"""

from __future__ import annotations

import typing
import warnings

import numpy as np
import xarray as xr

from .fortran.int2p_kernels import dinterp_pressure_1d as _dinterp_pressure_1d
from .fortran.vinth2p_kernels import (
    ddelta_pressure_hybrid_pa as _ddelta_pressure_hybrid_pa,
)
from .fortran.vinth2p_kernels import (
    ddelta_pressure_hybrid_pa_into as _ddelta_pressure_hybrid_pa_into,
)
from .fortran.vinth2p_kernels import (
    dpressure_at_hybrid_levels_pa as _dpressure_at_hybrid_levels_pa,
)
from .fortran.vinth2p_kernels import (
    dpressure_at_hybrid_levels_pa_into as _dpressure_at_hybrid_levels_pa_into,
)
from .fortran.vinth2p_kernels import dsigma2hybrid_nodes as _dsigma2hybrid_nodes
from .fortran.vinth2p_kernels import (
    dsigma2hybrid_nodes_corder_into as _dsigma2hybrid_nodes_corder_into,
)
from .fortran.vinth2p_kernels import (
    dsigma2hybrid_nodes_into as _dsigma2hybrid_nodes_into,
)
from .fortran.vinth2p_kernels import (
    dvinth2p_ecmwf_nodes_corder_pa_into as _dvinth2p_ecmwf_nodes_corder_pa_into,
)
from .fortran.vinth2p_kernels import dvinth2p_ecmwf_nodes_pa as _dvinth2p_ecmwf_nodes_pa
from .fortran.vinth2p_kernels import (
    dvinth2p_ecmwf_nodes_pa_into as _dvinth2p_ecmwf_nodes_pa_into,
)
from .fortran.vinth2p_kernels import (
    dvinth2p_nodes_corder_pa_into as _dvinth2p_nodes_corder_pa_into,
)
from .fortran.vinth2p_kernels import dvinth2p_nodes_pa as _dvinth2p_nodes_pa
from .fortran.vinth2p_kernels import dvinth2p_nodes_pa_into as _dvinth2p_nodes_pa_into

__all__ = [
    "interp_pressure_1d",
    "pressure_at_hybrid_levels",
    "delta_pressure_hybrid",
    "interp_hybrid_to_pressure",
    "interp_sigma_to_hybrid",
    "interp_multidim",
]

supported_types = typing.Union[xr.DataArray, np.ndarray]
_PRESSURE_INTERP_SPVL = np.float64(np.finfo(np.float64).max)
_VINTH2P_SPVL = np.float64(-9.96921e36)

__pres_lev_mandatory__ = np.array(
    [
        1000,
        925,
        850,
        700,
        500,
        400,
        300,
        250,
        200,
        150,
        100,
        70,
        50,
        30,
        20,
        10,
        7,
        5,
        3,
        2,
        1,
    ]
).astype(
    np.float32
)  # Mandatory pressure levels (mb)
__pres_lev_mandatory__ = __pres_lev_mandatory__ * 100.0  # Convert mb to Pa


def _normalize_interp_method(method: str) -> str:
    """Normalize public interpolation method aliases to one canonical label.

    In particular, accept ``"loglog"`` from user-facing calls and normalize it
    to the internal ``"log-log"`` spelling used by the vinth2p kernels.
    """

    normalized = method.lower().replace("_", "-")
    if normalized == "loglog":
        return "log-log"
    return normalized


def _rename_colliding_coeff_dim(target, hya, hyb):
    """Avoid accidental xarray alignment when hybrid coeff dims collide."""

    if not (
        isinstance(target, xr.DataArray)
        and isinstance(hya, xr.DataArray)
        and isinstance(hyb, xr.DataArray)
    ):
        return hya, hyb

    coeff_dim = hya.dims[0]
    if coeff_dim in target.dims and target.sizes[coeff_dim] != hya.sizes[coeff_dim]:
        new_dim = "lev"
        if new_dim in target.dims:
            new_dim = "__hybrid_lev__"
        hya = hya.rename({coeff_dim: new_dim})
        hyb = hyb.rename({hyb.dims[0]: new_dim})

    return hya, hyb


def _with_dataarray_metadata(template, data, coords=None, dims=None):
    """Build a new DataArray while preserving the template metadata."""

    return xr.DataArray(
        data,
        dims=template.dims if dims is None else dims,
        coords=template.coords if coords is None else coords,
        name=template.name,
        attrs=template.attrs.copy(),
    )


def _dimension_coord_or_default(array, dim, *, output_dim=None, size=None):
    """Return a stable 1D coordinate for ``dim`` or a default integer index."""

    output_dim = dim if output_dim is None else output_dim
    coord = array.coords.get(dim)
    if coord is None:
        length = array.sizes[dim] if size is None else size
        return xr.DataArray(np.arange(length), dims=(output_dim,))

    if size is not None:
        coord = coord.isel({dim: slice(None, size)})

    return xr.DataArray(
        np.asarray(coord.data),
        dims=(output_dim,),
        attrs=coord.attrs.copy(),
    )


def _finalize_hybrid_level_output(
    result: xr.DataArray,
    *,
    lev_name: str,
    lev_coord: xr.DataArray,
    lev_dim: str | None,
    output_dims,
) -> xr.DataArray:
    """Apply optional public xarray output naming and ordering controls."""

    target_lev_dim = lev_name if lev_dim is None else lev_dim

    target_lev_coord = xr.DataArray(
        np.asarray(lev_coord.data),
        dims=(target_lev_dim,),
        attrs=lev_coord.attrs.copy(),
    )

    if target_lev_dim != lev_name:
        result = result.rename({lev_name: target_lev_dim})

    result = result.assign_coords({target_lev_dim: target_lev_coord})

    if output_dims is not None:
        requested_dims = tuple(output_dims)
        output_dims = tuple(dim for dim in requested_dims if dim in result.dims)
        output_dims += tuple(dim for dim in result.dims if dim not in output_dims)
        if len(output_dims) != result.ndim or set(output_dims) != set(result.dims):
            raise ValueError(
                "`output_dims` must contain each output dimension exactly once: "
                f"expected a permutation of {result.dims}, got {requested_dims}"
            )
        result = result.transpose(*output_dims)

    return result


def interp_pressure_1d(
    values=None,
    source_pressure=None,
    target_pressure=None,
    *,
    method: str = "log",
    extrapolate: bool = False,
    missing_value=np.nan,
    **legacy_kwargs,
):
    """Interpolate a one-dimensional profile between pressure coordinates.

    Parameters
    ----------
    values : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        One-dimensional field values defined on ``source_pressure``.

    source_pressure : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        One-dimensional source pressure levels. Values must be strictly
        monotonic after missing levels are removed.

    target_pressure : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        One-dimensional target pressure levels. Values may be increasing or
        decreasing.

    method : {"linear", "log"}, optional
        Interpolate linearly in pressure or in log-pressure. Defaults to
        ``"log"``.

    extrapolate : bool, optional
        If True, use the end slopes to extrapolate outside the valid source
        pressure range. Defaults to False.

    missing_value : scalar, optional
        Public missing-value marker. Defaults to ``np.nan``.

    Returns
    -------
    :class:`xarray.DataArray`, :class:`numpy.ndarray`
        Interpolated values on ``target_pressure``. Xarray inputs return a
        one-dimensional DataArray with the target pressure coordinate.

    Notes
    -----
    The legacy keyword aliases ``x``, ``p_in``, and ``p_out`` are still
    accepted for backward compatibility, but the preferred public parameter
    names are now ``values``, ``source_pressure``, and ``target_pressure``.
    """

    legacy_name_map = {
        "x": "values",
        "p_in": "source_pressure",
        "p_out": "target_pressure",
    }
    for legacy_name, canonical_name in legacy_name_map.items():
        if legacy_name not in legacy_kwargs:
            continue

        legacy_value = legacy_kwargs.pop(legacy_name)
        if canonical_name == "values":
            if values is not None:
                raise TypeError(
                    "interp_pressure_1d() received both `values` and legacy "
                    f"`{legacy_name}`"
                )
            values = legacy_value
        elif canonical_name == "source_pressure":
            if source_pressure is not None:
                raise TypeError(
                    "interp_pressure_1d() received both `source_pressure` and "
                    f"legacy `{legacy_name}`"
                )
            source_pressure = legacy_value
        else:
            if target_pressure is not None:
                raise TypeError(
                    "interp_pressure_1d() received both `target_pressure` and "
                    f"legacy `{legacy_name}`"
                )
            target_pressure = legacy_value

    if legacy_kwargs:
        unexpected = ", ".join(f"`{name}`" for name in sorted(legacy_kwargs))
        raise TypeError(
            "interp_pressure_1d() got unexpected keyword argument(s): " f"{unexpected}"
        )

    missing_arguments = []
    if values is None:
        missing_arguments.append("values")
    if source_pressure is None:
        missing_arguments.append("source_pressure")
    if target_pressure is None:
        missing_arguments.append("target_pressure")
    if missing_arguments:
        missing_text = ", ".join(f"`{name}`" for name in missing_arguments)
        raise TypeError(
            "interp_pressure_1d() missing required argument(s): " f"{missing_text}"
        )

    _require_compiled_interp("interp_pressure_1d", _dinterp_pressure_1d)
    _reject_lazy_or_unit_backed_inputs(
        "interp_pressure_1d",
        values,
        source_pressure,
        target_pressure,
    )

    def is_nan_missing_value(value) -> bool:
        if value is None:
            return True

        try:
            return bool(np.isnan(value))
        except TypeError:
            return False

    def pressure_interp_missing_mask(values: np.ndarray) -> np.ndarray:
        mask = ~np.isfinite(values)
        if missing_value is not None and not is_nan_missing_value(missing_value):
            mask |= values == float(missing_value)
        return mask

    def require_strict_monotonic_pressure(name: str, values: np.ndarray) -> None:
        if values.size < 2:
            return

        diffs = np.diff(values)
        if not (np.all(diffs > 0.0) or np.all(diffs < 0.0)):
            raise ValueError(
                f"{name} must be strictly monotonic after missing values are removed"
            )

    def wrap_pressure_interp_output(output_values):
        if not isinstance(values, xr.DataArray):
            return output_values

        if isinstance(target_pressure, xr.DataArray):
            output_dim = target_pressure.dims[0]
            output_coord = xr.DataArray(
                np.asarray(target_pressure.data),
                dims=(output_dim,),
                attrs=target_pressure.attrs.copy(),
            )
        else:
            output_dim = values.dims[0]
            output_coord = xr.DataArray(np.asarray(target_pressure), dims=(output_dim,))

        return xr.DataArray(
            np.asarray(output_values),
            dims=(output_dim,),
            coords={output_dim: output_coord},
            name=values.name,
            attrs=values.attrs.copy(),
        )

    normalized_method = method.lower()
    if normalized_method == "linear":
        linlog = -1 if extrapolate else 1
    elif normalized_method == "log":
        linlog = -2 if extrapolate else 2
    else:
        raise ValueError("`method` must be either 'linear' or 'log'")

    values_array = np.asarray(
        values.data if isinstance(values, xr.DataArray) else values,
        dtype=np.float64,
    )
    source_pressure_array = np.asarray(
        (
            source_pressure.data
            if isinstance(source_pressure, xr.DataArray)
            else source_pressure
        ),
        dtype=np.float64,
    )
    target_pressure_array = np.asarray(
        (
            target_pressure.data
            if isinstance(target_pressure, xr.DataArray)
            else target_pressure
        ),
        dtype=np.float64,
    )

    if (
        values_array.ndim != 1
        or source_pressure_array.ndim != 1
        or target_pressure_array.ndim != 1
    ):
        raise ValueError(
            "`values`, `source_pressure`, and `target_pressure` must each be "
            "one-dimensional"
        )

    if values_array.shape != source_pressure_array.shape:
        raise ValueError(
            "`values` and `source_pressure` must have the same length for "
            "pressure interpolation"
        )

    output_missing = (
        np.nan if is_nan_missing_value(missing_value) else float(missing_value)
    )
    result = np.full(target_pressure_array.shape, output_missing, dtype=np.float64)
    if target_pressure_array.size == 0:
        return wrap_pressure_interp_output(result)

    input_missing_mask = pressure_interp_missing_mask(
        values_array
    ) | pressure_interp_missing_mask(source_pressure_array)
    valid_source_pressure = source_pressure_array[~input_missing_mask]
    if valid_source_pressure.size < 2:
        raise ValueError(
            "interp_pressure_1d requires at least two valid input levels after "
            "missing values are removed"
        )

    require_strict_monotonic_pressure("`source_pressure`", valid_source_pressure)

    output_missing_mask = pressure_interp_missing_mask(target_pressure_array)
    valid_target_pressure = target_pressure_array[~output_missing_mask]
    require_strict_monotonic_pressure("`target_pressure`", valid_target_pressure)

    if abs(linlog) != 1:
        if np.any(valid_source_pressure <= 0.0) or np.any(valid_target_pressure <= 0.0):
            raise ValueError(
                "log-pressure interpolation requires strictly positive pressures"
            )

    if valid_target_pressure.size == 0:
        return wrap_pressure_interp_output(result)

    source_pressure_work = np.array(source_pressure_array, dtype=np.float64, copy=True)
    values_work = np.array(values_array, dtype=np.float64, copy=True)
    target_pressure_work = np.array(valid_target_pressure, dtype=np.float64, copy=True)
    source_pressure_work[pressure_interp_missing_mask(source_pressure_work)] = (
        _PRESSURE_INTERP_SPVL
    )
    values_work[pressure_interp_missing_mask(values_work)] = _PRESSURE_INTERP_SPVL

    output_valid, ier = _dinterp_pressure_1d(
        source_pressure_work,
        values_work,
        target_pressure_work,
        linlog,
        _PRESSURE_INTERP_SPVL,
    )
    output_valid = np.asarray(output_valid, dtype=np.float64)
    output_valid[output_valid == _PRESSURE_INTERP_SPVL] = output_missing

    if ier != 0:
        raise RuntimeError(
            f"interp_pressure_1d Fortran backend returned ier={ier} for the "
            "requested pressure interpolation"
        )

    result[~output_missing_mask] = output_valid
    return wrap_pressure_interp_output(result)


def pressure_at_hybrid_levels(
    psfc,
    hya,
    hyb,
    p0=100000.0,
    lev_dim: str = "lev",
    output_dims=("time", "lev", "lat", "lon"),
):
    """Calculate pressure at hybrid levels.

    Parameters
    ----------
    psfc : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        Surface pressure in Pascals.

    hya, hyb : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        One-dimensional hybrid coefficients.

    p0 : float, optional
        Reference pressure in Pascals.

    lev_dim : str, optional
        Output vertical dimension name for xarray inputs. Defaults to ``"lev"``.

    output_dims : sequence of str, optional
        Preferred output dimension order for xarray inputs. Defaults to
        ``("time", "lev", "lat", "lon")``; any names not present in the
        output are ignored, and remaining dimensions keep their relative order.

    Returns
    -------
    :class:`xarray.DataArray`, :class:`numpy.ndarray`
        Pressure at hybrid levels in Pascals. The output shape is
        ``(len(hya), *psfc.shape)`` for eager NumPy inputs, or
        ``(lev, *psfc.dims)`` for eager xarray inputs.
    """

    if not all(isinstance(x, (xr.DataArray, np.ndarray)) for x in (psfc, hya, hyb)):
        raise TypeError("psfc, hya, and hyb must be xarray DataArrays or numpy arrays")

    if not (
        all(isinstance(x, np.ndarray) for x in (psfc, hya, hyb))
        or all(isinstance(x, xr.DataArray) for x in (psfc, hya, hyb))
    ):
        raise TypeError(
            "psfc, hya, and hyb must all be the same type (all numpy arrays or all xarray DataArrays)"
        )

    if hya.shape != hyb.shape:
        raise ValueError(f"dimension mismatch: hya: {hya.shape} hyb: {hyb.shape}")

    _require_compiled_interp(
        "pressure_at_hybrid_levels",
        _dpressure_at_hybrid_levels_pa,
        _dpressure_at_hybrid_levels_pa_into,
    )

    if isinstance(hya, np.ndarray):
        if hya.ndim != 1:
            raise ValueError(
                f"hya and hyb must be 1-dimensional if numpy inputs: {hya.shape}"
            )
        default_output_dims = ("time", "lev", "lat", "lon")
        if lev_dim != "lev" or tuple(output_dims) != default_output_dims:
            raise TypeError(
                "`lev_dim` and `output_dims` are supported only for xarray inputs"
            )
        output_columns = _pressure_at_hybrid_levels_flat(psfc, hya, hyb, p0)
        return output_columns.reshape((hya.shape[0],) + psfc.shape, order="C")

    _reject_lazy_or_unit_backed_inputs("pressure_at_hybrid_levels", psfc, hya, hyb)

    if hya.dims != hyb.dims:
        warnings.warn(
            "hya and hyb have different dimension names, attempting rename",
            stacklevel=2,
        )
        hyb = hyb.rename({b: a for a, b in zip(hya.dims, hyb.dims)})

    hya, hyb = _rename_colliding_coeff_dim(psfc, hya, hyb)
    lev_name = hya.dims[0]
    lev_coord = _dimension_coord_or_default(hya, lev_name)
    output_columns = _pressure_at_hybrid_levels_flat(psfc.data, hya.data, hyb.data, p0)
    output_values = output_columns.reshape((hya.shape[0],) + psfc.shape, order="C")
    coords = {name: coord for name, coord in psfc.coords.items()}
    coords[lev_name] = lev_coord
    pressure = xr.DataArray(
        output_values,
        dims=(lev_name, *psfc.dims),
        coords=coords,
    )
    return _finalize_hybrid_level_output(
        pressure,
        lev_name=lev_name,
        lev_coord=lev_coord,
        lev_dim=lev_dim,
        output_dims=output_dims,
    )


def delta_pressure_hybrid(
    ps,
    hya,
    hyb,
    p0=100000.0,
    lev_dim: str = "lev",
    output_dims=("time", "lev", "lat", "lon"),
):
    """Calculate pressure layer thickness for hybrid coordinates.

    Parameters
    ----------
    ps : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        Surface pressure in Pascals.

    hya, hyb : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        One-dimensional hybrid coefficients defined at layer interfaces. If
        ``hya`` and ``hyb`` have length ``nlev + 1``, the result contains
        ``nlev`` layer-thickness values.

    p0 : float, optional
        Reference pressure in Pascals. Use ``p0=1`` when ``hya`` already has
        pressure units, for example ERA5 ``hyai`` / ``hyam`` coefficients.

    lev_dim : str, optional
        Output vertical dimension name for xarray inputs. Defaults to ``"lev"``.

    output_dims : sequence of str, optional
        Preferred output dimension order for xarray inputs. Defaults to
        ``("time", "lev", "lat", "lon")``; any names not present in the
        output are ignored, and remaining dimensions keep their relative order.

    Returns
    -------
    :class:`xarray.DataArray`, :class:`numpy.ndarray`
        Pressure layer thickness in Pascals for each hybrid layer. The output
        shape is ``(len(hya) - 1, *ps.shape)`` for eager NumPy inputs, or
        ``(lev, *ps.dims)`` for eager xarray inputs.

    Notes
    -----
    This function expects hybrid coefficients defined at layer interfaces, not
    layer midpoints. For datasets such as ERA5, use ``hyai`` / ``hybi`` to
    compute layer-thickness values. If you need pressure at the layer
    midpoints instead, use :func:`pressure_at_hybrid_levels` with ``hyam`` /
    ``hybm``.

    For CESM-family data it is often convenient to keep using ``hyai`` /
    ``hybi`` for the calculation while leaving the default ``lev_dim="lev"``
    and ``output_dims=("time", "lev", "lat", "lon")`` so the result lines up
    naturally with variables on model midpoints such as
    ``V(time, lev, lat, lon)``.

    With xarray inputs, the default behavior is therefore to expose the output
    vertical dimension as ``"lev"`` and to prefer the common
    ``("time", "lev", "lat", "lon")`` ordering. If one or more of these
    dimension names are not present in the output, they are ignored and the
    remaining dimensions keep their relative order.
    """

    if not all(isinstance(x, (xr.DataArray, np.ndarray)) for x in (ps, hya, hyb)):
        raise TypeError("Inputs must be xarray DataArrays or numpy arrays")

    if not isinstance(p0, (float, int, np.floating, np.integer)):
        raise TypeError(f"p0 must be a scalar numeric value, received {type(p0)}")

    if hya.shape != hyb.shape:
        raise ValueError(f"dimension mismatch: hya: {hya.shape} hyb: {hyb.shape}")

    if np.ndim(hya) != 1:
        raise ValueError(f"hya and hyb must be 1-dimensional: {hya.shape}")

    _require_compiled_interp(
        "delta_pressure_hybrid",
        _ddelta_pressure_hybrid_pa,
        _ddelta_pressure_hybrid_pa_into,
    )

    if isinstance(ps, np.ndarray):
        default_output_dims = ("time", "lev", "lat", "lon")
        if lev_dim != "lev" or tuple(output_dims) != default_output_dims:
            raise TypeError(
                "`lev_dim` and `output_dims` are supported only for xarray inputs"
            )
        hya_values = np.asarray(hya.data if isinstance(hya, xr.DataArray) else hya)
        hyb_values = np.asarray(hyb.data if isinstance(hyb, xr.DataArray) else hyb)
        output_columns = _delta_pressure_hybrid_flat(ps, hya_values, hyb_values, p0)
        return output_columns.reshape((hya_values.shape[0] - 1,) + ps.shape, order="C")

    _reject_lazy_or_unit_backed_inputs("delta_pressure_hybrid", ps, hya, hyb)

    if isinstance(hya, np.ndarray):
        hya = xr.DataArray(hya, dims=("lev",))
        hyb = xr.DataArray(hyb, dims=("lev",))
    else:
        hya = _with_dataarray_metadata(hya, np.asarray(hya.data))
        hyb = _with_dataarray_metadata(hyb, np.asarray(hyb.data))
        if hya.dims != hyb.dims:
            warnings.warn(
                "hya and hyb have different dimension names, attempting rename",
                stacklevel=2,
            )
            hyb = hyb.rename({b: a for a, b in zip(hya.dims, hyb.dims)})

    hya, hyb = _rename_colliding_coeff_dim(ps, hya, hyb)

    lev_name = hya.dims[0]
    lev_coord = _dimension_coord_or_default(hya, lev_name, size=hya.shape[0] - 1)
    output_columns = _delta_pressure_hybrid_flat(ps.data, hya.data, hyb.data, p0)
    output_values = output_columns.reshape((hya.shape[0] - 1,) + ps.shape, order="C")
    coords = {name: coord for name, coord in ps.coords.items()}
    coords[lev_name] = lev_coord
    dph = xr.DataArray(
        output_values,
        dims=(lev_name, *ps.dims),
        coords=coords,
    )
    dph.name = "dph"
    dph.attrs = {
        "long_name": "pressure layer thickness",
        "units": "Pa",
    }
    return _finalize_hybrid_level_output(
        dph,
        lev_name=lev_name,
        lev_coord=lev_coord,
        lev_dim=lev_dim,
        output_dims=output_dims,
    )


def _pressure_from_hybrid(psfc, hya, hyb, p0=100000.0):
    """Backward-compatible wrapper for :func:`pressure_at_hybrid_levels`."""

    return pressure_at_hybrid_levels(
        psfc,
        hya,
        hyb,
        p0,
        lev_dim=None,
        output_dims=None,
    )


def _pre_interp_multidim(
    data_in: xr.DataArray,
    cyclic: bool,
    missing_val,
):
    """Helper Function: Handling missing data functionality and adding cyclic
    point if required.

    Parameters
    ----------
    data_in : :class:`xarray.DataArray`
        The data on which to operate

    cyclic : :class:`bool`
        Determines if cyclic point should be added or not.
        If true then add point, else do nothing.

    missing_val : int, float, optional
        Provides an alternative to NaN

    Returns
    -------
    data_in : :class:`xarray.DataArray`
       The data input with cyclic points added (if cyclic is true)
       and missing_val values replaced with np.nan

    Notes
    -------
    """
    # replace missing_val with np.nan
    if missing_val is not None:
        data_in = _with_dataarray_metadata(
            data_in, np.where(data_in.values == missing_val, np.nan, data_in.values)
        )

    # add cyclic points and create new data array
    if cyclic:
        padded_data = np.pad(data_in.values, ((0, 0), (1, 1)), mode="wrap")
        padded_longitudes = np.pad(
            data_in.coords[data_in.dims[-1]], (1, 1), mode="wrap"
        )
        padded_longitudes[0] -= 360
        padded_longitudes[-1] += 360

        data_in = _with_dataarray_metadata(
            data_in,
            padded_data,
            coords={
                data_in.dims[-2]: data_in.coords[data_in.dims[-2]].values,
                data_in.dims[-1]: padded_longitudes,
            },
        )

    return data_in


def _post_interp_multidim(data_in, missing_val):
    """Helper Function: Handling missing data functionality.

    Parameters
    ----------
    data_in : :class:`xarray.DataArray`
        The data on which to operate

    missing_val : int, float, optional
         Provides an alternative to NaN

    Returns
    -------
    data_in : :class:`xarray.DataArray`
       The data input with np.nan values replaced with missing_val
    """
    if missing_val is not None:
        data_in = _with_dataarray_metadata(
            data_in, np.where(np.isnan(data_in.values), missing_val, data_in.values)
        )

    return data_in


def _sigma_from_hybrid(psfc, hya, hyb, p0=100000.0):
    """Calculate sigma at the hybrid levels."""

    if isinstance(hya, xr.DataArray) and isinstance(hyb, xr.DataArray):
        if hya.dims != hyb.dims:
            warnings.warn(
                "hya and hyb have different dimension names, attempting rename",
                stacklevel=2,
            )
            hyb = hyb.rename({b: a for a, b in zip(hya.dims, hyb.dims)})
        hya, hyb = _rename_colliding_coeff_dim(psfc, hya, hyb)

    # sig(k) = hya(k) * p0 / psfc + hyb(k)

    # This will be in Pa
    return hya * p0 / psfc + hyb


def _is_dask_backed(array):
    """Return True when an xarray object is backed by a dask array."""

    return (
        isinstance(array, xr.DataArray)
        and getattr(array.data, "chunks", None) is not None
    )


def _is_pint_backed(array):
    """Return True when an xarray object is backed by pint quantities."""

    if not isinstance(array, xr.DataArray):
        return False

    module_name = getattr(array.data, "__module__", "")
    return module_name.startswith("pint") or (
        hasattr(array.data, "magnitude") and hasattr(array.data, "units")
    )


def _require_compiled_interp(function_name: str, *handles) -> None:
    """Raise a clear error when a compiled interpolation backend is unavailable."""

    if not any(handle is not None for handle in handles):
        raise RuntimeError(
            f"{function_name} requires the compiled Fortran backend in "
            "`skyborn.interp.fortran.vinth2p_kernels`."
        )


def _reject_lazy_or_unit_backed_inputs(function_name: str, *arrays) -> None:
    """Reject lazy or pint-backed arrays for compiled-only interpolation APIs."""

    active_arrays = [array for array in arrays if array is not None]

    if any(_is_dask_backed(array) for array in active_arrays):
        raise NotImplementedError(
            f"{function_name} no longer provides a Dask fallback path. "
            "Use eager NumPy/xarray inputs backed by the compiled Fortran kernels."
        )

    if any(_is_pint_backed(array) for array in active_arrays):
        raise NotImplementedError(
            f"{function_name} no longer provides a Pint/MetPy fallback path. "
            "Use eager NumPy/xarray inputs backed by the compiled Fortran kernels."
        )


def _align_hybrid_level_dimension(hyam, hybm, lev_dim):
    """Rename hybrid coefficient dimensions so they align with ``lev_dim``."""

    if hyam.shape != hybm.shape:
        raise ValueError(
            f"dimension mismatch between `hyam` and `hybm`: {hyam.shape} vs {hybm.shape}"
        )

    if hyam.ndim != 1 or hybm.ndim != 1:
        raise ValueError("`hyam` and `hybm` must be one-dimensional arrays")

    if hyam.dims != hybm.dims:
        warnings.warn(
            "hyam and hybm have different dimension names, attempting rename",
            stacklevel=2,
        )
        hybm = hybm.rename(
            {source: target for target, source in zip(hyam.dims, hybm.dims)}
        )

    coeff_dim = hyam.dims[0]
    if coeff_dim != lev_dim:
        hyam = hyam.rename({coeff_dim: lev_dim})
        hybm = hybm.rename({hybm.dims[0]: lev_dim})

    return hyam, hybm


def _vinth2p_intyp(method: str) -> int:
    """Translate public interpolation methods to legacy vinth2p flags."""

    method = _normalize_interp_method(method)

    if method == "linear":
        return 1
    if method == "log":
        return 2
    if method == "log-log":
        return 3

    raise ValueError(
        f"Unknown interpolation method: {method}. "
        f'Supported methods are: "linear", "log", and "log-log".'
    )


def _vinth2p_varflg(variable: str) -> int:
    """Translate public extrapolation variable labels to ECMWF vinth flags."""

    if variable == "temperature":
        return 1
    if variable == "geopotential":
        return -1
    return 0


def _compiled_interp_output_dtype(data: xr.DataArray) -> np.dtype:
    """Return the public result dtype while the compiled kernels stay float64."""

    dtype = np.dtype(np.asarray(data.data).dtype)
    if np.issubdtype(dtype, np.floating):
        if dtype.itemsize <= np.dtype(np.float32).itemsize:
            return np.dtype(np.float32)
        return np.dtype(np.float64)
    return np.dtype(np.float64)


def _restore_compiled_interp_output_dtype(
    output_values: np.ndarray, data: xr.DataArray
) -> np.ndarray:
    """Cast compiled float64 output back to the public result precision."""

    output_dtype = _compiled_interp_output_dtype(data)
    if output_values.dtype == output_dtype:
        return output_values
    return output_values.astype(output_dtype, copy=False)


def _compiled_float64_vector(values) -> np.ndarray:
    """Materialize a 1D float64 array for the compiled interpolation boundary."""

    return np.asarray(values, dtype=np.float64).reshape(-1)


def _compiled_float64_flat(values) -> np.ndarray:
    """Materialize a contiguous 1D float64 buffer for the compiled boundary."""

    return np.ascontiguousarray(np.asarray(values, dtype=np.float64)).reshape(-1)


def _compiled_float64_columns(values, ncol: int, nlev: int) -> np.ndarray:
    """Materialize a Fortran-ordered float64 (nlev, ncol) matrix for kernels."""

    return np.asfortranarray(np.asarray(values, dtype=np.float64).reshape(ncol, nlev).T)


def _compiled_float64_output(shape, order: str = "F") -> np.ndarray:
    """Allocate a float64 output buffer for compiled interpolation kernels."""

    return np.empty(shape, dtype=np.float64, order=order)


def _delta_pressure_hybrid_flat(ps_values, hya_values, hyb_values, p0) -> np.ndarray:
    """Run eager hybrid layer-thickness calculation through the compiled backend."""

    if _ddelta_pressure_hybrid_pa is None:
        raise RuntimeError("delta-pressure backend is not available")

    ps_flat = _compiled_float64_flat(ps_values)
    hya_vector = _compiled_float64_vector(hya_values)
    hyb_vector = _compiled_float64_vector(hyb_values)
    ncol = ps_flat.size
    nlev = hya_vector.size
    nlevo = nlev - 1
    output_columns = _compiled_float64_output((nlevo, ncol), order="F")

    if _ddelta_pressure_hybrid_pa_into is not None:
        try:
            _ddelta_pressure_hybrid_pa_into(
                ps_flat,
                output_columns,
                hya_vector,
                hyb_vector,
                float(p0),
                nlevo=nlevo,
                ncol=ncol,
                nlev=nlev,
            )
        except Exception as exc:
            if "ddelta_pressure_hybrid_pa_into" not in str(exc):
                raise
            _ddelta_pressure_hybrid_pa_into(
                ps_flat,
                output_columns,
                hya_vector,
                hyb_vector,
                float(p0),
                nlevo,
            )
        return output_columns

    try:
        return _ddelta_pressure_hybrid_pa(
            ps_flat,
            hya_vector,
            hyb_vector,
            float(p0),
            nlevo=nlevo,
            ncol=ncol,
            nlev=nlev,
        )
    except Exception as exc:
        if "ddelta_pressure_hybrid_pa" not in str(exc):
            raise
        return _ddelta_pressure_hybrid_pa(
            ps_flat,
            hya_vector,
            hyb_vector,
            float(p0),
            nlevo,
        )


def _pressure_at_hybrid_levels_flat(
    ps_values, hya_values, hyb_values, p0
) -> np.ndarray:
    """Run eager hybrid-pressure calculation through the compiled backend."""

    if _dpressure_at_hybrid_levels_pa is None:
        raise RuntimeError("pressure-at-hybrid-levels backend is not available")

    ps_flat = _compiled_float64_flat(ps_values)
    hya_vector = _compiled_float64_vector(hya_values)
    hyb_vector = _compiled_float64_vector(hyb_values)
    ncol = ps_flat.size
    nlev = hya_vector.size
    output_columns = _compiled_float64_output((nlev, ncol), order="F")

    if _dpressure_at_hybrid_levels_pa_into is not None:
        try:
            _dpressure_at_hybrid_levels_pa_into(
                ps_flat,
                output_columns,
                hya_vector,
                hyb_vector,
                float(p0),
                ncol=ncol,
                nlev=nlev,
            )
        except Exception as exc:
            if "dpressure_at_hybrid_levels_pa_into" not in str(exc):
                raise
            _dpressure_at_hybrid_levels_pa_into(
                ps_flat,
                output_columns,
                hya_vector,
                hyb_vector,
                float(p0),
            )
        return output_columns

    try:
        return _dpressure_at_hybrid_levels_pa(
            ps_flat,
            hya_vector,
            hyb_vector,
            float(p0),
            ncol=ncol,
            nlev=nlev,
        )
    except Exception as exc:
        if "dpressure_at_hybrid_levels_pa" not in str(exc):
            raise
        return _dpressure_at_hybrid_levels_pa(
            ps_flat,
            hya_vector,
            hyb_vector,
            float(p0),
        )


def _as_c_contiguous_compiled_flat(array: xr.DataArray, dims, shape):
    """Return a float64 flat buffer for aligned C-order floating inputs."""

    if array is None or not isinstance(array, xr.DataArray):
        return None

    if tuple(array.dims) != tuple(dims) or tuple(array.shape) != tuple(shape):
        return None

    values = array.data
    if not isinstance(values, np.ndarray):
        return None
    if not values.flags.c_contiguous or not np.issubdtype(values.dtype, np.floating):
        return None

    if values.dtype == np.float64:
        return values.reshape(-1)

    return _compiled_float64_flat(values)


def _as_broadcast_float64_flat(array: xr.DataArray, template: xr.DataArray):
    """Return a flattened float64 array after lightweight xarray broadcasting."""

    if array is None or not isinstance(array, xr.DataArray):
        return None

    try:
        broadcast = array.broadcast_like(template)
    except Exception:
        return None

    return _compiled_float64_flat(broadcast.data)


def _build_vinth2p_output(data, interp_axis, new_levels, output_values, base_template):
    """Wrap a NumPy output array in the public xarray result."""

    output_values = _restore_compiled_interp_output_dtype(output_values, data)
    coords = {k: v for k, v in base_template.coords.items()}
    coords["plev"] = new_levels
    output_dims = tuple(
        dim if idx != interp_axis else "plev" for idx, dim in enumerate(data.dims)
    )
    return xr.DataArray(
        output_values,
        dims=output_dims,
        coords=coords,
        name=data.name,
        attrs=data.attrs,
    )


def _interp_hybrid_to_pressure_corder(
    data: xr.DataArray,
    ps: xr.DataArray,
    hyam: xr.DataArray,
    hybm: xr.DataArray,
    p0: float,
    new_levels: np.ndarray,
    lev_dim: str,
    method: str,
    extrapolate: bool,
    variable: str,
    t_bot: xr.DataArray,
    phi_sfc: xr.DataArray,
) -> xr.DataArray | None:
    """Use a no-transpose fast path for NumPy C-order arrays with aligned dims."""

    if _dvinth2p_nodes_corder_pa_into is None:
        return None

    interp_axis = data.dims.index(lev_dim)
    raw_data = data.data
    if not isinstance(raw_data, np.ndarray):
        return None
    if not raw_data.flags.c_contiguous or not np.issubdtype(
        raw_data.dtype, np.floating
    ):
        return None

    shape_before = data.shape[:interp_axis]
    shape_after = data.shape[interp_axis + 1 :]
    lead_dims = data.dims[:interp_axis] + data.dims[interp_axis + 1 :]
    lead_shape = shape_before + shape_after
    nouter = int(np.prod(shape_before, dtype=np.int64)) if shape_before else 1
    ninner = int(np.prod(shape_after, dtype=np.int64)) if shape_after else 1
    nlevi = data.shape[interp_axis]
    nlevo = new_levels.size
    raw_data_flat = (
        raw_data.reshape(-1)
        if raw_data.dtype == np.float64
        else _compiled_float64_flat(raw_data)
    )

    ps_flat = _as_c_contiguous_compiled_flat(ps, lead_dims, lead_shape)
    if ps_flat is None or not np.isfinite(ps_flat).all():
        return None

    hyam_values = _compiled_float64_vector(hyam.data)
    hybm_values = _compiled_float64_vector(hybm.data)
    new_level_values = _compiled_float64_vector(new_levels)
    intyp = _vinth2p_intyp(method)
    base_template = data.isel({lev_dim: 0}, drop=True)

    output_values = _compiled_float64_output(
        (*shape_before, nlevo, *shape_after), order="C"
    )
    output_flat = output_values.reshape(-1)

    if extrapolate:
        varflg = _vinth2p_varflg(variable)
        if _dvinth2p_ecmwf_nodes_corder_pa_into is None:
            return None
        if varflg == 0:
            tbot_flat = np.zeros(nouter * ninner, dtype=np.float64)
            phi_flat = np.zeros(nouter * ninner, dtype=np.float64)
        else:
            tbot_flat = _as_c_contiguous_compiled_flat(t_bot, lead_dims, lead_shape)
            if tbot_flat is None:
                tbot_flat = _as_broadcast_float64_flat(t_bot, base_template)
            phi_flat = _as_c_contiguous_compiled_flat(phi_sfc, lead_dims, lead_shape)
            if phi_flat is None:
                phi_flat = _as_broadcast_float64_flat(phi_sfc, base_template)
            if tbot_flat is None or phi_flat is None:
                return None
        _dvinth2p_ecmwf_nodes_corder_pa_into(
            raw_data_flat,
            output_flat,
            hyam_values,
            hybm_values,
            float(p0),
            new_level_values,
            intyp,
            ps_flat,
            _VINTH2P_SPVL,
            1,
            nouter,
            ninner,
            varflg,
            tbot_flat,
            phi_flat,
        )
    else:
        _dvinth2p_nodes_corder_pa_into(
            raw_data_flat,
            output_flat,
            hyam_values,
            hybm_values,
            float(p0),
            new_level_values,
            intyp,
            ps_flat,
            _VINTH2P_SPVL,
            0,
            nouter,
            ninner,
        )

    output_values[output_values == _VINTH2P_SPVL] = np.nan
    return _build_vinth2p_output(
        data=data,
        interp_axis=interp_axis,
        new_levels=new_levels,
        output_values=output_values,
        base_template=base_template,
    )


def _interp_hybrid_to_pressure(
    data: xr.DataArray,
    ps: xr.DataArray,
    hyam: xr.DataArray,
    hybm: xr.DataArray,
    p0: float,
    new_levels: np.ndarray,
    lev_dim: str,
    method: str,
    extrapolate: bool,
    variable: str,
    t_bot: xr.DataArray,
    phi_sfc: xr.DataArray,
) -> xr.DataArray:
    """Run eager hybrid-to-pressure interpolation through the compiled vinth2p backend."""

    if _dvinth2p_nodes_pa is None:
        raise RuntimeError("vinth2p backend is not available")

    corder_output = _interp_hybrid_to_pressure_corder(
        data=data,
        ps=ps,
        hyam=hyam,
        hybm=hybm,
        p0=p0,
        new_levels=new_levels,
        lev_dim=lev_dim,
        method=method,
        extrapolate=extrapolate,
        variable=variable,
        t_bot=t_bot,
        phi_sfc=phi_sfc,
    )
    if corder_output is not None:
        return corder_output

    base_template = data.isel({lev_dim: 0}, drop=True)
    lead_dims = base_template.dims
    interp_axis = data.dims.index(lev_dim)
    data_view = data.transpose(*lead_dims, lev_dim)

    nlevi = data_view.shape[-1]
    nlevo = new_levels.size
    lead_shape = data_view.shape[:-1]
    ncol = int(np.prod(lead_shape, dtype=np.int64)) if lead_shape else 1

    data_columns = _compiled_float64_columns(data_view.data, ncol=ncol, nlev=nlevi)
    ps_columns = _compiled_float64_flat(
        ps.broadcast_like(base_template).transpose(*lead_dims).data
    )
    ps_columns = np.where(np.isfinite(ps_columns), ps_columns, _VINTH2P_SPVL)

    hyam_values = _compiled_float64_vector(hyam.data)
    hybm_values = _compiled_float64_vector(hybm.data)
    new_level_values = _compiled_float64_vector(new_levels)
    intyp = _vinth2p_intyp(method)
    output_columns = _compiled_float64_output((nlevo, ncol), order="F")

    if extrapolate:
        varflg = _vinth2p_varflg(variable)
        if varflg == 0:
            t_bot_columns = np.zeros(ncol, dtype=np.float64)
            phi_columns = np.zeros(ncol, dtype=np.float64)
        else:
            t_bot_columns = _compiled_float64_flat(
                t_bot.broadcast_like(base_template).transpose(*lead_dims).data
            )
            phi_columns = _compiled_float64_flat(
                phi_sfc.broadcast_like(base_template).transpose(*lead_dims).data
            )
        if _dvinth2p_ecmwf_nodes_pa_into is not None:
            _dvinth2p_ecmwf_nodes_pa_into(
                data_columns,
                output_columns,
                hyam_values,
                hybm_values,
                float(p0),
                new_level_values,
                intyp,
                ps_columns,
                _VINTH2P_SPVL,
                1,
                varflg,
                t_bot_columns,
                phi_columns,
            )
        else:
            output_columns = _dvinth2p_ecmwf_nodes_pa(
                data_columns,
                hyam_values,
                hybm_values,
                float(p0),
                new_level_values,
                intyp,
                ps_columns,
                _VINTH2P_SPVL,
                1,
                varflg,
                t_bot_columns,
                phi_columns,
            )
    else:
        if _dvinth2p_nodes_pa_into is not None:
            _dvinth2p_nodes_pa_into(
                data_columns,
                output_columns,
                hyam_values,
                hybm_values,
                float(p0),
                new_level_values,
                intyp,
                ps_columns,
                _VINTH2P_SPVL,
                0,
            )
        else:
            output_columns = _dvinth2p_nodes_pa(
                data_columns,
                hyam_values,
                hybm_values,
                float(p0),
                new_level_values,
                intyp,
                ps_columns,
                _VINTH2P_SPVL,
                0,
            )

    output_values = output_columns.T.reshape((*lead_shape, nlevo))
    output_values[output_values == _VINTH2P_SPVL] = np.nan
    output_values = _restore_compiled_interp_output_dtype(output_values, data)

    coords = {k: v for k, v in base_template.coords.items()}
    coords["plev"] = new_levels
    output = xr.DataArray(
        output_values,
        dims=(*lead_dims, "plev"),
        coords=coords,
        name=data.name,
        attrs=data.attrs,
    )

    dims = [data.dims[i] if i != interp_axis else "plev" for i in range(data.ndim)]
    return output.transpose(*dims)


def _interp_sigma_to_hybrid(
    data: xr.DataArray,
    sig_coords: xr.DataArray,
    ps: xr.DataArray,
    hyam: xr.DataArray,
    hybm: xr.DataArray,
    p0: float,
    lev_dim: str,
    method: str,
) -> xr.DataArray:
    """Run eager sigma-to-hybrid interpolation through the compiled backend."""

    if _dsigma2hybrid_nodes is None:
        raise RuntimeError("sigma2hybrid backend is not available")

    sigma_source_values = _compiled_float64_vector(
        sig_coords.data if isinstance(sig_coords, xr.DataArray) else sig_coords
    )
    sigma_diffs = np.diff(sigma_source_values)
    if not (
        np.all(sigma_diffs >= 0.0)
        or np.all(sigma_diffs <= 0.0)
        or sigma_source_values.size <= 1
    ):
        raise ValueError("sigma2hybrid backend requires monotonic sigma coordinates")

    corder_output = _interp_sigma_to_hybrid_corder(
        data=data,
        ps=ps,
        hyam=hyam,
        hybm=hybm,
        p0=p0,
        lev_dim=lev_dim,
        method=method,
        sigma_source_values=sigma_source_values,
    )
    if corder_output is not None:
        return corder_output

    sigma_target = _sigma_from_hybrid(ps, hyam, hybm, p0)
    target_dim = sigma_target.dims[0]
    base_template = data.isel({lev_dim: 0}, drop=True)
    lead_dims = base_template.dims
    interp_axis = data.dims.index(lev_dim)

    data_view = data.transpose(*lead_dims, lev_dim)
    sigma_target_view = sigma_target.transpose(target_dim, *lead_dims)

    nlevi = data_view.shape[-1]
    nlevo = sigma_target_view.shape[0]
    lead_shape = data_view.shape[:-1]
    ncol = int(np.prod(lead_shape, dtype=np.int64)) if lead_shape else 1

    data_columns = _compiled_float64_columns(data_view.data, ncol=ncol, nlev=nlevi)
    sigma_target_columns = np.asfortranarray(
        np.asarray(sigma_target_view.data, dtype=np.float64).reshape(nlevo, ncol)
    )

    intyp = _vinth2p_intyp(method)
    output_columns = _compiled_float64_output((nlevo, ncol), order="F")

    if _dsigma2hybrid_nodes_into is not None:
        _dsigma2hybrid_nodes_into(
            data_columns,
            output_columns,
            sigma_source_values,
            sigma_target_columns,
            intyp,
            _VINTH2P_SPVL,
            nlevo,
            nlevi=nlevi,
            ncol=ncol,
        )
    else:
        output_columns = _dsigma2hybrid_nodes(
            data_columns,
            sigma_source_values,
            sigma_target_columns,
            intyp,
            _VINTH2P_SPVL,
            nlevo,
            nlevi=nlevi,
            ncol=ncol,
        )

    output_values = output_columns.T.reshape((*lead_shape, nlevo))
    output_values[output_values == _VINTH2P_SPVL] = np.nan
    output_values = _restore_compiled_interp_output_dtype(output_values, data)

    coords = {k: v for k, v in base_template.coords.items()}
    coords["hlev"] = _dimension_coord_or_default(hyam, hyam.dims[0], output_dim="hlev")
    output = xr.DataArray(
        output_values,
        dims=(*lead_dims, "hlev"),
        coords=coords,
        name=data.name,
        attrs=data.attrs,
    )

    dims = [data.dims[i] if i != interp_axis else "hlev" for i in range(data.ndim)]
    return output.transpose(*dims)


def _interp_sigma_to_hybrid_corder(
    data: xr.DataArray,
    ps: xr.DataArray,
    hyam: xr.DataArray,
    hybm: xr.DataArray,
    p0: float,
    lev_dim: str,
    method: str,
    sigma_source_values: np.ndarray,
) -> xr.DataArray | None:
    """Use a no-transpose fast path for NumPy C-order sigma->hybrid inputs."""

    if _dsigma2hybrid_nodes_corder_into is None:
        return None

    interp_axis = data.dims.index(lev_dim)
    raw_data = data.data
    if not isinstance(raw_data, np.ndarray):
        return None
    if not raw_data.flags.c_contiguous or not np.issubdtype(
        raw_data.dtype, np.floating
    ):
        return None

    shape_before = data.shape[:interp_axis]
    shape_after = data.shape[interp_axis + 1 :]
    lead_dims = data.dims[:interp_axis] + data.dims[interp_axis + 1 :]
    lead_shape = shape_before + shape_after
    nouter = int(np.prod(shape_before, dtype=np.int64)) if shape_before else 1
    ninner = int(np.prod(shape_after, dtype=np.int64)) if shape_after else 1
    nlevi = data.shape[interp_axis]
    nlevo = hyam.shape[0]
    base_template = data.isel({lev_dim: 0}, drop=True)
    raw_data_flat = (
        raw_data.reshape(-1)
        if raw_data.dtype == np.float64
        else _compiled_float64_flat(raw_data)
    )

    ps_flat = _as_c_contiguous_compiled_flat(ps, lead_dims, lead_shape)
    if ps_flat is None:
        ps_flat = _as_broadcast_float64_flat(ps, base_template)
    if ps_flat is None or not np.isfinite(ps_flat).all():
        return None

    hyam_values = _compiled_float64_vector(hyam.data)
    hybm_values = _compiled_float64_vector(hybm.data)
    intyp = _vinth2p_intyp(method)
    output_values = _compiled_float64_output(
        (*shape_before, nlevo, *shape_after), order="C"
    )
    output_flat = output_values.reshape(-1)

    _dsigma2hybrid_nodes_corder_into(
        raw_data_flat,
        output_flat,
        sigma_source_values,
        hyam_values,
        hybm_values,
        float(p0),
        ps_flat,
        intyp,
        _VINTH2P_SPVL,
        nouter,
        ninner,
        nlevo,
        nlevi=nlevi,
    )

    output_values[output_values == _VINTH2P_SPVL] = np.nan
    output_values = _restore_compiled_interp_output_dtype(output_values, data)
    coords = {k: v for k, v in base_template.coords.items()}
    coords["hlev"] = _dimension_coord_or_default(hyam, hyam.dims[0], output_dim="hlev")
    output = xr.DataArray(
        output_values,
        dims=tuple(
            dim if idx != interp_axis else "hlev" for idx, dim in enumerate(data.dims)
        ),
        coords=coords,
        name=data.name,
        attrs=data.attrs,
    )
    return output


def interp_hybrid_to_pressure(
    data: xr.DataArray,
    ps: xr.DataArray,
    hyam: xr.DataArray,
    hybm: xr.DataArray,
    p0: float = 100000.0,
    new_levels: np.ndarray = __pres_lev_mandatory__,
    lev_dim: str = None,
    method: str = "linear",
    extrapolate: bool = False,
    variable: str = None,
    t_bot: xr.DataArray = None,
    phi_sfc: xr.DataArray = None,
) -> xr.DataArray:
    """Interpolate and extrapolate data from hybrid-sigma levels to isobaric levels.

    This function interpolates atmospheric data from hybrid-sigma coordinate levels
    to constant pressure levels, with optional extrapolation below ground using
    ECMWF formulations. Preserves all metadata from the input data.

    Notes
    -----
    Atmosphere hybrid-sigma pressure coordinates are commonly defined in two different
    ways as described below and in CF Conventions. This particular function expects the
    first formulation. However, with some minor adjustments on the user side it can
    support datasets leveraging the second formulation as well. In this case, you can
    set the input parameters p0=1 and hyam=ap to adapt the function to meet your needs.

    Formulation 1: p(n,k,j,i) = a(k)*p0 + b(k)*ps(n,j,i)
    Formulation 2: p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)

    Parameters
    ----------
    data : :class:`xarray.DataArray`
        Multidimensional data array with hybrid-sigma levels and a ``lev_dim`` coordinate.

    ps : :class:`xarray.DataArray`
        A multi-dimensional array of surface pressures (Pa), same time/space shape as data.

    hyam, hybm : :class:`xarray.DataArray`
        One-dimensional arrays containing the hybrid A and B coefficients. Must have the same
        dimension size as the ``lev_dim`` dimension of data.

    p0 : float, optional
        Scalar numeric value equal to surface reference pressure (Pa). Defaults to 100000 Pa.

    new_levels : ndarray, optional
        A one-dimensional array of output pressure levels (Pa). If not given, the mandatory
        list of 21 pressure levels is used.

    lev_dim : str, optional
        String that is the name of level dimension in data. If None, attempts to detect
        automatically using CF conventions.

    method : str, optional
        Interpolation method passed to the compiled vinth2p kernels.
        Supported values are ``"linear"``, ``"log"``, and ``"log-log"``
        (or ``"loglog"``). Defaults to ``"linear"``.

    extrapolate : bool, optional
        If True, below ground extrapolation for ``variable`` will be done using
        an `ECMWF formulation <https://dx.doi.org/10.5065/D6HX19NH>`__. Defaults
        to False.

    variable : str, optional
        String representing what variable is extrapolated below surface level.
        Temperature extrapolation = "temperature". Geopotential height
        extrapolation = "geopotential". All other variables = "other". If
        "other", the value of ``data`` at the lowest model level will be used
        as the below ground fill value. Required if extrapolate is True.

    t_bot : :class:`xarray.DataArray`, optional
        Temperature in Kelvin at the lowest layer of the model. Not necessarily
        the same as surface temperature. Required if ``extrapolate`` is True
        and ``variable`` is not ``'other'``

    phi_sfc: :class:`xarray.DataArray`, optional
        Geopotential in J/kg at the lowest layer of the model. Not necessarily
        the same as surface geopotential. Required if ``extrapolate`` is True
        and ``variable`` is not ``'other'``.

    Returns
    -------
    output : :class:`xarray.DataArray`
        Interpolated data with isobaric levels as the new vertical coordinate

    Examples
    --------
    Basic interpolation from hybrid-sigma to pressure levels:

    >>> import skyborn.interp as si
    >>> import xarray as xr
    >>> import numpy as np
    >>>
    >>> # Interpolate temperature to standard pressure levels
    >>> temp_p = si.interp_hybrid_to_pressure(
    ...     data=temperature,
    ...     ps=surface_pressure,
    ...     hyam=hybrid_a,
    ...     hybm=hybrid_b,
    ...     new_levels=np.array([100000, 85000, 70000])  # Pa
    ... )

    See Also
    --------
    interp_sigma_to_hybrid : Interpolate from sigma to hybrid coordinates

    Related NCL Functions:
    `vinth2p <https://www.ncl.ucar.edu/Document/Functions/Built-in/vinth2p.shtml>`__,
    `vinth2p_ecmwf <https://www.ncl.ucar.edu/Document/Functions/Built-in/vinth2p_ecmwf.shtml>`__
    """

    if not all(isinstance(x, xr.DataArray) for x in (data, ps, hyam, hybm)):
        raise TypeError("data, ps, hyam, and hybm must be xarray DataArray objects")

    _require_compiled_interp(
        "interp_hybrid_to_pressure",
        _dvinth2p_nodes_pa,
        _dvinth2p_nodes_pa_into,
        _dvinth2p_nodes_corder_pa_into,
        _dvinth2p_ecmwf_nodes_pa,
        _dvinth2p_ecmwf_nodes_pa_into,
        _dvinth2p_ecmwf_nodes_corder_pa_into,
    )
    _reject_lazy_or_unit_backed_inputs(
        "interp_hybrid_to_pressure", data, ps, hyam, hybm, t_bot, phi_sfc
    )

    new_levels = np.asarray(new_levels)

    # Check inputs
    if extrapolate and (variable is None):
        raise ValueError("If `extrapolate` is True, `variable` must be provided.")

    if variable in ["geopotential", "temperature"] and (
        t_bot is None or phi_sfc is None
    ):
        raise ValueError(
            "If `variable` is 'geopotential' or 'temperature', both `t_bot` and `phi_sfc` must be provided"
        )

    if variable not in ["geopotential", "temperature", "other", None]:
        raise ValueError(
            "The value of `variable` is "
            + variable
            + ", but the accepted values are 'temperature', 'geopotential', 'other', or None."
        )

    # Determine the level dimension and then the interpolation axis
    if lev_dim is None:
        try:
            if hasattr(data.cf, "guess_coord_axis"):
                data = data.cf.guess_coord_axis()
            lev_dim = data.cf["vertical"].name
        except Exception as exc:
            raise ValueError(
                "Unable to determine vertical dimension name. Please specify the name via `lev_dim` argument."
            ) from exc

    hyam, hybm = _align_hybrid_level_dimension(hyam, hybm, lev_dim)
    method = _normalize_interp_method(method)
    return _interp_hybrid_to_pressure(
        data=data,
        ps=ps,
        hyam=hyam,
        hybm=hybm,
        p0=p0,
        new_levels=new_levels,
        lev_dim=lev_dim,
        method=method,
        extrapolate=extrapolate,
        variable=variable,
        t_bot=t_bot,
        phi_sfc=phi_sfc,
    )


def interp_sigma_to_hybrid(
    data: xr.DataArray,
    sig_coords: xr.DataArray,
    ps: xr.DataArray,
    hyam: xr.DataArray,
    hybm: xr.DataArray,
    p0: float = 100000.0,
    lev_dim: str = None,
    method: str = "linear",
) -> xr.DataArray:
    """Interpolate data from sigma to hybrid coordinates.

    This function interpolates atmospheric data from sigma coordinate levels
    to hybrid-sigma coordinate levels. Preserves all metadata from the input data.

    Parameters
    ----------
    data : :class:`xarray.DataArray`
        Multidimensional data array with sigma levels and a ``lev_dim`` coordinate.

    sig_coords : :class:`xarray.DataArray`
        A one-dimensional array of sigma coordinates corresponding to the ``lev_dim``
        dimension of ``data``.

    ps : :class:`xarray.DataArray`
        A multi-dimensional array of surface pressures (Pa), same time/space shape as data.

    hyam, hybm : :class:`xarray.DataArray`
        One-dimensional arrays containing the hybrid A and B coefficients. Must have the same
        dimension as the desired output hybrid levels.

    p0 : float, optional
        Scalar numeric value equal to surface reference pressure (Pa). Defaults to 100000 Pa.

    lev_dim : str, optional
        String that is the name of level dimension in data. If None, attempts to detect
        automatically using CF conventions.

    method : str, optional
        Interpolation method passed to the compiled sigma2hybrid kernels.
        Supported values are ``"linear"``, ``"log"``, and ``"log-log"``
        (or ``"loglog"``). Defaults to ``"linear"``.

    Returns
    -------
    output : :class:`xarray.DataArray`
        Interpolated data with hybrid levels as the new vertical coordinate

    Examples
    --------
    Basic interpolation from sigma to hybrid coordinates:

    >>> import skyborn.interp as si
    >>> import xarray as xr
    >>> import numpy as np
    >>>
    >>> # Interpolate data from sigma to hybrid levels
    >>> data_hybrid = si.interp_sigma_to_hybrid(
    ...     data=sigma_data,
    ...     sig_coords=sigma_levels,
    ...     ps=surface_pressure,
    ...     hyam=hybrid_a,
    ...     hybm=hybrid_b
    ... )

    See Also
    --------
    interp_hybrid_to_pressure : Interpolate from hybrid to pressure coordinates

    Related NCL Function:
    `sigma2hybrid <https://www.ncl.ucar.edu/Document/Functions/Built-in/sigma2hybrid.shtml>`__
    """

    _require_compiled_interp(
        "interp_sigma_to_hybrid",
        _dsigma2hybrid_nodes,
        _dsigma2hybrid_nodes_into,
        _dsigma2hybrid_nodes_corder_into,
    )
    _reject_lazy_or_unit_backed_inputs(
        "interp_sigma_to_hybrid", data, sig_coords, ps, hyam, hybm
    )

    # Determine the level dimension and then the interpolation axis
    if lev_dim is None:
        try:
            lev_dim = data.cf["vertical"].name
        except Exception:
            raise ValueError(
                "Unable to determine vertical dimension name. Please specify the name via `lev_dim` argument.'"
            )

    method = _normalize_interp_method(method)
    return _interp_sigma_to_hybrid(
        data=data,
        sig_coords=sig_coords,
        ps=ps,
        hyam=hyam,
        hybm=hybm,
        p0=p0,
        lev_dim=lev_dim,
        method=method,
    )


def interp_multidim(
    data_in: supported_types,
    lat_out: np.ndarray,
    lon_out: np.ndarray,
    lat_in: np.ndarray = None,
    lon_in: np.ndarray = None,
    cyclic: bool = False,
    missing_val: np.number = None,
    method: str = "linear",
    fill_value: typing.Union[str, np.number] = np.nan,
) -> supported_types:
    """Multidimensional interpolation of variables. Uses ``xarray.interp`` to
    perform interpolation. Will not perform extrapolation by default, returns
    missing values if any surrounding points contain missing values.

    .. warning::
        The output data type may be promoted to that of the coordinate data.

    Parameters
    ----------
    data_in : :class:`xarray.DataArray`, ndarray
        Data array with data to be interpolated and associated coords. If
        it is a np array, then ``lat_in`` and ``lon_in`` must be provided. Length must
        be coordinated with given coordinates.

    lat_out: ndarray
        List of latitude coordinates to be interpolated to.

    lon_out: ndarray
        List of longitude coordinates to be interpolated to.

    lat_in: ndarray
        List of latitude coordinates corresponding to ``data_in``. Must be
        given if ``data_in`` is not an xarray.

    lon_in: ndarray
        List of longitude coordinates corresponding to ``data_in``. Must be
        given if ``data_in`` is not an xarray.

    cyclic: bool, optional
        Set as true if lon values are cyclical but do not fully wrap around
        the globe
        (0, 1.5, 3, ..., 354, 355.5) Default is false

    missing_val : :class:`np.number`, optional
        Provide a number to represent missing data. Alternative to using ``np.nan``

    method: str, optional
        Provide specific method of interpolation. Default is "linear"
        “linear” or “nearest” for multidimensional array

    fill_value: str, optional
        Set as 'extrapolate' to allow extrapolation of data. Default is
        no extrapolation.

    Returns
    -------
    data_out : ndarray, :class:`xarray.DataArray`
       Returns the same type of object as input ``data_in``. However, the type of
       the data in the array may be promoted to that of the coordinates. Shape
       will be the same as input array except for last two dimensions which will
       be equal to the coordinates given in ``data_out``.

    Examples
    --------
    Basic multidimensional interpolation:

    >>> import skyborn.interp as si
    >>> import xarray as xr
    >>> import numpy as np
    >>>
    >>> # Create sample data
    >>> data = np.asarray([[1, 2, 3, 4, 5, 99],
    ...                   [2, 4, 6, 8, 10, 12]])
    >>> lat_in = [0, 1]
    >>> lon_in = [0, 50, 100, 250, 300, 350]
    >>> data_in = xr.DataArray(data,
    ...                        dims=['lat', 'lon'],
    ...                        coords={'lat':lat_in,
    ...                                'lon': lon_in})
    >>>
    >>> # Interpolate to new coordinates with cyclic boundary
    >>> do = si.interp_multidim(data_in,
    ...                         [0, 1],
    ...                         [0, 50, 360],
    ...                         cyclic=True,
    ...                         missing_val=99)
    >>> print(do)
    <xarray.DataArray (lat: 2, lon: 3)>
    array([[ 1.,  2., 99.],
       [ 2.,  4., 99.]])
    Coordinates:
      * lat      (lat) int64 0 1
      * lon      (lon) int64 0 50 360

    See Also
    --------
    Related External Functions:
    `xarray.DataArray.interp <https://docs.xarray.dev/en/stable/generated/xarray.DataArray.interp.html>`__,
    `cartopy.util.add_cyclic_point <https://scitools.org.uk/cartopy/docs/latest/reference/generated/cartopy.util.add_cyclic_point.html>`__

    Related NCL Function:
    `NCL linint2 <https://www.ncl.ucar.edu/Document/Functions/Built-in/linint2.shtml>`__
    """
    # check for xarray/numpy
    if not isinstance(data_in, xr.DataArray):
        if lat_in is None or lon_in is None:
            raise ValueError(
                "Argument lat_in and lon_in must be provided if data_in is not an xarray"
            )
        data_in = xr.DataArray(
            data_in, dims=["lat", "lon"], coords={"lat": lat_in, "lon": lon_in}
        )

    output_coords = {
        data_in.dims[-1]: lon_out,
        data_in.dims[-2]: lat_out,
    }

    data_in_modified = _pre_interp_multidim(data_in, cyclic, missing_val)
    data_out = data_in_modified.interp(
        output_coords, method=method, kwargs={"fill_value": fill_value}
    )
    data_out_modified = _post_interp_multidim(data_out, missing_val=missing_val)

    return data_out_modified
