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

import metpy.interpolate
import numpy as np
import xarray as xr

try:
    from .fortran.vinth2p_backend import dsigma2hybrid_nodes as _dsigma2hybrid_nodes
    from .fortran.vinth2p_backend import (
        dsigma2hybrid_nodes_into as _dsigma2hybrid_nodes_into,
    )
    from .fortran.vinth2p_backend import (
        dvinth2p_ecmwf_nodes_corder_pa_into as _dvinth2p_ecmwf_nodes_corder_pa_into,
    )
    from .fortran.vinth2p_backend import (
        dvinth2p_ecmwf_nodes_pa as _dvinth2p_ecmwf_nodes_pa,
    )
    from .fortran.vinth2p_backend import (
        dvinth2p_ecmwf_nodes_pa_into as _dvinth2p_ecmwf_nodes_pa_into,
    )
    from .fortran.vinth2p_backend import (
        dvinth2p_nodes_corder_pa_into as _dvinth2p_nodes_corder_pa_into,
    )
    from .fortran.vinth2p_backend import dvinth2p_nodes_pa as _dvinth2p_nodes_pa
    from .fortran.vinth2p_backend import (
        dvinth2p_nodes_pa_into as _dvinth2p_nodes_pa_into,
    )
except Exception:
    _dsigma2hybrid_nodes = None
    _dsigma2hybrid_nodes_into = None
    _dvinth2p_nodes_pa = None
    _dvinth2p_ecmwf_nodes_pa = None
    _dvinth2p_nodes_pa_into = None
    _dvinth2p_ecmwf_nodes_pa_into = None
    _dvinth2p_nodes_corder_pa_into = None
    _dvinth2p_ecmwf_nodes_corder_pa_into = None

__all__ = [
    "pressure_at_hybrid_levels",
    "delta_pressure_hybrid",
    "interp_hybrid_to_pressure",
    "interp_sigma_to_hybrid",
    "interp_multidim",
]

supported_types = typing.Union[xr.DataArray, np.ndarray]
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


def _func_interpolate(method="linear"):
    """Define interpolation function."""

    if method == "linear":
        func_interpolate = metpy.interpolate.interpolate_1d
    elif method == "log":
        func_interpolate = metpy.interpolate.log_interpolate_1d
    else:
        raise ValueError(
            f"Unknown interpolation method: {method}. "
            f'Supported methods are: "log" and "linear".'
        )

    return func_interpolate


def _interpolate_mb(data, curr_levels, new_levels, axis, method="linear"):
    """Wrapper used by ``xarray.map_blocks`` for vertical interpolation."""

    func_interpolate = _func_interpolate(method)
    return func_interpolate(new_levels, curr_levels, data, axis=axis)


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


def pressure_at_hybrid_levels(psfc, hya, hyb, p0=100000.0):
    """Calculate pressure at hybrid levels.

    Parameters
    ----------
    psfc : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        Surface pressure in Pascals.

    hya, hyb : :class:`xarray.DataArray`, :class:`numpy.ndarray`
        One-dimensional hybrid coefficients.

    p0 : float, optional
        Reference pressure in Pascals.

    Returns
    -------
    :class:`xarray.DataArray`, :class:`numpy.ndarray`
        Pressure at hybrid levels in Pascals.
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

    if isinstance(hya, np.ndarray):
        if hya.ndim != 1:
            raise ValueError(
                f"hya and hyb must be 1-dimensional if numpy inputs: {hya.shape}"
            )
        reshape = (hya.shape[0],) + (1,) * np.ndim(psfc)
        hya = hya.reshape(reshape)
        hyb = hyb.reshape(reshape)
        return hya * p0 + hyb * psfc

    if hya.dims != hyb.dims:
        warnings.warn(
            "hya and hyb have different dimension names, attempting rename",
            stacklevel=2,
        )
        hyb = hyb.rename({b: a for a, b in zip(hya.dims, hyb.dims)})

    hya, hyb = _rename_colliding_coeff_dim(psfc, hya, hyb)

    # p(k) = hya(k) * p0 + hyb(k) * psfc

    # This will be in Pa
    return hya * p0 + hyb * psfc


def delta_pressure_hybrid(ps, hya, hyb, p0=100000.0):
    """Calculate pressure layer thickness for hybrid coordinates."""

    if not all(isinstance(x, (xr.DataArray, np.ndarray)) for x in (ps, hya, hyb)):
        raise TypeError("Inputs must be xarray DataArrays or numpy arrays")

    if not isinstance(p0, (float, int, np.floating, np.integer)):
        raise TypeError(f"p0 must be a scalar numeric value, received {type(p0)}")

    if hya.shape != hyb.shape:
        raise ValueError(f"dimension mismatch: hya: {hya.shape} hyb: {hyb.shape}")

    if np.ndim(hya) != 1:
        raise ValueError(f"hya and hyb must be 1-dimensional: {hya.shape}")

    if isinstance(ps, np.ndarray):
        hya_values = np.asarray(hya.data if isinstance(hya, xr.DataArray) else hya)
        hyb_values = np.asarray(hyb.data if isinstance(hyb, xr.DataArray) else hyb)
        reshape = (hya_values.shape[0] - 1,) + (1,) * ps.ndim
        pa = (
            p0 * hya_values[:-1].reshape(reshape)
            + hyb_values[:-1].reshape(reshape) * ps
        )
        pb = p0 * hya_values[1:].reshape(reshape) + hyb_values[1:].reshape(reshape) * ps
        return np.abs(pa - pb)

    if isinstance(hya, np.ndarray):
        hya = xr.DataArray(hya, dims=("lev",))
        hyb = xr.DataArray(hyb, dims=("lev",))
    else:
        hya = xr.DataArray(np.asarray(hya.data), dims=hya.dims)
        hyb = xr.DataArray(np.asarray(hyb.data), dims=hyb.dims)
        if hya.dims != hyb.dims:
            warnings.warn(
                "hya and hyb have different dimension names, attempting rename",
                stacklevel=2,
            )
            hyb = hyb.rename({b: a for a, b in zip(hya.dims, hyb.dims)})

    hya, hyb = _rename_colliding_coeff_dim(ps, hya, hyb)

    lev_name = hya.dims[0]
    pa = (
        p0 * hya.isel({lev_name: slice(None, -1)})
        + hyb.isel({lev_name: slice(None, -1)}) * ps
    )
    pb = (
        p0 * hya.isel({lev_name: slice(1, None)})
        + hyb.isel({lev_name: slice(1, None)}) * ps
    )

    dph = abs(pa - pb)
    dph.name = "dph"
    dph.attrs = {
        "long_name": "pressure layer thickness",
        "units": "Pa",
    }
    return dph


def _pressure_from_hybrid(psfc, hya, hyb, p0=100000.0):
    """Backward-compatible wrapper for :func:`pressure_at_hybrid_levels`."""

    return pressure_at_hybrid_levels(psfc, hya, hyb, p0)


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
        data_in = xr.DataArray(
            np.where(data_in.values == missing_val, np.nan, data_in.values),
            dims=data_in.dims,
            coords=data_in.coords,
        )

    # add cyclic points and create new data array
    if cyclic:
        padded_data = np.pad(data_in.values, ((0, 0), (1, 1)), mode="wrap")
        padded_longitudes = np.pad(
            data_in.coords[data_in.dims[-1]], (1, 1), mode="wrap"
        )
        padded_longitudes[0] -= 360
        padded_longitudes[-1] += 360

        data_in = xr.DataArray(
            padded_data,
            coords={
                data_in.dims[-2]: data_in.coords[data_in.dims[-2]].values,
                data_in.dims[-1]: padded_longitudes,
            },
            dims=data_in.dims,
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
        data_in = xr.DataArray(
            np.where(np.isnan(data_in.values), missing_val, data_in.values),
            dims=data_in.dims,
            coords=data_in.coords,
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


def _vertical_remap(func_interpolate, new_levels, xcoords, data, interp_axis=0):
    """Execute the defined interpolation function on data."""

    if (
        not isinstance(xcoords, xr.DataArray)
        and np.ndim(xcoords) == 1
        and np.ndim(data) > 1
    ):
        reshape = [1] * np.ndim(data)
        reshape[interp_axis] = np.shape(xcoords)[0]
        xcoords = np.reshape(xcoords, tuple(reshape))

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", r"Interpolation point out of data bounds encountered"
        )
        return func_interpolate(new_levels, xcoords, data, axis=interp_axis)


def _temp_extrapolate(data_or_t_bot, *args):
    r"""Extrapolate temperature below ground.

    This helper accepts both the upstream-style call
    ``_temp_extrapolate(t_bot, lev, p_sfc, ps, phi_sfc)`` and the legacy
    Skyborn call ``_temp_extrapolate(data, lev_dim, lev, p_sfc, ps, phi_sfc)``.

    .. math::
        T = T_* \left( 1 + \alpha ln \frac{p}{p_s} + \frac{1}{2}\left( \alpha ln \frac{p}{p_s} \right)^2 + \frac{1}{6} \left( \alpha ln \frac{p}{p_s} \right)^3 \right)
    """
    if len(args) == 4:
        lev, p_sfc, ps, phi_sfc = args
        t_bot = data_or_t_bot
    elif len(args) == 5 and isinstance(args[0], str):
        lev_dim, lev, p_sfc, ps, phi_sfc = args
        t_bot = data_or_t_bot.isel({lev_dim: -1}, drop=True)
    else:
        raise TypeError(
            "_temp_extrapolate accepts either "
            "(t_bot, lev, p_sfc, ps, phi_sfc) or "
            "(data, lev_dim, lev, p_sfc, ps, phi_sfc)"
        )

    R_d = 287.04  # dry air gas constant
    g_inv = 1 / 9.80616  # inverse of gravity
    alpha = 0.0065 * R_d * g_inv

    tstar = t_bot * (1 + alpha * (ps / p_sfc - 1))
    hgt = phi_sfc * g_inv
    t0 = tstar + 0.0065 * hgt
    tplat = xr.apply_ufunc(np.minimum, 298, t0, dask="parallelized")

    tprime0 = xr.where(
        (2000 <= hgt) & (hgt <= 2500),
        0.002 * ((2500 - hgt) * t0 + ((hgt - 2000) * tplat)),
        np.nan,
    )
    tprime0 = xr.where(2500 < hgt, tplat, tprime0)

    alnp = xr.where(
        hgt < 2000,
        alpha * np.log(lev / ps),
        R_d * (tprime0 - tstar) / phi_sfc * np.log(lev / ps),
    )
    alnp = xr.where(tprime0 < tstar, 0, alnp)

    return tstar * (1 + alnp + (0.5 * (alnp**2)) + (1 / 6 * (alnp**3)))


def _geo_height_extrapolate(t_bot, lev, p_sfc, ps, phi_sfc):
    r"""This helper function extrapolates geopotential height below ground using
    the ECMWF formulation described in `Vertical Interpolation and Truncation
    of Model-Coordinate Data <https://dx.doi.org/10.5065/D6HX19NH>`__ by
    Trenberth, Berry, & Buja [NCAR/TN-396, 1993]. Specifically equation 15 is
    used:

    .. math::
        \Phi = \Phi_s - R_d T_* ln \frac{p}{p_s} \left[ 1 + \frac{1}{2}\alpha ln\frac{p}{p_s} + \frac{1}{6} \left( \alpha ln \frac{p}{p_s} \right)^2 \right]

    Parameters
    ----------
    t_bot: :class:`xarray.DataArray`
        Temperature at the lowest (bottom) level of the model.

    lev: int
        The pressure level of interest. Must be in the same units as ``ps`` and ``p_sfc``

    p_sfc: :class:`xarray.DataArray`
        The pressure at the lowest level of the model. Must be in the same units as ``lev`` and ``ps``

    ps : :class:`xarray.DataArray`
        An array of surface pressures. Must be in the same units as ``lev`` and ``p_sfc``

    phi_sfc:
        The geopotential at the lowest level of the model.

    Returns
    -------
    result: :class:`xarray.DataArray`
        The extrapolated geopotential height in geopotential meters at the provided pressure levels.
    """
    R_d = 287.04  # dry air gas constant
    g_inv = 1 / 9.80616  # inverse of gravity
    alpha = 0.0065 * R_d * g_inv

    tstar = t_bot * (1 + alpha * (ps / p_sfc - 1))
    hgt = phi_sfc * g_inv
    t0 = tstar + 0.0065 * hgt

    alph = xr.where(
        (tstar <= 290.5) & (t0 > 290.5), R_d / phi_sfc * (290.5 - tstar), alpha
    )

    alph = xr.where((tstar > 290.5) & (t0 > 290.5), 0, alph)
    tstar = xr.where((tstar > 290.5) & (t0 > 290.5), 0.5 * (290.5 + tstar), tstar)

    tstar = xr.where((tstar < 255), 0.5 * (tstar + 255), tstar)

    alnp = alph * np.log(lev / ps)
    return hgt - R_d * tstar * g_inv * np.log(lev / ps) * (
        1 + 0.5 * alnp + 1 / 6 * alnp**2
    )


def _vertical_remap_extrap(
    new_levels, lev_dim, data, output, pressure, ps, variable, t_bot, phi_sfc
):
    """A helper function to call the appropriate extrapolation function based
    on the user's inputs.

    Parameters
    ----------
    new_levels: array-like
        The desired pressure levels for extrapolation in Pascals.

    lev_dim: str
        The name of the vertical dimension.

    data: :class:`xarray.DataArray`
        The data to extrapolate

    output: :class:`xarray.DataArray`
        An array to hold the output data

    pressure: :class:`xarray.DataArray`
        The pressure at the lowest level of the model. Must be in the same units as ``lev`` and ``ps``

    ps : :class:`xarray.DataArray`
        An array of surface pressures. Must be in the same units as ``lev`` and ``p_sfc``

    variable : str, optional
        String representing what variable is extrapolated below surface level.
        Temperature extrapolation = "temperature". Geopotential height
        extrapolation = "geopotential". All other variables = "other". If
        "other", the value of ``data`` at the lowest model level will be used
        as the below ground fill value. Required if extrapolate is True.

    t_bot: :class:`xarray.DataArray`
        Temperature at the lowest (bottom) level of the model.

    phi_sfc:
        The geopotential at the lowest level of the model.

    Returns
    -------
    output: :class:`xarray.DataArray`
        A DataArray containing the data after extrapolation.
    """

    sfc_index = pressure[lev_dim].argmax(dim=lev_dim)  # index of the model surface
    p_sfc = pressure.isel(
        {lev_dim: sfc_index}, drop=True
    )  # extract pressure at lowest level

    if variable == "temperature":
        output = output.where(
            output.plev <= p_sfc,
            _temp_extrapolate(t_bot, output.plev, p_sfc, ps, phi_sfc),
        )
    elif variable == "geopotential":
        output = output.where(
            output.plev <= p_sfc,
            _geo_height_extrapolate(t_bot, output.plev, p_sfc, ps, phi_sfc),
        )
    else:
        output = output.where(
            output.plev <= p_sfc, data.isel({lev_dim: sfc_index}, drop=True)
        )

    return output


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


def _strip_unexpected_pint_units(output, in_pint):
    """Remove unexpected pint wrapping introduced by MetPy."""

    if in_pint:
        return output

    if hasattr(output.data, "magnitude"):
        output.data = output.data.magnitude

    return output


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

    if method == "linear":
        return 1
    if method == "log":
        return 2

    raise ValueError(
        f"Unknown interpolation method: {method}. "
        f'Supported methods are: "log" and "linear".'
    )


def _vinth2p_varflg(variable: str) -> int:
    """Translate public extrapolation variable labels to ECMWF vinth flags."""

    if variable == "temperature":
        return 1
    if variable == "geopotential":
        return -1
    return 0


def _as_c_contiguous_float64_view(array: xr.DataArray, dims, shape):
    """Return a zero-copy 1D float64 view when the DataArray matches the fast path."""

    if array is None or not isinstance(array, xr.DataArray):
        return None

    if tuple(array.dims) != tuple(dims) or tuple(array.shape) != tuple(shape):
        return None

    values = array.data
    if not isinstance(values, np.ndarray):
        return None
    if values.dtype != np.float64 or not values.flags.c_contiguous:
        return None

    return values.reshape(-1)


def _as_broadcast_float64_flat(array: xr.DataArray, template: xr.DataArray):
    """Return a flattened float64 array after lightweight xarray broadcasting."""

    if array is None or not isinstance(array, xr.DataArray):
        return None

    try:
        broadcast = array.broadcast_like(template)
    except Exception:
        return None

    values = np.asarray(broadcast.data, dtype=np.float64)
    return np.ascontiguousarray(values).reshape(-1)


def _build_vinth2p_output(data, interp_axis, new_levels, output_values, base_template):
    """Wrap a NumPy output array in the public xarray result."""

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


def _interp_hybrid_to_pressure_fortran_corder(
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
    if raw_data.dtype != np.float64 or not raw_data.flags.c_contiguous:
        return None

    shape_before = data.shape[:interp_axis]
    shape_after = data.shape[interp_axis + 1 :]
    lead_dims = data.dims[:interp_axis] + data.dims[interp_axis + 1 :]
    lead_shape = shape_before + shape_after
    nouter = int(np.prod(shape_before, dtype=np.int64)) if shape_before else 1
    ninner = int(np.prod(shape_after, dtype=np.int64)) if shape_after else 1
    nlevi = data.shape[interp_axis]
    nlevo = new_levels.size

    ps_flat = _as_c_contiguous_float64_view(ps, lead_dims, lead_shape)
    if ps_flat is None or not np.isfinite(ps_flat).all():
        return None

    hyam_values = np.asarray(hyam.data, dtype=np.float64).reshape(nlevi)
    hybm_values = np.asarray(hybm.data, dtype=np.float64).reshape(nlevi)
    new_level_values = np.asarray(new_levels, dtype=np.float64).reshape(nlevo)
    intyp = _vinth2p_intyp(method)
    base_template = data.isel({lev_dim: 0}, drop=True)

    output_values = np.empty(
        (*shape_before, nlevo, *shape_after), dtype=np.float64, order="C"
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
            tbot_flat = _as_c_contiguous_float64_view(t_bot, lead_dims, lead_shape)
            if tbot_flat is None:
                tbot_flat = _as_broadcast_float64_flat(t_bot, base_template)
            phi_flat = _as_c_contiguous_float64_view(phi_sfc, lead_dims, lead_shape)
            if phi_flat is None:
                phi_flat = _as_broadcast_float64_flat(phi_sfc, base_template)
            if tbot_flat is None or phi_flat is None:
                return None
        _dvinth2p_ecmwf_nodes_corder_pa_into(
            raw_data.reshape(-1),
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
            raw_data.reshape(-1),
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


def _interp_hybrid_to_pressure_fortran(
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

    corder_output = _interp_hybrid_to_pressure_fortran_corder(
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

    data_columns = np.asfortranarray(
        np.asarray(data_view.data, dtype=np.float64).reshape(ncol, nlevi).T
    )
    ps_columns = np.asarray(
        ps.broadcast_like(base_template).transpose(*lead_dims).data,
        dtype=np.float64,
    ).reshape(ncol)
    ps_columns = np.where(np.isfinite(ps_columns), ps_columns, _VINTH2P_SPVL)

    hyam_values = np.asarray(hyam.data, dtype=np.float64).reshape(nlevi)
    hybm_values = np.asarray(hybm.data, dtype=np.float64).reshape(nlevi)
    new_level_values = np.asarray(new_levels, dtype=np.float64).reshape(nlevo)
    intyp = _vinth2p_intyp(method)
    output_columns = np.empty((nlevo, ncol), dtype=np.float64, order="F")

    if extrapolate:
        varflg = _vinth2p_varflg(variable)
        if varflg == 0:
            t_bot_columns = np.zeros(ncol, dtype=np.float64)
            phi_columns = np.zeros(ncol, dtype=np.float64)
        else:
            t_bot_columns = np.asarray(
                t_bot.broadcast_like(base_template).transpose(*lead_dims).data,
                dtype=np.float64,
            ).reshape(ncol)
            phi_columns = np.asarray(
                phi_sfc.broadcast_like(base_template).transpose(*lead_dims).data,
                dtype=np.float64,
            ).reshape(ncol)
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

    output_values = np.asarray(output_columns, dtype=np.float64).T.reshape(
        (*lead_shape, nlevo)
    )
    output_values[output_values == _VINTH2P_SPVL] = np.nan

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


def _interp_sigma_to_hybrid_fortran(
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

    sigma_source_values = np.asarray(
        sig_coords.data if isinstance(sig_coords, xr.DataArray) else sig_coords,
        dtype=np.float64,
    ).reshape(-1)
    sigma_diffs = np.diff(sigma_source_values)
    if not (
        np.all(sigma_diffs >= 0.0)
        or np.all(sigma_diffs <= 0.0)
        or sigma_source_values.size <= 1
    ):
        raise ValueError("sigma2hybrid backend requires monotonic sigma coordinates")

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

    data_columns = np.asfortranarray(
        np.asarray(data_view.data, dtype=np.float64).reshape(ncol, nlevi).T
    )
    sigma_target_columns = np.asfortranarray(
        np.asarray(sigma_target_view.data, dtype=np.float64).reshape(nlevo, ncol)
    )

    intyp = _vinth2p_intyp(method)
    output_columns = np.empty((nlevo, ncol), dtype=np.float64, order="F")

    if _dsigma2hybrid_nodes_into is not None:
        _dsigma2hybrid_nodes_into(
            data_columns,
            output_columns,
            sigma_source_values,
            sigma_target_columns,
            intyp,
            _VINTH2P_SPVL,
        )
    else:
        output_columns = _dsigma2hybrid_nodes(
            data_columns,
            sigma_source_values,
            sigma_target_columns,
            intyp,
            _VINTH2P_SPVL,
        )

    output_values = np.asarray(output_columns, dtype=np.float64).T.reshape(
        (*lead_shape, nlevo)
    )
    output_values[output_values == _VINTH2P_SPVL] = np.nan

    h_coords = (
        sigma_target_columns[:, 0].copy() if ncol else np.asarray([], dtype=np.float64)
    )
    coords = {k: v for k, v in base_template.coords.items()}
    coords["hlev"] = h_coords
    output = xr.DataArray(
        output_values,
        dims=(*lead_dims, "hlev"),
        coords=coords,
        name=data.name,
        attrs=data.attrs,
    )

    dims = [data.dims[i] if i != interp_axis else "hlev" for i in range(data.ndim)]
    return output.transpose(*dims)


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
        String that is the interpolation method; can be either "linear" or "log".
        Defaults to "linear".

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

    new_levels = np.asarray(new_levels)
    in_pint = any(_is_pint_backed(arr) for arr in (data, ps, hyam, hybm))
    in_dask = any(_is_dask_backed(arr) for arr in (data, ps, hyam, hybm))

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

    try:
        func_interpolate = _func_interpolate(method)
    except ValueError as vexc:
        raise ValueError(vexc.args[0])

    if (
        not in_dask
        and not in_pint
        and _dvinth2p_nodes_pa is not None
        and method in {"linear", "log"}
    ):
        return _interp_hybrid_to_pressure_fortran(
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

    interp_axis = data.dims.index(lev_dim)

    # Calculate pressure levels at the hybrid levels
    pressure = pressure_at_hybrid_levels(ps, hyam, hybm, p0)  # Pa

    # Make pressure shape same as data shape
    pressure = pressure.transpose(*data.dims)
    output = None

    chunk_sizes = getattr(data, "chunksizes", None) or {}
    if _is_dask_backed(data) and lev_dim in chunk_sizes:
        if len(chunk_sizes[lev_dim]) == 1:
            try:
                output = xr.map_blocks(
                    _interpolate_mb,
                    data,
                    args=(pressure, new_levels, interp_axis, method),
                )
            except Exception:
                output = None
        else:
            warnings.warn(
                f"Chunking along {lev_dim} is not recommended for performance reasons.",
                stacklevel=2,
            )

    if in_dask and output is None:
        from dask.array.core import map_blocks

        if not _is_dask_backed(data):
            data = data.chunk({dim: size for dim, size in data.sizes.items()})

        pressure = pressure.chunk(data.chunksizes)
        out_chunks = list(data.chunks)
        out_chunks[interp_axis] = (new_levels.size,)
        out_chunks = tuple(out_chunks)

        output = map_blocks(
            _vertical_remap,
            func_interpolate,
            new_levels,
            pressure.data,
            data.data,
            interp_axis,
            chunks=out_chunks,
            dtype=data.dtype,
            drop_axis=[interp_axis],
            new_axis=[interp_axis],
        )
    elif output is None:
        output = _vertical_remap(
            func_interpolate, new_levels, pressure.data, data.data, interp_axis
        )

    if isinstance(output, xr.DataArray):
        output = output.copy(deep=False)
        output.name = data.name
        output.attrs = data.attrs
    else:
        output = xr.DataArray(output, name=data.name, attrs=data.attrs)

    output = _strip_unexpected_pint_units(output, in_pint)

    # Set output dims and coords
    dims = [data.dims[i] if i != interp_axis else "plev" for i in range(data.ndim)]

    # Rename output dims. This is only needed with above workaround block
    dims_dict = {output.dims[i]: dims[i] for i in range(len(output.dims))}
    output = output.rename(dims_dict)

    coords = {}
    for k, v in data.coords.items():
        if k != lev_dim:
            coords.update({k: v})
        else:
            coords.update({"plev": new_levels})

    output = output.transpose(*dims).assign_coords(coords)

    if extrapolate:
        output = _vertical_remap_extrap(
            new_levels, lev_dim, data, output, pressure, ps, variable, t_bot, phi_sfc
        )
        output = _strip_unexpected_pint_units(output, in_pint)

    return output


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
        String that is the interpolation method; can be either "linear" or "log".
        Defaults to "linear".

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

    # Determine the level dimension and then the interpolation axis
    if lev_dim is None:
        try:
            lev_dim = data.cf["vertical"].name
        except Exception:
            raise ValueError(
                "Unable to determine vertical dimension name. Please specify the name via `lev_dim` argument.'"
            )

    try:
        func_interpolate = _func_interpolate(method)
    except ValueError as vexc:
        raise ValueError(vexc.args[0])

    in_dask = _is_dask_backed(data) or _is_dask_backed(ps)
    in_pint = _is_pint_backed(data) or _is_pint_backed(ps)

    if (
        not in_dask
        and not in_pint
        and _dsigma2hybrid_nodes is not None
        and method in {"linear", "log"}
    ):
        sigma_source_values = np.asarray(
            sig_coords.data if isinstance(sig_coords, xr.DataArray) else sig_coords,
            dtype=np.float64,
        ).reshape(-1)
        sigma_diffs = np.diff(sigma_source_values)
        if (
            np.all(sigma_diffs >= 0.0)
            or np.all(sigma_diffs <= 0.0)
            or sigma_diffs.size == 0
        ):
            return _interp_sigma_to_hybrid_fortran(
                data=data,
                sig_coords=sig_coords,
                ps=ps,
                hyam=hyam,
                hybm=hybm,
                p0=p0,
                lev_dim=lev_dim,
                method=method,
            )

    # Calculate sigma levels at the hybrid levels
    sigma = _sigma_from_hybrid(ps, hyam, hybm, p0)  # Pa
    sig_coord_values = (
        sig_coords.data if isinstance(sig_coords, xr.DataArray) else sig_coords
    )

    non_lev_dims = list(data.dims)
    if data.ndim > 1:
        non_lev_dims.remove(lev_dim)
        data_stacked = data.stack(combined=non_lev_dims).transpose()
        sigma_stacked = sigma.stack(combined=non_lev_dims).transpose()

        h_coords = sigma_stacked[0, :].copy()

        output = data_stacked[:, : len(hyam)].copy()

        for idx, (d, s) in enumerate(zip(data_stacked, sigma_stacked)):
            output[idx, :] = xr.DataArray(
                _vertical_remap(func_interpolate, s.data, sig_coord_values, d.data),
                dims=[lev_dim],
            )

        # Make output shape same as data shape
        output = output.unstack().transpose(*data.dims)
    else:
        h_coords = sigma

        output = data[: len(hyam)].copy()
        output[: len(hyam)] = xr.DataArray(
            _vertical_remap(func_interpolate, sigma.data, sig_coord_values, data.data),
            dims=[lev_dim],
        )

    # Set output dims and coords
    output = output.rename({lev_dim: "hlev"})
    output = output.assign_coords({"hlev": h_coords.data})

    return output


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
