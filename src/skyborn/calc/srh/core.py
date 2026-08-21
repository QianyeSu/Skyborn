"""Storm-relative helicity calculations backed by a compiled Fortran/C extension.

Inputs (plain floats/arrays, all MKS or angle-free):

* ``height`` in meters (converted to AGL internally)
* ``u``, ``v`` in m/s
* ``depth`` in meters (layer depth, required)
* ``bottom`` in meters AGL (default 0)
* ``storm_u``, ``storm_v`` in m/s (default 0)

Outputs (m^2/s^2): positive SRH, negative SRH, and total SRH. The algorithm
follows metpy's ``storm_relative_helicity``: heights are converted to AGL,
the layer [bottom, bottom+depth] is subsampled with the boundary heights
linearly interpolated (NaN outside the data range, matching metpy), and the
hodograph cross terms are summed separately for positive/negative/zero.

1D input returns three floats; 3D input (level, lat, lon) returns three 2D
fields, wrapped as xarray DataArrays when given xarray input.
"""

from __future__ import annotations

from importlib import import_module
from typing import Union

import numpy as np
import xarray as xr

__all__ = [
    "calculate_storm_relative_helicity",
    "srh_profile",
    "srh_grid",
]

_BACKEND = import_module(f"{__package__}.srh_core")


def _as_numpy(values):
    return getattr(values, "values", values)


def _wrap_xarray_result(height: xr.DataArray, values: np.ndarray) -> xr.DataArray:
    dims = list(height.dims[-2:])
    coords = {}
    for dim in dims:
        if dim in height.coords:
            coords[dim] = height.coords[dim]
    return xr.DataArray(
        values,
        dims=dims or ["y", "x"],
        coords=coords,
        attrs={"units": "m2 s-2", "long_name": "storm-relative helicity"},
    )


def srh_profile(
    height: Union[np.ndarray, xr.DataArray],
    u: Union[np.ndarray, xr.DataArray],
    v: Union[np.ndarray, xr.DataArray],
    depth: float,
    bottom: float = 0.0,
    storm_u: float = 0.0,
    storm_v: float = 0.0,
):
    """Compute (positive, negative, total) SRH for a single profile."""
    h = np.asarray(_as_numpy(height), dtype=np.float64)
    ua = np.asarray(_as_numpy(u), dtype=np.float64)
    va = np.asarray(_as_numpy(v), dtype=np.float64)
    if h.ndim != 1:
        raise ValueError(f"height must be 1D, got shape {h.shape}")
    if ua.shape != h.shape or va.shape != h.shape:
        raise ValueError("height, u, and v must share a shape")
    return tuple(
        float(x)
        for x in _BACKEND.srh_profile(
            np.ascontiguousarray(h),
            np.ascontiguousarray(ua),
            np.ascontiguousarray(va),
            float(depth),
            float(bottom),
            float(storm_u),
            float(storm_v),
        )
    )


def srh_grid(
    height: Union[np.ndarray, xr.DataArray],
    u: Union[np.ndarray, xr.DataArray],
    v: Union[np.ndarray, xr.DataArray],
    depth: float,
    bottom: float = 0.0,
    storm_u: float = 0.0,
    storm_v: float = 0.0,
):
    """Compute (positive, negative, total) SRH for a 3D grid of profiles."""
    h = np.asarray(_as_numpy(height), dtype=np.float64)
    ua = np.asarray(_as_numpy(u), dtype=np.float64)
    va = np.asarray(_as_numpy(v), dtype=np.float64)
    if h.ndim != 3:
        raise ValueError(f"height must be 3D, got shape {h.shape}")
    if ua.shape != h.shape or va.shape != h.shape:
        raise ValueError("height, u, and v must share a shape")
    return tuple(
        np.asarray(arr, dtype=np.float64)
        for arr in _BACKEND.srh_grid(
            np.asfortranarray(h),
            np.asfortranarray(ua),
            np.asfortranarray(va),
            float(depth),
            float(bottom),
            float(storm_u),
            float(storm_v),
        )
    )


def calculate_storm_relative_helicity(
    height: Union[np.ndarray, xr.DataArray],
    u: Union[np.ndarray, xr.DataArray],
    v: Union[np.ndarray, xr.DataArray],
    depth: float,
    *,
    bottom: float = 0.0,
    storm_u: float = 0.0,
    storm_v: float = 0.0,
):
    """Calculate storm-relative helicity.

    Parameters
    ----------
    height : np.ndarray or xr.DataArray
        Atmospheric height in meters. Either 1D (level,) or 3D (level, lat, lon).
    u : np.ndarray or xr.DataArray
        U-component wind in m/s, same shape as ``height``.
    v : np.ndarray or xr.DataArray
        V-component wind in m/s, same shape as ``height``.
    depth : float
        Depth of the layer in meters.
    bottom : float, optional
        Height of the layer bottom AGL in meters (default 0, i.e. surface).
    storm_u, storm_v : float, optional
        U/V components of storm motion in m/s (default 0).

    Returns
    -------
    tuple of (positive_srh, negative_srh, total_srh)
        For 1D input: three floats in m^2/s^2. For 3D input: three 2D arrays
        (lat, lon), wrapped as ``xr.DataArray`` for xarray input.
    """
    h = _as_numpy(height)
    ua = _as_numpy(u)
    va = _as_numpy(v)
    is_xarray = hasattr(height, "attrs")

    if np.ndim(h) == 1:
        return srh_profile(h, ua, va, depth, bottom, storm_u, storm_v)

    if np.ndim(h) == 3:
        out = srh_grid(h, ua, va, depth, bottom, storm_u, storm_v)
        if is_xarray:
            out = tuple(_wrap_xarray_result(height, arr) for arr in out)
        return out

    raise ValueError(f"height must be 1D or 3D, got shape {np.shape(h)}")
