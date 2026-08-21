"""CAPE/CIN calculations backed by the compiled Fortran/C extension.

This module mirrors the ``skyborn.calc.dcape`` layout: a compiled backend
handles the per-profile numerics, while this Python layer handles shape
dispatch (1D profile -> scalars, 3D grid -> 2D fields) and optional
xarray metadata preservation.

Input convention (same as the DCAPE backend):

* ``pressure``    in hPa
* ``temperature`` in Celsius
* ``dewpoint``    in Celsius

Outputs:

* ``cape`` / ``cin``  in J kg^-1 (CIN is reported <= 0, matching metpy)
* ``lfc_p`` / ``el_p`` in hPa
* parcel profile temperatures in Celsius

The parcel is a surface-based parcel: dry adiabatic ascent to the
Bolton-approximation LCL, then a fixed-step RK4 moist pseudo-adiabat.
CAPE/CIN follow metpy's convention (``which_lfc='bottom'``,
``which_el='top'``), integrated against ln(p) with zero crossings
inserted by log-pressure linear interpolation.
"""

from __future__ import annotations

from importlib import import_module
from typing import Union

import numpy as np
import xarray as xr

__all__ = [
    "calculate_cape_cin",
    "cape_profile",
    "cape_grid",
    "calculate_parcel_profile",
    "parcel_profile",
    "parcel_profile_grid",
    "calculate_most_unstable_parcel",
    "most_unstable_parcel",
    "most_unstable_parcel_grid",
    "calculate_most_unstable_cape_cin",
    "most_unstable_cape_cin",
    "most_unstable_cape_cin_grid",
]

_BACKEND = import_module(f"{__package__}.cape_core")


def _as_numpy(values):
    """Extract a NumPy array from xarray-like or array-like inputs."""
    return getattr(values, "values", values)


def _wrap_xarray_result(pressure: xr.DataArray, values: np.ndarray) -> xr.DataArray:
    """Preserve xarray metadata for a 2D CAPE result field."""
    dims = list(pressure.dims[-2:])
    coords = {}
    for dim in dims:
        if dim in pressure.coords:
            coords[dim] = pressure.coords[dim]
    return xr.DataArray(
        values,
        dims=dims or ["y", "x"],
        coords=coords,
        attrs={"units": "J kg-1", "long_name": "CAPE/CIN"},
    )


def _wrap_xarray_3d(pressure: xr.DataArray, values: np.ndarray) -> xr.DataArray:
    """Preserve xarray metadata for a 3D parcel-profile field."""
    dims = list(pressure.dims)
    coords = {}
    for dim in dims:
        if dim in pressure.coords:
            coords[dim] = pressure.coords[dim]
    return xr.DataArray(
        values,
        dims=dims,
        coords=coords,
        attrs={"units": "degree_Celsius", "long_name": "parcel temperature"},
    )


def cape_profile(
    pressure: Union[np.ndarray, xr.DataArray],
    temperature: Union[np.ndarray, xr.DataArray],
    dewpoint: Union[np.ndarray, xr.DataArray],
):
    """Compute (cape, cin, lfc_p, el_p) for a single vertical profile.

    Returns four Python floats: CAPE (J/kg), CIN (J/kg, <= 0), LFC
    pressure (hPa), and EL pressure (hPa). Missing levels return NaN.
    """
    p = np.asarray(_as_numpy(pressure), dtype=np.float64)
    t = np.asarray(_as_numpy(temperature), dtype=np.float64)
    td = np.asarray(_as_numpy(dewpoint), dtype=np.float64)
    if p.ndim != 1:
        raise ValueError(f"pressure must be 1D, got shape {p.shape}")
    if t.shape != p.shape or td.shape != p.shape:
        raise ValueError("pressure, temperature, and dewpoint must share a shape")
    return tuple(
        float(v)
        for v in _BACKEND.cape_profile(
            np.ascontiguousarray(p),
            np.ascontiguousarray(t),
            np.ascontiguousarray(td),
        )
    )


def cape_grid(
    pressure: Union[np.ndarray, xr.DataArray],
    temperature: Union[np.ndarray, xr.DataArray],
    dewpoint: Union[np.ndarray, xr.DataArray],
):
    """Compute (cape, cin, lfc_p, el_p) for a 3D grid of vertical profiles.

    Expects (level, lat, lon) layout; returns four 2D (lat, lon) arrays.
    """
    p = np.asarray(_as_numpy(pressure), dtype=np.float64)
    t = np.asarray(_as_numpy(temperature), dtype=np.float64)
    td = np.asarray(_as_numpy(dewpoint), dtype=np.float64)
    if p.ndim != 3:
        raise ValueError(f"pressure must be 3D, got shape {p.shape}")
    if t.shape != p.shape or td.shape != p.shape:
        raise ValueError("pressure, temperature, and dewpoint must share a shape")
    return tuple(
        np.asarray(arr, dtype=np.float64)
        for arr in _BACKEND.cape_grid(
            np.asfortranarray(p),
            np.asfortranarray(t),
            np.asfortranarray(td),
        )
    )


def parcel_profile(
    pressure: Union[np.ndarray, xr.DataArray],
    temperature: Union[np.ndarray, xr.DataArray],
    dewpoint: Union[np.ndarray, xr.DataArray],
) -> np.ndarray:
    """Compute the surface-based parcel temperature profile (Celsius)."""
    p = np.asarray(_as_numpy(pressure), dtype=np.float64)
    t = np.asarray(_as_numpy(temperature), dtype=np.float64)
    td = np.asarray(_as_numpy(dewpoint), dtype=np.float64)
    if p.ndim != 1:
        raise ValueError(f"pressure must be 1D, got shape {p.shape}")
    if t.shape != p.shape or td.shape != p.shape:
        raise ValueError("pressure, temperature, and dewpoint must share a shape")
    return np.asarray(
        _BACKEND.parcel_profile(
            np.ascontiguousarray(p),
            np.ascontiguousarray(t),
            np.ascontiguousarray(td),
        ),
        dtype=np.float64,
    )


def parcel_profile_grid(
    pressure: Union[np.ndarray, xr.DataArray],
    temperature: Union[np.ndarray, xr.DataArray],
    dewpoint: Union[np.ndarray, xr.DataArray],
) -> np.ndarray:
    """Compute parcel temperature profiles for a 3D grid (level, lat, lon)."""
    p = np.asarray(_as_numpy(pressure), dtype=np.float64)
    t = np.asarray(_as_numpy(temperature), dtype=np.float64)
    td = np.asarray(_as_numpy(dewpoint), dtype=np.float64)
    if p.ndim != 3:
        raise ValueError(f"pressure must be 3D, got shape {p.shape}")
    if t.shape != p.shape or td.shape != p.shape:
        raise ValueError("pressure, temperature, and dewpoint must share a shape")
    return np.asarray(
        _BACKEND.parcel_profile_grid(
            np.asfortranarray(p),
            np.asfortranarray(t),
            np.asfortranarray(td),
        ),
        dtype=np.float64,
    )


def calculate_cape_cin(
    pressure: Union[np.ndarray, xr.DataArray],
    temperature: Union[np.ndarray, xr.DataArray],
    dewpoint: Union[np.ndarray, xr.DataArray],
):
    """Compute surface-based CAPE and CIN.

    Parameters
    ----------
    pressure : np.ndarray or xr.DataArray
        Pressure in hPa. Either 1D (level,) for a single profile, or
        3D (level, lat, lon) for a spatial grid.
    temperature : np.ndarray or xr.DataArray
        Temperature in Celsius, same shape as ``pressure``.
    dewpoint : np.ndarray or xr.DataArray
        Dewpoint temperature in Celsius, same shape as ``pressure``.

    Returns
    -------
    tuple of (cape, cin, lfc_p, el_p)
        For 1D input: four floats (cape and cin in J/kg, lfc_p and el_p
        in hPa). For 3D input: four 2D arrays (lat, lon), wrapped as
        ``xr.DataArray`` when the input is an xarray object.
    """
    p = _as_numpy(pressure)
    t = _as_numpy(temperature)
    td = _as_numpy(dewpoint)
    is_xarray = hasattr(pressure, "attrs")

    if np.ndim(p) == 1:
        return cape_profile(p, t, td)

    if np.ndim(p) == 3:
        out = cape_grid(p, t, td)
        if is_xarray:
            out = tuple(_wrap_xarray_result(pressure, arr) for arr in out)
        return out

    raise ValueError(f"pressure must be 1D or 3D, got shape {np.shape(p)}")


def calculate_parcel_profile(
    pressure: Union[np.ndarray, xr.DataArray],
    temperature: Union[np.ndarray, xr.DataArray],
    dewpoint: Union[np.ndarray, xr.DataArray],
):
    """Compute the surface-based parcel temperature profile.

    For 1D input returns a 1D array (Celsius); for 3D input returns a
    3D array (level, lat, lon), wrapped as ``xr.DataArray`` when the
    input is an xarray object.
    """
    p = _as_numpy(pressure)
    t = _as_numpy(temperature)
    td = _as_numpy(dewpoint)
    is_xarray = hasattr(pressure, "attrs")

    if np.ndim(p) == 1:
        return parcel_profile(p, t, td)

    if np.ndim(p) == 3:
        out = parcel_profile_grid(p, t, td)
        if is_xarray:
            out = _wrap_xarray_3d(pressure, out)
        return out

    raise ValueError(f"pressure must be 1D or 3D, got shape {np.shape(p)}")


# ---------------------------------------------------------------------------
# Most unstable parcel
# ---------------------------------------------------------------------------
def most_unstable_parcel(
    pressure: Union[np.ndarray, xr.DataArray],
    temperature: Union[np.ndarray, xr.DataArray],
    dewpoint: Union[np.ndarray, xr.DataArray],
    depth: float = 300.0,
):
    """Find the most unstable parcel (maximum theta-e) within the bottom layer.

    Returns ``(mup_p, mup_t, mup_td, mup_idx)`` for a single profile:
    pressure (hPa), temperature (Celsius), dewpoint (Celsius) of the most
    unstable level, and its 0-based index within the cleaned (NaN-free,
    pressure-decreasing) profile. ``depth`` is the layer depth in hPa,
    matching metpy's default of 300 hPa.
    """
    p = np.asarray(_as_numpy(pressure), dtype=np.float64)
    t = np.asarray(_as_numpy(temperature), dtype=np.float64)
    td = np.asarray(_as_numpy(dewpoint), dtype=np.float64)
    if p.ndim != 1:
        raise ValueError(f"pressure must be 1D, got shape {p.shape}")
    if t.shape != p.shape or td.shape != p.shape:
        raise ValueError("pressure, temperature, and dewpoint must share a shape")
    return tuple(
        _BACKEND.most_unstable_parcel(
            np.ascontiguousarray(p),
            np.ascontiguousarray(t),
            np.ascontiguousarray(td),
            float(depth),
        )
    )


def most_unstable_parcel_grid(
    pressure: Union[np.ndarray, xr.DataArray],
    temperature: Union[np.ndarray, xr.DataArray],
    dewpoint: Union[np.ndarray, xr.DataArray],
    depth: float = 300.0,
):
    """Most unstable parcel for every column of a 3D (level, lat, lon) grid.

    Returns four 3D arrays: most-unstable pressure, temperature, dewpoint,
    and level index (NaN / -1 where no unstable level exists). Only the
    level selected for each column is populated; other levels are NaN.
    """
    p = np.asarray(_as_numpy(pressure), dtype=np.float64)
    t = np.asarray(_as_numpy(temperature), dtype=np.float64)
    td = np.asarray(_as_numpy(dewpoint), dtype=np.float64)
    if p.ndim != 3:
        raise ValueError(f"pressure must be 3D, got shape {p.shape}")
    if t.shape != p.shape or td.shape != p.shape:
        raise ValueError("pressure, temperature, and dewpoint must share a shape")
    return tuple(
        np.asarray(arr)
        for arr in _BACKEND.most_unstable_parcel_grid(
            np.asfortranarray(p),
            np.asfortranarray(t),
            np.asfortranarray(td),
            float(depth),
        )
    )


def most_unstable_cape_cin(
    pressure: Union[np.ndarray, xr.DataArray],
    temperature: Union[np.ndarray, xr.DataArray],
    dewpoint: Union[np.ndarray, xr.DataArray],
    depth: float = 300.0,
):
    """Most unstable CAPE/CIN for a single profile (parcel from the max-theta-e level)."""
    p = np.asarray(_as_numpy(pressure), dtype=np.float64)
    t = np.asarray(_as_numpy(temperature), dtype=np.float64)
    td = np.asarray(_as_numpy(dewpoint), dtype=np.float64)
    if p.ndim != 1:
        raise ValueError(f"pressure must be 1D, got shape {p.shape}")
    if t.shape != p.shape or td.shape != p.shape:
        raise ValueError("pressure, temperature, and dewpoint must share a shape")
    return tuple(
        float(v)
        for v in _BACKEND.mucape_profile(
            np.ascontiguousarray(p),
            np.ascontiguousarray(t),
            np.ascontiguousarray(td),
            float(depth),
        )
    )


def most_unstable_cape_cin_grid(
    pressure: Union[np.ndarray, xr.DataArray],
    temperature: Union[np.ndarray, xr.DataArray],
    dewpoint: Union[np.ndarray, xr.DataArray],
    depth: float = 300.0,
):
    """Most unstable CAPE/CIN for a 3D grid; returns four 2D (lat, lon) arrays."""
    p = np.asarray(_as_numpy(pressure), dtype=np.float64)
    t = np.asarray(_as_numpy(temperature), dtype=np.float64)
    td = np.asarray(_as_numpy(dewpoint), dtype=np.float64)
    if p.ndim != 3:
        raise ValueError(f"pressure must be 3D, got shape {p.shape}")
    if t.shape != p.shape or td.shape != p.shape:
        raise ValueError("pressure, temperature, and dewpoint must share a shape")
    return tuple(
        np.asarray(arr, dtype=np.float64)
        for arr in _BACKEND.mucape_grid(
            np.asfortranarray(p),
            np.asfortranarray(t),
            np.asfortranarray(td),
            float(depth),
        )
    )


def calculate_most_unstable_parcel(
    pressure: Union[np.ndarray, xr.DataArray],
    temperature: Union[np.ndarray, xr.DataArray],
    dewpoint: Union[np.ndarray, xr.DataArray],
    depth: float = 300.0,
):
    """Compute the most unstable parcel for 1D or 3D input.

    For 1D input returns ``(mup_p, mup_t, mup_td, mup_idx)`` floats. For 3D
    input returns four 3D arrays, wrapped as ``xr.DataArray`` when the input
    is an xarray object.
    """
    p = _as_numpy(pressure)
    t = _as_numpy(temperature)
    td = _as_numpy(dewpoint)
    is_xarray = hasattr(pressure, "attrs")

    if np.ndim(p) == 1:
        return most_unstable_parcel(p, t, td, depth)

    if np.ndim(p) == 3:
        out = most_unstable_parcel_grid(p, t, td, depth)
        if is_xarray:
            out = tuple(_wrap_xarray_3d(pressure, arr) for arr in out[:3]) + (out[3],)
        return out

    raise ValueError(f"pressure must be 1D or 3D, got shape {np.shape(p)}")


def calculate_most_unstable_cape_cin(
    pressure: Union[np.ndarray, xr.DataArray],
    temperature: Union[np.ndarray, xr.DataArray],
    dewpoint: Union[np.ndarray, xr.DataArray],
    depth: float = 300.0,
):
    """Compute most unstable CAPE/CIN for 1D or 3D input.

    For 1D input returns ``(mucape, mucin, mulfc_p, muel_p)`` floats. For 3D
    input returns four 2D arrays (lat, lon), wrapped as ``xr.DataArray`` when
    the input is an xarray object.
    """
    p = _as_numpy(pressure)
    t = _as_numpy(temperature)
    td = _as_numpy(dewpoint)
    is_xarray = hasattr(pressure, "attrs")

    if np.ndim(p) == 1:
        return most_unstable_cape_cin(p, t, td, depth)

    if np.ndim(p) == 3:
        out = most_unstable_cape_cin_grid(p, t, td, depth)
        if is_xarray:
            out = tuple(_wrap_xarray_result(pressure, arr) for arr in out)
        return out

    raise ValueError(f"pressure must be 1D or 3D, got shape {np.shape(p)}")
