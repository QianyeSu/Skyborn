"""DCAPE calculations backed by the compiled Fortran/C extension."""

from __future__ import annotations

from importlib import import_module
from typing import Union

import numpy as np
import xarray as xr

__all__ = ["calculate_dcape", "dcape_profile", "dcape_grid"]


_BACKEND = import_module(f"{__package__}.dcape_core")


def _as_numpy(values):
    """Extract a NumPy array from xarray-like or array-like inputs."""

    return getattr(values, "values", values)


def dcape_profile(
    pressure: Union[np.ndarray, xr.DataArray],
    temperature: Union[np.ndarray, xr.DataArray],
    dewpoint: Union[np.ndarray, xr.DataArray],
) -> float:
    """Compute DCAPE for a single vertical profile."""

    p = np.asarray(_as_numpy(pressure), dtype=np.float64)
    t = np.asarray(_as_numpy(temperature), dtype=np.float64)
    td = np.asarray(_as_numpy(dewpoint), dtype=np.float64)

    if p.ndim != 1:
        raise ValueError(f"pressure must be 1D, got shape {p.shape}")

    return float(
        _BACKEND.dcape_profile(
            np.asfortranarray(p),
            np.asfortranarray(t),
            np.asfortranarray(td),
        )
    )


def dcape_grid(
    pressure: Union[np.ndarray, xr.DataArray],
    temperature: Union[np.ndarray, xr.DataArray],
    dewpoint: Union[np.ndarray, xr.DataArray],
) -> np.ndarray:
    """Compute DCAPE for a 3D grid of vertical profiles."""

    p = np.asarray(_as_numpy(pressure), dtype=np.float64)
    t = np.asarray(_as_numpy(temperature), dtype=np.float64)
    td = np.asarray(_as_numpy(dewpoint), dtype=np.float64)

    if p.ndim != 3:
        raise ValueError(f"pressure must be 3D, got shape {p.shape}")
    if t.shape != p.shape or td.shape != p.shape:
        raise ValueError("pressure, temperature, and dewpoint must share a shape")

    return np.asarray(
        _BACKEND.dcape_grid(
            np.asfortranarray(p),
            np.asfortranarray(t),
            np.asfortranarray(td),
        ),
        dtype=np.float64,
    )


def _wrap_xarray_result(pressure: xr.DataArray, values: np.ndarray) -> xr.DataArray:
    """Preserve xarray metadata for the DCAPE grid result."""

    dims = list(pressure.dims[-2:])
    coords = {}
    for dim in dims:
        if dim in pressure.coords:
            coords[dim] = pressure.coords[dim]

    return xr.DataArray(
        values,
        dims=dims or ["y", "x"],
        coords=coords,
        attrs={"units": "J kg-1", "long_name": "DCAPE"},
    )


def calculate_dcape(
    pressure: Union[np.ndarray, xr.DataArray],
    temperature: Union[np.ndarray, xr.DataArray],
    dewpoint: Union[np.ndarray, xr.DataArray],
) -> Union[np.ndarray, xr.DataArray, float]:
    """
    Calculate Downdraft Convective Available Potential Energy (DCAPE).

    This implementation uses the compiled Skyborn DCAPE backend. It finds the
    minimum theta-e layer between 700 and 500 hPa, computes the wet-bulb
    temperature along a pseudoadiabat down to the surface, and integrates the
    negative buoyancy.

    Parameters
    ----------
    pressure : np.ndarray or xr.DataArray
        Pressure in hPa. Can be 1D (level,) for a single profile, or 3D
        (level, lat, lon) for a spatial grid.
    temperature : np.ndarray or xr.DataArray
        Temperature in Celsius. Same shape as ``pressure``.
    dewpoint : np.ndarray or xr.DataArray
        Dewpoint temperature in Celsius. Same shape as ``pressure``.

    Returns
    -------
    np.ndarray or xr.DataArray or float
        DCAPE in J kg^-1. For 3D input returns a 2D array (lat, lon).
        For 1D input returns a scalar float.

    Notes
    -----
    - The compiled backend is required; importing this module raises an
      ``ImportError`` if the extension is unavailable.
    - Profiles with fewer than 4 valid levels, or without a detectable
      theta-e minimum between 500--700 hPa, return NaN.
    """

    p = _as_numpy(pressure)
    t = _as_numpy(temperature)
    td = _as_numpy(dewpoint)
    is_xarray = hasattr(pressure, "attrs")

    if np.ndim(p) == 1:
        return float(dcape_profile(p, t, td))

    if np.ndim(p) == 3:
        out = dcape_grid(p, t, td)
        if is_xarray:
            return _wrap_xarray_result(pressure, out)
        return out

    raise ValueError(f"pressure must be 1D or 3D, got shape {np.shape(p)}")
