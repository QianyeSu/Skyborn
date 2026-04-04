"""xarray dataset extraction helpers for Skyborn curly-vector plots."""

from __future__ import annotations

import warnings
from collections.abc import Hashable
from typing import Any

import numpy as np
import xarray as xr

from .._shared.coords import _filled_scalar_field_array

_ISSUED_PLOT_WARNINGS: set[str] = set()


def _warn_plot_once(key: str, message: str) -> None:
    if key in _ISSUED_PLOT_WARNINGS:
        return
    _ISSUED_PLOT_WARNINGS.add(key)
    warnings.warn(message, UserWarning, stacklevel=3)


def _apply_dataset_isel(da: xr.DataArray, isel):
    if not isel:
        return da
    requested = dict(isel)
    indexers = {dim: value for dim, value in requested.items() if dim in da.dims}
    unused = tuple(sorted(dim for dim in requested if dim not in da.dims))
    if unused:
        _warn_plot_once(
            f"unused_isel:{unused}",
            f"Ignoring isel indexers {list(unused)} because they are not present "
            f"in DataArray dims {tuple(da.dims)}.",
        )
    return da.isel(indexers) if indexers else da


def _get_plot_dataarray(ds: xr.Dataset, value, *, isel=None, role: str) -> xr.DataArray:
    if isinstance(value, xr.DataArray):
        da = value
    else:
        da = ds[value]
    da = _apply_dataset_isel(da, isel)
    da = da.squeeze(drop=True)
    if da.ndim == 0:
        raise ValueError(f"{role} must resolve to at least one dimension")
    if da.ndim > 2:
        remaining_dims = [dim for dim, size in da.sizes.items() if size > 1]
        raise ValueError(
            f"{role} must resolve to 1D or 2D data after selection; remaining dims: "
            f"{remaining_dims}. Pass a 2D slice or use isel=..."
        )
    return da


def _transpose_2d_dataarray_to_dims(
    da: xr.DataArray, dims: tuple[Hashable, Hashable], *, role: str
) -> xr.DataArray:
    if da.ndim != 2:
        raise ValueError(f"{role} must resolve to a 2D array before plotting")
    if da.dims == dims:
        return da
    if set(da.dims) != set(dims):
        raise ValueError(
            f"{role} dims {da.dims} do not match the target vector dims {dims}; "
            "2D vector inputs must already live on the same physical grid"
        )
    return da.transpose(*dims)


def _extract_curly_vector_dataset_source(
    ds, x, y, u, v, *, isel=None, return_metadata: bool = False
):
    x_da = _get_plot_dataarray(ds, x, isel=isel, role="x")
    y_da = _get_plot_dataarray(ds, y, isel=isel, role="y")
    u_da = _get_plot_dataarray(ds, u, isel=isel, role="u")
    v_da = _get_plot_dataarray(ds, v, isel=isel, role="v")

    if u_da.ndim != 2 or v_da.ndim != 2:
        raise ValueError("u and v must each resolve to 2D arrays before plotting")
    v_da = _transpose_2d_dataarray_to_dims(v_da, u_da.dims, role="v")
    if u_da.shape != v_da.shape:
        raise ValueError(
            "u and v must share the same 2D shape. If your vector components live on "
            "different staggered or unmatched grids, align them onto the same physical "
            "grid before calling curly_vector()."
        )

    metadata = {
        "vector_dims": tuple(u_da.dims),
        "x_descending": False,
        "y_descending": False,
    }

    if x_da.ndim == 1 and y_da.ndim == 1:
        x_values = np.asarray(x_da.data, dtype=float)
        y_values = np.asarray(y_da.data, dtype=float)
        u_values = np.asarray(u_da.data, dtype=float)
        v_values = np.asarray(v_da.data, dtype=float)

        expected_shape = (y_da.size, x_da.size)
        if u_da.shape != expected_shape:
            raise ValueError(
                f"u/v shape {u_da.shape} does not match the rectilinear x/y grid "
                f"shape {expected_shape}"
            )

        if x_values.size > 1 and x_values[0] > x_values[-1]:
            x_values = x_values[::-1]
            u_values = u_values[:, ::-1]
            v_values = v_values[:, ::-1]
            metadata["x_descending"] = True
        if y_values.size > 1 and y_values[0] > y_values[-1]:
            y_values = y_values[::-1]
            u_values = u_values[::-1, :]
            v_values = v_values[::-1, :]
            metadata["y_descending"] = True

        if return_metadata:
            return x_values, y_values, u_values, v_values, metadata
        return x_values, y_values, u_values, v_values

    if x_da.ndim == 2 and y_da.ndim == 2:
        x_da = _transpose_2d_dataarray_to_dims(x_da, u_da.dims, role="x")
        y_da = _transpose_2d_dataarray_to_dims(y_da, u_da.dims, role="y")
        if x_da.shape != y_da.shape:
            raise ValueError("2D x and y coordinates must have the same shape")
        if x_da.shape != u_da.shape:
            raise ValueError(
                f"2D x/y shape {x_da.shape} must match the u/v shape {u_da.shape}"
            )
    else:
        raise ValueError("x and y must both be 1D or both be 2D")

    values = (
        np.asarray(x_da.data, dtype=float),
        np.asarray(y_da.data, dtype=float),
        np.asarray(u_da.data, dtype=float),
        np.asarray(v_da.data, dtype=float),
    )
    if return_metadata:
        return (*values, metadata)
    return values


def _prepare_dataset_style_field(
    value: Any,
    *,
    isel,
    expected_shape: tuple[int, ...],
    vector_dims: tuple[Hashable, Hashable],
    x_descending: bool,
    y_descending: bool,
    role: str,
) -> Any:
    if value is None or isinstance(value, str) or np.isscalar(value):
        return value

    if isinstance(value, xr.DataArray):
        da = _apply_dataset_isel(value, isel).squeeze(drop=True)
        if da.ndim == 2:
            da = _transpose_2d_dataarray_to_dims(da, vector_dims, role=role)
            array = _filled_scalar_field_array(da.data)
        else:
            array = np.asarray(da)
    else:
        array = np.asarray(value)

    if array.shape != expected_shape:
        return value

    field = _filled_scalar_field_array(array)
    if y_descending:
        field = field[::-1, :]
    if x_descending:
        field = field[:, ::-1]
    return field
