"""
Comprehensive tests for skyborn.interp.triple_to_grid module.

Tests grid_to_triple and triple_to_grid conversions, validation rules,
chunking behavior, metadata warnings, and backward compatibility wrappers.
Target: 95%+ coverage for triple_to_grid.py
"""

import importlib
import warnings
from types import SimpleNamespace

import dask.array as da
import numpy as np
import pytest
import xarray as xr

from skyborn.interp import grid_to_triple, triple_to_grid
from skyborn.interp.errors import ChunkError, CoordinateError, DimensionError
from skyborn.interp.triple_to_grid import (
    _default_missing_value_for_dtype,
    _grid_to_triple,
    _triple_to_grid,
    _triple_to_grid_2d,
    grid2triple,
    triple2grid,
    triple_to_grid_2d,
)

triple_to_grid_module = importlib.import_module("skyborn.interp.triple_to_grid")


class _FakeCoord:
    """Minimal coordinate stub for xarray-validation branch tests."""

    def __init__(self, ndim, shape):
        self.ndim = ndim
        self.shape = shape


class _FakeGridToTripleDataArray:
    """Minimal DataArray-like stub used to isolate coordinate validation logic."""

    def __init__(self, shape, x_coord, y_coord):
        self.dims = ("y", "x")
        self.shape = shape
        self.ndim = len(shape)
        self.coords = {"x": x_coord, "y": y_coord}
        self.attrs = {}


@pytest.fixture
def small_rect_grid():
    """Create a small 2D rectilinear grid (y, x) with deterministic values."""
    x = np.array([0.0, 1.0, 2.0], dtype=np.float64)  # mx=3
    y = np.array([10.0, 11.0, 12.0, 13.0], dtype=np.float64)  # ny=4
    ny, mx = (y.size, x.size)
    grid = (np.arange(ny)[:, None] * 100 + np.arange(mx)[None, :]).astype(np.float64)
    da = xr.DataArray(
        grid, dims=("y", "x"), coords={"y": y, "x": x}, attrs={"name": "Z"}
    )
    return x, y, da


class TestGridToTripleBasic:
    """Test grid_to_triple -> triple_to_grid round-trip with xarray input."""

    def test_roundtrip_xarray(self, small_rect_grid):
        x_out, y_out, da = small_rect_grid
        triple = grid_to_triple(da)
        assert isinstance(triple, xr.DataArray)
        assert triple.ndim == 2 and triple.shape[0] == 3

        x_in = np.asarray(triple[0].values)
        y_in = np.asarray(triple[1].values)
        z_in = np.asarray(triple[2].values)

        grid = triple_to_grid(
            z_in, x_in=x_in, y_in=y_in, x_out=x_out, y_out=y_out, method=1, domain=1.0
        )
        # Output may be numpy when input is numpy; accept both
        assert isinstance(grid, (np.ndarray, xr.DataArray))

        # Compare as multisets of values to avoid ordering assumptions
        grid_arr = grid.values if isinstance(grid, xr.DataArray) else grid
        np.testing.assert_array_equal(
            np.sort(grid_arr.ravel()), np.sort(da.values.ravel())
        )

    def test_numpy_happy_path(self, small_rect_grid):
        x_out, y_out, da = small_rect_grid
        # numpy data with explicit coords
        triple = grid_to_triple(da.values, x_in=x_out, y_in=y_out)
        assert isinstance(triple, np.ndarray)
        assert triple.shape[0] == 3
        x_in = triple[0]
        y_in = triple[1]
        z_in = triple[2]
        grid = triple_to_grid(
            z_in, x_in=x_in, y_in=y_in, x_out=x_out, y_out=y_out, method=1, domain=1.0
        )
        grid_arr = grid.values if isinstance(grid, xr.DataArray) else grid
        np.testing.assert_array_equal(
            np.sort(grid_arr.ravel()), np.sort(da.values.ravel())
        )

    def test_numpy_nan_missing_input_not_mutated(self, small_rect_grid):
        x_out, y_out, da = small_rect_grid
        data = da.values.copy()
        data[1, 1] = np.nan
        original = data.copy()

        triple = grid_to_triple(data, x_in=x_out, y_in=y_out)

        assert isinstance(triple, np.ndarray)
        np.testing.assert_allclose(data, original, rtol=0.0, atol=0.0, equal_nan=True)
        assert triple.shape[1] == data.size - 1

    def test_numpy_custom_missing_input_not_mutated(self, small_rect_grid):
        x_out, y_out, da = small_rect_grid
        msg = np.float64(-9999.0)
        data = da.values.copy()
        data[2, 1] = msg
        original = data.copy()

        triple = grid_to_triple(data, x_in=x_out, y_in=y_out, msg_py=msg)

        assert isinstance(triple, np.ndarray)
        np.testing.assert_array_equal(data, original)
        assert triple.shape[1] == data.size - 1


class TestTripleToGridBasic:
    """Test triple_to_grid with multiple leftmost dims and correctness."""

    def test_multi_left_dims(self, small_rect_grid):
        x_out, y_out, da = small_rect_grid
        X, Y = np.meshgrid(x_out, y_out, indexing="xy")
        x_in = X.ravel()
        y_in = Y.ravel()
        z_base = da.values.ravel()

        data = np.stack([z_base, z_base + 10.0], axis=0)  # (time=2, points)
        data_da = xr.DataArray(data, dims=("time", "points"))

        grid = triple_to_grid(
            data_da,
            x_in=x_in,
            y_in=y_in,
            x_out=x_out,
            y_out=y_out,
            method=1,
            domain=1.0,
        )
        assert isinstance(grid, xr.DataArray)
        assert grid.shape == (2, y_out.size, x_out.size)

        # Compare only multisets of values for each slice
        np.testing.assert_array_equal(
            np.sort(grid.values[0].ravel()), np.sort(da.values.ravel())
        )
        np.testing.assert_array_equal(
            np.sort(grid.values[1].ravel()), np.sort((da.values + 10.0).ravel())
        )

    def test_method0_happy_path(self, small_rect_grid):
        x_out, y_out, da = small_rect_grid
        X, Y = np.meshgrid(x_out, y_out, indexing="xy")
        x_in = X.ravel()
        y_in = Y.ravel()
        z = da.values.ravel()
        grid = triple_to_grid(
            z, x_in=x_in, y_in=y_in, x_out=x_out, y_out=y_out, method=0, domain=1.0
        )
        grid_arr = grid.values if isinstance(grid, xr.DataArray) else grid
        np.testing.assert_array_equal(
            np.sort(grid_arr.ravel()), np.sort(da.values.ravel())
        )

    def test_numpy_nan_missing_input_not_mutated(self, small_rect_grid):
        x_out, y_out, da = small_rect_grid
        X, Y = np.meshgrid(x_out, y_out, indexing="xy")
        x_in = X.ravel()
        y_in = Y.ravel()
        z = da.values.ravel().copy()
        z[3] = np.nan
        original = z.copy()

        grid = triple_to_grid(
            z, x_in=x_in, y_in=y_in, x_out=x_out, y_out=y_out, method=1, domain=1.0
        )

        assert isinstance(grid, np.ndarray)
        np.testing.assert_allclose(z, original, rtol=0.0, atol=0.0, equal_nan=True)

    def test_numpy_custom_missing_input_not_mutated(self, small_rect_grid):
        x_out, y_out, da = small_rect_grid
        X, Y = np.meshgrid(x_out, y_out, indexing="xy")
        x_in = X.ravel()
        y_in = Y.ravel()
        msg = np.float64(-9999.0)
        z = da.values.ravel().copy()
        z[5] = msg
        original = z.copy()

        grid = triple_to_grid(
            z,
            x_in=x_in,
            y_in=y_in,
            x_out=x_out,
            y_out=y_out,
            method=1,
            domain=1.0,
            missing_value=msg,
        )

        assert isinstance(grid, np.ndarray)
        np.testing.assert_array_equal(z, original)

    def test_ndim_check_after_len_equality(self, small_rect_grid):
        """Craft x_in as 2D with shape (N,1) so length equality passes but ndim>1 triggers."""
        x_out, y_out, da = small_rect_grid
        X, Y = np.meshgrid(x_out, y_out, indexing="xy")
        x_in = X.ravel().reshape(-1, 1)  # shape (N,1)
        y_in = Y.ravel()  # shape (N,)
        z = da.values.ravel()
        with pytest.raises(DimensionError):
            triple_to_grid(z, x_in=x_in, y_in=y_in, x_out=x_out, y_out=y_out)


class TestPrivateWrapperBranches:
    """Cover dtype and dask-specific wrapper branches."""

    def test_default_missing_value_for_dtype_handles_complex_and_integer(self):
        complex_missing = _default_missing_value_for_dtype(np.complex128)
        integer_missing = _default_missing_value_for_dtype(np.int32)

        assert np.isnan(complex_missing.real)
        assert np.isnan(complex_missing.imag)
        assert integer_missing == triple_to_grid_module.msg_dtype[np.int32]

    def test_grid_to_triple_integer_wrapper_uses_non_nan_missing_path(
        self, small_rect_grid
    ):
        x_out, y_out, grid_da = small_rect_grid
        msg = np.int32(-9999)
        data = grid_da.values.astype(np.int32)
        data[1, 1] = msg
        original = data.copy()

        triple = _grid_to_triple(x_out, y_out, data, msg_py=None)

        assert triple.shape[1] == data.size - 1
        np.testing.assert_array_equal(data, original)

    def test_triple_to_grid_integer_wrapper_uses_non_nan_missing_path(
        self, small_rect_grid
    ):
        x_out, y_out, grid_da = small_rect_grid
        X, Y = np.meshgrid(x_out, y_out, indexing="xy")
        x_in = X.ravel()
        y_in = Y.ravel()
        msg = np.int32(-9999)
        values = grid_da.values.astype(np.int32).ravel()
        values[4] = msg
        original = values.copy()

        grid = _triple_to_grid(
            values,
            x_in=x_in,
            y_in=y_in,
            x_out=x_out,
            y_out=y_out,
            shape=(y_out.size, x_out.size),
            method=1,
            domain=1.0,
            msg_py=msg,
        )

        assert grid.shape == (y_out.size, x_out.size)
        np.testing.assert_array_equal(values, original)

    def test_triple_to_grid_chunked_multidim_warns_and_preserves_dask(
        self, small_rect_grid
    ):
        x_out, y_out, grid_da = small_rect_grid
        X, Y = np.meshgrid(x_out, y_out, indexing="xy")
        x_in = X.ravel()
        y_in = Y.ravel()
        points = x_in.size

        data_da = xr.DataArray(
            np.stack([grid_da.values.ravel(), grid_da.values.ravel() + 100.0], axis=0),
            dims=("time", "points"),
        ).chunk({"time": 1, "points": points})

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            grid = triple_to_grid(
                data_da,
                x_in,
                y_in,
                x_out,
                y_out,
                method=1,
                meta=True,
            )

        assert any("metadata is not yet supported" in str(w.message) for w in caught)
        assert isinstance(grid, xr.DataArray)
        assert hasattr(grid.data, "chunks")
        assert grid.shape == (2, y_out.size, x_out.size)

    def test_triple_to_grid_dask_array_input_materializes_numpy(self, small_rect_grid):
        x_out, y_out, grid_da = small_rect_grid
        X, Y = np.meshgrid(x_out, y_out, indexing="xy")
        x_in = X.ravel()
        y_in = Y.ravel()
        values = da.from_array(grid_da.values.ravel(), chunks=(x_in.size,))

        grid = triple_to_grid(
            values,
            x_in=x_in,
            y_in=y_in,
            x_out=x_out,
            y_out=y_out,
            method=1,
            domain=1.0,
        )

        assert isinstance(grid, np.ndarray)
        assert grid.shape == (y_out.size, x_out.size)


class TestValidation:
    """Validation rules and error messages for triple_to_grid."""

    def test_required_coords_for_numpy(self, small_rect_grid):
        x_out, y_out, da = small_rect_grid
        X, Y = np.meshgrid(x_out, y_out, indexing="xy")
        x_in = X.ravel()
        y_in = Y.ravel()
        z = da.values.ravel()
        with pytest.raises(CoordinateError):
            triple_to_grid(z, None, y_in, x_out, y_out)
        with pytest.raises(CoordinateError):
            triple_to_grid(z, x_in, None, x_out, y_out)

    def test_ndim_and_length_checks(self, small_rect_grid):
        x_out, y_out, da = small_rect_grid
        X, Y = np.meshgrid(x_out, y_out, indexing="xy")
        x_in = X.ravel()
        y_in = Y.ravel()
        z = da.values.ravel()
        with pytest.raises(DimensionError):
            triple_to_grid(z, x_in.reshape(1, -1), y_in, x_out, y_out)
        with pytest.raises(DimensionError):
            triple_to_grid(z, x_in, y_in.reshape(1, -1), x_out, y_out)
        with pytest.raises(DimensionError):
            triple_to_grid(z, x_in, y_in, x_out.reshape(1, -1), y_out)
        with pytest.raises(DimensionError):
            triple_to_grid(z, x_in, y_in, x_out, y_out.reshape(1, -1))
        with pytest.raises(DimensionError):
            triple_to_grid(z[:-1], x_in, y_in, x_out, y_out)

    def test_method_domain_distmx_rules(self, small_rect_grid):
        x_out, y_out, da = small_rect_grid
        X, Y = np.meshgrid(x_out, y_out, indexing="xy")
        x_in = X.ravel()
        y_in = Y.ravel()
        z = da.values.ravel()
        with pytest.raises(TypeError):
            triple_to_grid(z, x_in, y_in, x_out, y_out, method="1")
        with pytest.raises(TypeError):
            triple_to_grid(z, x_in, y_in, x_out, y_out, method=2)
        with pytest.raises(ValueError):
            triple_to_grid(z, x_in, y_in, x_out, y_out, domain=np.array([1.0, 2.0]))
        with pytest.raises(ValueError):
            triple_to_grid(z, x_in, y_in, x_out, y_out, method=0, distmx=10.0)
        with pytest.raises(ValueError):
            triple_to_grid(z, x_in, y_in, x_out, y_out, method=1, distmx=[10.0, 20.0])


class TestChunking:
    """Chunking behavior and errors for triple_to_grid."""

    def test_wrong_last_dim_chunk_raises(self, small_rect_grid):
        x_out, y_out, da = small_rect_grid
        X, Y = np.meshgrid(x_out, y_out, indexing="xy")
        x_in = X.ravel()
        y_in = Y.ravel()
        z = da.values.ravel()
        # Wrong last-dim chunk -> should raise ChunkError
        data_da = xr.DataArray(z, dims=("points",)).chunk(
            {"points": max(1, z.size // 2)}
        )
        with pytest.raises(ChunkError):
            triple_to_grid(data_da, x_in, y_in, x_out, y_out, method=1)


class TestMetaAndWrappers:
    """Metadata warnings and backward compatible wrappers."""

    def test_meta_warning(self, small_rect_grid):
        x_out, y_out, da = small_rect_grid
        X, Y = np.meshgrid(x_out, y_out, indexing="xy")
        x_in = X.ravel()
        y_in = Y.ravel()
        z = da.values.ravel()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            grid = triple_to_grid(z, x_in, y_in, x_out, y_out, method=1, meta=True)
            assert any("metadata is not yet supported" in str(ww.message) for ww in w)
            assert isinstance(grid, (np.ndarray, xr.DataArray))

    def test_grid2triple_triple2grid_compat(self, small_rect_grid):
        x_out, y_out, da = small_rect_grid
        with warnings.catch_warnings(record=True) as w1:
            warnings.simplefilter("always")
            trip = grid2triple(da.coords["x"], da.coords["y"], da, None)
            assert any("deprecated" in str(ww.message).lower() for ww in w1)
        assert isinstance(trip, xr.DataArray)
        assert trip.ndim == 2 and trip.shape[0] == 3
        with warnings.catch_warnings(record=True) as w2:
            warnings.simplefilter("always")
            grid = triple2grid(
                trip[0].values, trip[1].values, trip[2].values, x_out, y_out
            )
            assert any("deprecated" in str(ww.message).lower() for ww in w2)
        grid_arr = grid.values if isinstance(grid, xr.DataArray) else grid
        # Compare as multisets of values only
        np.testing.assert_array_equal(
            np.sort(grid_arr.ravel()), np.sort(da.values.ravel())
        )

    def test_triple_to_grid_2d_warns(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            triple_to_grid_2d(None, None, None, None, None, None)
            assert any("not yet implemented" in str(ww.message).lower() for ww in w)

    def test_private_stub_triple_to_grid_2d_noop(self):
        # Call the private stub to cover its placeholder line
        _triple_to_grid_2d(None, None, None, None, None, None)


class TestGridToTripleValidation:
    """Validation errors for grid_to_triple."""

    def test_numpy_without_coords_raises(self, small_rect_grid):
        x, y, da = small_rect_grid
        with pytest.raises(CoordinateError):
            # numpy input requires x_in/y_in
            grid_to_triple(da.values)

    def test_wrong_dims_and_sizes(self, small_rect_grid):
        x, y, da = small_rect_grid
        # data must be 2D
        with pytest.raises(DimensionError):
            grid_to_triple(da.values[np.newaxis, ...], x_in=x, y_in=y)
        # x_in must be 1D
        with pytest.raises(DimensionError):
            grid_to_triple(da.values, x_in=x.reshape(1, -1), y_in=y)
        # y_in must be 1D
        with pytest.raises(DimensionError):
            grid_to_triple(da.values, x_in=x, y_in=y.reshape(1, -1))
        # x_in size must match right dim
        with pytest.raises(DimensionError):
            grid_to_triple(da.values, x_in=x[:-1], y_in=y)
        # y_in size must match left dim
        with pytest.raises(DimensionError):
            grid_to_triple(da.values, x_in=x, y_in=y[:-1])

    def test_xarray_data_not_2d_raises(self, small_rect_grid):
        x, y, da = small_rect_grid
        da3 = da.expand_dims({"z": 2})  # make it 3D
        with pytest.raises(DimensionError):
            grid_to_triple(da3)

    def test_xarray_x_coord_must_be_1d(self, monkeypatch):
        """Test xarray-style inputs reject multi-dimensional x coordinates."""
        monkeypatch.setattr(
            triple_to_grid_module,
            "xr",
            SimpleNamespace(DataArray=_FakeGridToTripleDataArray),
        )
        data = _FakeGridToTripleDataArray(
            shape=(4, 3),
            x_coord=_FakeCoord(ndim=2, shape=(4, 3)),
            y_coord=_FakeCoord(ndim=1, shape=(4,)),
        )

        with pytest.raises(DimensionError, match="`x_in` must have one dimension"):
            grid_to_triple(data)

    def test_xarray_x_coord_size_must_match_right_dimension(self, monkeypatch):
        """Test xarray-style inputs reject x coordinates with the wrong length."""
        monkeypatch.setattr(
            triple_to_grid_module,
            "xr",
            SimpleNamespace(DataArray=_FakeGridToTripleDataArray),
        )
        data = _FakeGridToTripleDataArray(
            shape=(4, 3),
            x_coord=_FakeCoord(ndim=1, shape=(2,)),
            y_coord=_FakeCoord(ndim=1, shape=(4,)),
        )

        with pytest.raises(
            DimensionError,
            match="`x_in` must have the same size",
        ):
            grid_to_triple(data)

    def test_xarray_y_coord_must_be_1d(self, monkeypatch):
        """Test xarray-style inputs reject multi-dimensional y coordinates."""
        monkeypatch.setattr(
            triple_to_grid_module,
            "xr",
            SimpleNamespace(DataArray=_FakeGridToTripleDataArray),
        )
        data = _FakeGridToTripleDataArray(
            shape=(4, 3),
            x_coord=_FakeCoord(ndim=1, shape=(3,)),
            y_coord=_FakeCoord(ndim=2, shape=(4, 1)),
        )

        with pytest.raises(DimensionError, match="`y_in` must have one dimension"):
            grid_to_triple(data)

    def test_xarray_y_coord_size_must_match_left_dimension(self, monkeypatch):
        """Test xarray-style inputs reject y coordinates with the wrong length."""
        monkeypatch.setattr(
            triple_to_grid_module,
            "xr",
            SimpleNamespace(DataArray=_FakeGridToTripleDataArray),
        )
        data = _FakeGridToTripleDataArray(
            shape=(4, 3),
            x_coord=_FakeCoord(ndim=1, shape=(3,)),
            y_coord=_FakeCoord(ndim=1, shape=(3,)),
        )

        with pytest.raises(
            DimensionError,
            match="`y_in` must have the same size",
        ):
            grid_to_triple(data)


class TestTripleToGridDaskPath:
    """Exercise the dask chunked input branch for triple_to_grid."""

    def test_non_chunked_xarray_path_returns_materialized_array(self, small_rect_grid):
        x_out, y_out, da = small_rect_grid
        X, Y = np.meshgrid(x_out, y_out, indexing="xy")
        x_in = X.ravel()
        y_in = Y.ravel()
        data_da = xr.DataArray(da.values.ravel(), dims=("points",))

        grid = triple_to_grid(data_da, x_in, y_in, x_out, y_out, method=1)

        assert isinstance(grid, xr.DataArray)
        assert not hasattr(grid.data, "chunks")

    def test_chunked_ok_path_returns_xarray(self, small_rect_grid):
        x_out, y_out, da = small_rect_grid
        # Build triple inputs covering full grid
        X, Y = np.meshgrid(x_out, y_out, indexing="xy")
        x_in = X.ravel()
        y_in = Y.ravel()
        z = da.values.ravel()

        # Prepare chunked DataArray with correct last-dim chunk (equal to points)
        points = z.size
        data_da = xr.DataArray(z, dims=("points",)).chunk({"points": points})

        grid = triple_to_grid(data_da, x_in, y_in, x_out, y_out, method=1)
        # When input is xarray and chunked, implementation returns xr.DataArray without compute
        assert isinstance(grid, xr.DataArray)
        # Heuristic: dask Array has attribute chunks
        assert hasattr(grid.data, "chunks")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
