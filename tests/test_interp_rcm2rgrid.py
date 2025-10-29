"""
Comprehensive tests for skyborn.interp.rcm2rgrid module.

Tests rcm2rgrid and rgrid2rcm bidirectional interpolation with numpy/xarray
inputs, validation errors, chunking checks, and round-trip consistency.
Target: 95%+ coverage for rcm2rgrid.py
"""

import numpy as np
import pytest
import xarray as xr

from skyborn.interp import rcm2rgrid, rgrid2rcm
from skyborn.interp.errors import ChunkError, CoordinateError


@pytest.fixture
def rect_grid():
    """Create a small rectilinear grid and its 2D mesh for testing."""
    lat1d = np.linspace(-10.0, 10.0, 5, dtype=np.float64)  # ny=5
    lon1d = np.linspace(100.0, 120.0, 6, dtype=np.float64)  # nx=6
    lat2d, lon2d = np.meshgrid(lat1d, lon1d, indexing="ij")  # (ny, nx)
    return lat1d, lon1d, lat2d, lon2d


@pytest.fixture
def curv_field(rect_grid):
    """Create a deterministic 3D field (time, lat, lon) on the rect grid."""
    lat1d, lon1d, lat2d, lon2d = rect_grid
    ny, nx = lat2d.shape
    time = np.arange(2)
    base = (np.arange(ny)[:, None] * 10 + np.arange(nx)[None, :]).astype(np.float64)
    fi = np.stack([base + t * 1000 for t in time], axis=0)  # (time, ny, nx)
    da = xr.DataArray(
        fi,
        dims=("time", "lat", "lon"),
        coords={"time": time, "lat": lat1d, "lon": lon1d},
        attrs={"units": "arb"},
    )
    return da


class TestRcm2rgridBasic:
    """Basic behavior and identity cases for rcm2rgrid."""

    def test_identity_xarray(self, rect_grid, curv_field):
        lat1d, lon1d, lat2d, lon2d = rect_grid
        fi = curv_field
        fo = rcm2rgrid(lat2d, lon2d, fi, lat1d, lon1d)

        assert isinstance(fo, xr.DataArray)
        assert fo.dims == fi.dims
        assert fo.attrs == fi.attrs
        assert np.allclose(fo.coords[fi.dims[-2]], lat1d)
        assert np.allclose(fo.coords[fi.dims[-1]], lon1d)
        np.testing.assert_allclose(fo.values, fi.values, rtol=0, atol=0)

    def test_identity_numpy(self, rect_grid, curv_field):
        lat1d, lon1d, lat2d, lon2d = rect_grid
        fi_np = curv_field.values
        fo_np = rcm2rgrid(lat2d, lon2d, fi_np, lat1d, lon1d)
        assert isinstance(fo_np, np.ndarray)
        np.testing.assert_allclose(fo_np, fi_np, rtol=0, atol=0)

    def test_chunk_validation_raises(self, rect_grid, curv_field):
        lat1d, lon1d, lat2d, lon2d = rect_grid
        fi = curv_field
        with pytest.raises(ChunkError):
            rcm2rgrid(lat2d, lon2d, fi.chunk({"lat": 2, "lon": 3}), lat1d, lon1d)

    def test_missing_coords_raises(self, curv_field, rect_grid):
        lat1d, lon1d, lat2d, lon2d = rect_grid
        fi = curv_field
        with pytest.raises(CoordinateError):
            rcm2rgrid(None, lon2d, fi, lat1d, lon1d)
        with pytest.raises(CoordinateError):
            rcm2rgrid(lat2d, None, fi, lat1d, lon1d)


class TestRgrid2rcmBasic:
    """Basic behavior and identity cases for rgrid2rcm."""

    def test_identity_xarray(self, rect_grid):
        lat1d, lon1d, lat2d, lon2d = rect_grid
        ny, nx = lat2d.shape
        time = np.arange(2)
        base = (np.arange(ny)[:, None] * 10 + np.arange(nx)[None, :]).astype(np.float64)
        fi_rect = np.stack([base + t * 1000 for t in time], axis=0)
        da_rect = xr.DataArray(
            fi_rect,
            dims=("time", "lat", "lon"),
            coords={"time": time, "lat": lat1d, "lon": lon1d},
            attrs={"units": "arb"},
        )
        fo_curv = rgrid2rcm(lat1d, lon1d, da_rect, lat2d, lon2d)
        assert isinstance(fo_curv, xr.DataArray)
        assert fo_curv.dims == da_rect.dims
        np.testing.assert_allclose(fo_curv.values, da_rect.values, rtol=0, atol=0)

    def test_identity_numpy(self, rect_grid):
        lat1d, lon1d, lat2d, lon2d = rect_grid
        ny, nx = lat2d.shape
        time = np.arange(2)
        base = (np.arange(ny)[:, None] * 10 + np.arange(nx)[None, :]).astype(np.float64)
        fi_rect = np.stack([base + t * 1000 for t in time], axis=0)
        fo_curv = rgrid2rcm(lat1d, lon1d, fi_rect, lat2d, lon2d)
        assert isinstance(fo_curv, np.ndarray)
        np.testing.assert_allclose(fo_curv, fi_rect, rtol=0, atol=0)


class TestRgrid2rcmValidation:
    """Validation and error handling for rgrid2rcm."""

    def test_numpy_missing_latlon_raises(self, rect_grid):
        lat1d, lon1d, lat2d, lon2d = rect_grid
        fi = np.zeros((1, lat2d.shape[0], lon2d.shape[1]), dtype=np.float64)
        with pytest.raises(CoordinateError):
            rgrid2rcm(None, lon1d, fi, lat2d, lon2d)
        with pytest.raises(CoordinateError):
            rgrid2rcm(lat1d, None, fi, lat2d, lon2d)

    def test_chunk_validation_raises(self, rect_grid):
        lat1d, lon1d, lat2d, lon2d = rect_grid
        ny, nx = lat2d.shape
        da = xr.DataArray(
            np.zeros((1, ny, nx), dtype=np.float64),
            dims=("time", "lat", "lon"),
            coords={"lat": lat1d, "lon": lon1d},
        )
        with pytest.raises(Exception):
            rgrid2rcm(lat1d, lon1d, da.chunk({"lat": 2, "lon": 3}), lat2d, lon2d)


class TestRoundTrip:
    """Round-trip consistency rgrid2rcm -> rcm2rgrid."""

    def test_rect_to_curv_to_rect(self, rect_grid):
        lat1d, lon1d, lat2d, lon2d = rect_grid
        ny, nx = lat2d.shape
        base = (np.arange(ny)[:, None] * 10 + np.arange(nx)[None, :]).astype(np.float64)
        fi = xr.DataArray(
            base, dims=("lat", "lon"), coords={"lat": lat1d, "lon": lon1d}
        )
        fi3 = fi.expand_dims({"time": [0]})  # Add singleton time dim
        curv = rgrid2rcm(lat1d, lon1d, fi3, lat2d, lon2d)
        rect = rcm2rgrid(lat2d, lon2d, curv, lat1d, lon1d)
        np.testing.assert_allclose(rect.values, fi3.values, rtol=0, atol=0)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
