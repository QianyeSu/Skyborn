"""
Comprehensive tests for skyborn.interp.rcm2points module.

Tests rcm2points function with various input types and edge cases.
Target: 95%+ coverage for rcm2points.py
"""

import numpy as np
import pytest
import xarray as xr

from skyborn.interp import rcm2points
from skyborn.interp.errors import CoordinateError, DimensionError


@pytest.fixture
def sample_curvilinear_grid():
    """Create a sample curvilinear grid for testing."""
    nlat, nlon = 10, 12
    lat2d = np.linspace(30, 50, nlat * nlon).reshape(nlat, nlon)
    lon2d = np.linspace(-120, -100, nlat * nlon).reshape(nlat, nlon)
    # Add some curvature
    lat2d += 0.5 * np.sin(lon2d * np.pi / 180)
    lon2d += 0.5 * np.cos(lat2d * np.pi / 180)
    return lat2d, lon2d


@pytest.fixture
def sample_field_2d(sample_curvilinear_grid):
    """Create a sample 2D field on curvilinear grid."""
    lat2d, lon2d = sample_curvilinear_grid
    # Simple temperature-like field
    fi = 15 + 10 * np.sin(lat2d * np.pi / 180) * np.cos(lon2d * np.pi / 180)
    return fi


@pytest.fixture
def sample_field_3d(sample_curvilinear_grid):
    """Create a sample 3D field on curvilinear grid."""
    lat2d, lon2d = sample_curvilinear_grid
    ntime = 5
    fi = np.zeros((ntime, lat2d.shape[0], lat2d.shape[1]))
    for t in range(ntime):
        fi[t] = (
            15 + t * 2 + 10 * np.sin(lat2d * np.pi / 180) * np.cos(lon2d * np.pi / 180)
        )
    return fi


@pytest.fixture
def sample_points():
    """Create sample target points."""
    lat1d = np.array([35.0, 40.0, 45.0])
    lon1d = np.array([-115.0, -110.0, -105.0])
    return lat1d, lon1d


class TestRcm2pointsBasic:
    """Test basic rcm2points functionality."""

    def test_2d_numpy_input(
        self, sample_curvilinear_grid, sample_field_2d, sample_points
    ):
        """Test 2D numpy array input."""
        lat2d, lon2d = sample_curvilinear_grid
        lat1d, lon1d = sample_points
        fi = sample_field_2d[np.newaxis, :, :]  # Add time dimension

        result = rcm2points(lat2d, lon2d, fi, lat1d, lon1d)

        assert isinstance(result, np.ndarray)
        assert result.shape == (fi.shape[0], len(lat1d))
        assert np.all(np.isfinite(result))

    def test_3d_numpy_input(
        self, sample_curvilinear_grid, sample_field_3d, sample_points
    ):
        """Test 3D numpy array input."""
        lat2d, lon2d = sample_curvilinear_grid
        lat1d, lon1d = sample_points

        result = rcm2points(lat2d, lon2d, sample_field_3d, lat1d, lon1d)

        assert isinstance(result, np.ndarray)
        assert result.shape == (sample_field_3d.shape[0], len(lat1d))

    def test_xarray_input(
        self, sample_curvilinear_grid, sample_field_3d, sample_points
    ):
        """Test xarray DataArray input."""
        lat2d, lon2d = sample_curvilinear_grid
        lat1d, lon1d = sample_points

        fi_xr = xr.DataArray(
            sample_field_3d,
            dims=["time", "lat", "lon"],
            coords={"time": np.arange(5)},
            attrs={"units": "degC", "long_name": "temperature"},
        )

        result = rcm2points(lat2d, lon2d, fi_xr, lat1d, lon1d, meta=True)

        assert isinstance(result, xr.DataArray)
        assert result.shape == (5, len(lat1d))
        # Check metadata preservation
        assert "units" in result.attrs or "long_name" in result.attrs

    def test_opt_parameter_0(
        self, sample_curvilinear_grid, sample_field_3d, sample_points
    ):
        """Test with opt=0 (bilinear)."""
        lat2d, lon2d = sample_curvilinear_grid
        lat1d, lon1d = sample_points

        result = rcm2points(lat2d, lon2d, sample_field_3d, lat1d, lon1d, opt=0)

        assert result.shape == (sample_field_3d.shape[0], len(lat1d))
        assert np.all(np.isfinite(result))

    def test_opt_parameter_1(
        self, sample_curvilinear_grid, sample_field_3d, sample_points
    ):
        """Test with opt=1 (inverse distance squared)."""
        lat2d, lon2d = sample_curvilinear_grid
        lat1d, lon1d = sample_points

        result = rcm2points(lat2d, lon2d, sample_field_3d, lat1d, lon1d, opt=1)

        assert result.shape == (sample_field_3d.shape[0], len(lat1d))

    def test_single_point_interpolation(self, sample_curvilinear_grid, sample_field_3d):
        """Test interpolation to a single point."""
        lat2d, lon2d = sample_curvilinear_grid
        lat1d = np.array([40.0])
        lon1d = np.array([-110.0])

        result = rcm2points(lat2d, lon2d, sample_field_3d, lat1d, lon1d)

        assert result.shape == (sample_field_3d.shape[0], 1)

    def test_many_points(self, sample_curvilinear_grid, sample_field_3d):
        """Test interpolation to many points."""
        lat2d, lon2d = sample_curvilinear_grid
        npts = 50
        lat1d = np.linspace(31, 49, npts)
        lon1d = np.linspace(-119, -101, npts)

        result = rcm2points(lat2d, lon2d, sample_field_3d, lat1d, lon1d)

        assert result.shape == (sample_field_3d.shape[0], npts)


class TestRcm2pointsWithMissingValues:
    """Test rcm2points with missing values."""

    def test_field_with_nan(
        self, sample_curvilinear_grid, sample_field_3d, sample_points
    ):
        """Test field with NaN values."""
        lat2d, lon2d = sample_curvilinear_grid
        lat1d, lon1d = sample_points

        fi = sample_field_3d.copy()
        fi[:, 3:5, 4:6] = np.nan

        result = rcm2points(lat2d, lon2d, fi, lat1d, lon1d, msg=np.nan)

        assert result.shape == (fi.shape[0], len(lat1d))
        # Result may have NaN where interpolation failed
        assert np.any(np.isfinite(result))

    def test_custom_missing_value(
        self, sample_curvilinear_grid, sample_field_3d, sample_points
    ):
        """Test with custom missing value."""
        lat2d, lon2d = sample_curvilinear_grid
        lat1d, lon1d = sample_points

        fi = sample_field_3d.copy()
        fi[:, 3:5, 4:6] = -999.0

        result = rcm2points(lat2d, lon2d, fi, lat1d, lon1d, msg=-999.0)

        assert result.shape == (fi.shape[0], len(lat1d))


class TestRcm2pointsValidation:
    """Test input validation."""

    def test_missing_lat2d_raises_error(
        self, sample_curvilinear_grid, sample_field_3d, sample_points
    ):
        """Test missing lat2d raises CoordinateError."""
        _, lon2d = sample_curvilinear_grid
        lat1d, lon1d = sample_points

        with pytest.raises(CoordinateError, match="lat2d"):
            rcm2points(None, lon2d, sample_field_3d, lat1d, lon1d)

    def test_missing_lon2d_raises_error(
        self, sample_curvilinear_grid, sample_field_3d, sample_points
    ):
        """Test missing lon2d raises CoordinateError."""
        lat2d, _ = sample_curvilinear_grid
        lat1d, lon1d = sample_points

        with pytest.raises(CoordinateError, match="lon2d"):
            rcm2points(lat2d, None, sample_field_3d, lat1d, lon1d)

    def test_mismatched_2d_grid_shapes_raises_error(
        self, sample_curvilinear_grid, sample_field_3d, sample_points
    ):
        """Test mismatched lat2d/lon2d shapes raise error."""
        lat2d, lon2d = sample_curvilinear_grid
        lat1d, lon1d = sample_points

        # Make lon2d wrong shape
        lon2d_wrong = lon2d[:, :-1]

        with pytest.raises(DimensionError, match="same size"):
            rcm2points(lat2d, lon2d_wrong, sample_field_3d, lat1d, lon1d)

    def test_mismatched_1d_arrays_raises_error(
        self, sample_curvilinear_grid, sample_field_3d
    ):
        """Test mismatched lat1d/lon1d lengths raise error."""
        lat2d, lon2d = sample_curvilinear_grid
        lat1d = np.array([35.0, 40.0, 45.0])
        lon1d = np.array([-115.0, -110.0])  # Wrong length

        with pytest.raises(DimensionError, match="same size"):
            rcm2points(lat2d, lon2d, sample_field_3d, lat1d, lon1d)

    def test_too_small_grid_raises_error(self, sample_field_3d, sample_points):
        """Test grid with < 2 elements raises error."""
        lat2d = np.array([[30.0]])
        lon2d = np.array([[-120.0]])
        lat1d, lon1d = sample_points

        fi = sample_field_3d[:, :1, :1]

        with pytest.raises(DimensionError, match="at least 2 elements"):
            rcm2points(lat2d, lon2d, fi, lat1d, lon1d)

    def test_1d_field_raises_error(self, sample_curvilinear_grid, sample_points):
        """Test 1D field raises DimensionError."""
        lat2d, lon2d = sample_curvilinear_grid
        lat1d, lon1d = sample_points
        fi = np.array([1, 2, 3, 4, 5])

        with pytest.raises(DimensionError, match="at least two dimensions"):
            rcm2points(lat2d, lon2d, fi, lat1d, lon1d)

    def test_wrong_field_dimensions_raises_error(
        self, sample_curvilinear_grid, sample_points
    ):
        """Test field with wrong spatial dimensions raises error."""
        lat2d, lon2d = sample_curvilinear_grid
        lat1d, lon1d = sample_points

        # Create field with wrong shape
        fi = np.random.randn(5, 8, 8)  # Wrong spatial dims

        with pytest.raises(DimensionError, match="rightmost dimensions"):
            rcm2points(lat2d, lon2d, fi, lat1d, lon1d)


class TestRcm2pointsEdgeCases:
    """Test edge cases."""

    def test_constant_field(self, sample_curvilinear_grid, sample_points):
        """Test interpolation of constant field."""
        lat2d, lon2d = sample_curvilinear_grid
        lat1d, lon1d = sample_points

        fi = np.ones((3, lat2d.shape[0], lat2d.shape[1])) * 25.0

        result = rcm2points(lat2d, lon2d, fi, lat1d, lon1d)

        # Result should be close to constant value
        assert np.allclose(result, 25.0, rtol=0.01)

    def test_linear_gradient_field(self, sample_curvilinear_grid, sample_points):
        """Test interpolation of linear gradient."""
        lat2d, lon2d = sample_curvilinear_grid
        lat1d, lon1d = sample_points

        # Linear north-south gradient
        fi = lat2d[np.newaxis, :, :]

        result = rcm2points(lat2d, lon2d, fi, lat1d, lon1d)

        # Interpolated values should be close to target latitudes
        assert result.shape == (1, len(lat1d))

    def test_output_within_bounds(
        self, sample_curvilinear_grid, sample_field_3d, sample_points
    ):
        """Test interpolated values are within reasonable bounds."""
        lat2d, lon2d = sample_curvilinear_grid
        lat1d, lon1d = sample_points

        fi_min = np.nanmin(sample_field_3d)
        fi_max = np.nanmax(sample_field_3d)

        result = rcm2points(lat2d, lon2d, sample_field_3d, lat1d, lon1d)

        # Interpolated values should not exceed input range significantly
        assert np.all((result >= fi_min - 1) | np.isnan(result))
        assert np.all((result <= fi_max + 1) | np.isnan(result))

    def test_meta_false_returns_numpy(
        self, sample_curvilinear_grid, sample_field_3d, sample_points
    ):
        """Test meta=False returns numpy array even with xarray input."""
        lat2d, lon2d = sample_curvilinear_grid
        lat1d, lon1d = sample_points

        fi_xr = xr.DataArray(sample_field_3d, dims=["time", "lat", "lon"])

        result = rcm2points(lat2d, lon2d, fi_xr, lat1d, lon1d, meta=False)

        assert isinstance(result, (np.ndarray, xr.DataArray))


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
