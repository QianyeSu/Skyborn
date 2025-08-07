"""
Test the gridfill xarray interface.

This module provides comprehensive tests for the xarray interface to gridfill,
ensuring proper handling of coordinate metadata, missing value detection,
and attribute preservation.
"""

import pytest
import numpy as np
import numpy.ma as ma

try:
    import xarray as xr

    HAS_XARRAY = True
except ImportError:
    HAS_XARRAY = False

from skyborn.gridfill.xarray import (
    fill,
    fill_multiple,
    validate_grid_coverage,
    _find_spatial_coordinates,
    _detect_cyclic_longitude,
)


@pytest.mark.skipif(not HAS_XARRAY, reason="xarray not available")
class TestXArrayInterface:
    """Test xarray interface for gridfill."""

    @pytest.fixture
    def sample_data(self):
        """Create sample xarray DataArray with missing values."""
        # Create coordinate arrays
        lons = np.linspace(0, 360, 36, endpoint=False)
        lats = np.linspace(-80, 80, 16)
        time = np.arange(3)

        # Create data array
        np.random.seed(42)  # For reproducible tests
        data = np.random.rand(3, 16, 36)

        # Add some missing values
        data[0, 5:8, 10:15] = np.nan
        data[1, 12:14, 20:25] = np.nan
        data[2, 2:5, 30:35] = np.nan

        # Create DataArray with proper coordinates
        da = xr.DataArray(
            data,
            coords={
                "time": time,
                "lat": (
                    "lat",
                    lats,
                    {"standard_name": "latitude", "units": "degrees_north"},
                ),
                "lon": (
                    "lon",
                    lons,
                    {"standard_name": "longitude", "units": "degrees_east"},
                ),
            },
            dims=["time", "lat", "lon"],
            attrs={
                "units": "K",
                "long_name": "temperature",
                "standard_name": "air_temperature",
            },
        )

        return da

    @pytest.fixture
    def global_data(self):
        """Create global dataset with cyclic longitude."""
        lons = np.linspace(0, 360, 72, endpoint=False)
        lats = np.linspace(-90, 90, 36)

        np.random.seed(123)
        data = np.random.rand(36, 72)

        # Add missing values
        data[15:20, 30:40] = np.nan

        da = xr.DataArray(
            data,
            coords={
                "lat": (
                    "lat",
                    lats,
                    {"standard_name": "latitude", "units": "degrees_north"},
                ),
                "lon": (
                    "lon",
                    lons,
                    {"standard_name": "longitude", "units": "degrees_east"},
                ),
            },
            dims=["lat", "lon"],
            attrs={"units": "C", "long_name": "sea_surface_temperature"},
        )

        return da

    def test_coordinate_detection(self, sample_data):
        """Test automatic spatial coordinate detection."""
        y_name, x_name, y_dim, x_dim = _find_spatial_coordinates(sample_data)

        assert y_name == "lat"
        assert x_name == "lon"
        assert y_dim == 1
        assert x_dim == 2

    def test_cyclic_detection(self, global_data):
        """Test automatic cyclic boundary detection."""
        lon_coord = global_data.coords["lon"]
        cyclic = _detect_cyclic_longitude(lon_coord)
        assert cyclic is True

    def test_fill_basic(self, sample_data):
        """Test basic fill operation."""
        filled = fill(sample_data, eps=1e-3)

        # Check that result is DataArray
        assert isinstance(filled, xr.DataArray)

        # Check that no NaN values remain
        assert not np.any(np.isnan(filled.values))

        # Check that coordinates are preserved
        assert filled.coords.keys() == sample_data.coords.keys()

        # Check that original attributes are preserved (excluding added history)
        original_attrs = {k: v for k, v in sample_data.attrs.items()}
        filled_attrs = {k: v for k, v in filled.attrs.items() if k != "history"}
        assert filled_attrs == original_attrs

    def test_fill_with_explicit_dimensions(self, sample_data):
        """Test fill with explicitly specified dimensions."""
        filled = fill(sample_data, eps=1e-3, x_dim="lon", y_dim="lat")

        assert isinstance(filled, xr.DataArray)
        assert not np.any(np.isnan(filled.values))

    def test_fill_cyclic(self, global_data):
        """Test fill with cyclic boundary conditions."""
        filled = fill(global_data, eps=1e-3, cyclic=True)

        assert isinstance(filled, xr.DataArray)
        assert not np.any(np.isnan(filled.values))

    def test_fill_no_missing_values(self, sample_data):
        """Test behavior when no missing values present."""
        # Create data without missing values
        clean_data = sample_data.fillna(0)

        with pytest.warns(UserWarning, match="No missing values found"):
            filled = fill(clean_data, eps=1e-3)

        # Should return essentially the same data
        xr.testing.assert_allclose(filled, clean_data)

    def test_fill_attributes_disabled(self, sample_data):
        """Test fill with keep_attrs=False."""
        filled = fill(sample_data, eps=1e-3, keep_attrs=False)

        # Should not have original attributes (except possibly history)
        assert "units" not in filled.attrs
        assert "long_name" not in filled.attrs

    def test_fill_multiple(self, sample_data):
        """Test fill_multiple function."""
        # Create second dataset
        data2 = sample_data.copy()
        data2.values[0, 8:12, 15:20] = np.nan

        filled_list = fill_multiple([sample_data, data2], eps=1e-3)

        assert len(filled_list) == 2
        assert all(isinstance(da, xr.DataArray) for da in filled_list)
        assert all(not np.any(np.isnan(da.values)) for da in filled_list)

    def test_validation_good_data(self, sample_data):
        """Test validation with good quality data."""
        validation = validate_grid_coverage(sample_data, min_coverage=0.1)

        assert validation["valid"] is True
        assert validation["coverage"] > 0.8
        assert len(validation["messages"]) == 0

    def test_validation_sparse_data(self):
        """Test validation with sparse data."""
        # Create very sparse data
        lons = np.linspace(0, 360, 10, endpoint=False)
        lats = np.linspace(-80, 80, 8)

        data = np.full((8, 10), np.nan)
        data[2:4, 3:6] = 1.0  # Only small region has data

        da = xr.DataArray(
            data,
            coords={
                "lat": ("lat", lats, {"standard_name": "latitude"}),
                "lon": ("lon", lons, {"standard_name": "longitude"}),
            },
            dims=["lat", "lon"],
        )

        validation = validate_grid_coverage(da, min_coverage=0.5)

        assert validation["valid"] is False
        assert validation["coverage"] < 0.5
        assert len(validation["messages"]) > 0

    def test_invalid_input_type(self):
        """Test error handling for invalid input types."""
        with pytest.raises(TypeError, match="data must be xarray.DataArray"):
            fill(np.random.rand(10, 10), eps=1e-3)

    def test_invalid_dimensions(self, sample_data):
        """Test error handling for invalid dimension names."""
        with pytest.raises(ValueError, match="x_dim 'invalid' not found"):
            fill(sample_data, eps=1e-3, x_dim="invalid", y_dim="lat")

    def test_missing_coordinates(self):
        """Test error handling when spatial coordinates cannot be found."""
        # Create DataArray without recognizable spatial coordinates
        data = xr.DataArray(np.random.rand(10, 15), dims=["dim1", "dim2"])

        with pytest.raises(ValueError, match="Could not identify spatial coordinates"):
            fill(data, eps=1e-3)

    def test_coordinate_name_variations(self):
        """Test detection of various coordinate naming conventions."""
        lons = np.linspace(0, 360, 10, endpoint=False)
        lats = np.linspace(-80, 80, 8)

        # Test different coordinate naming conventions
        naming_variations = [
            ({"latitude": lats, "longitude": lons}, ["latitude", "longitude"]),
            ({"lat": lats, "lon": lons}, ["lat", "lon"]),
            ({"y": lats, "x": lons}, ["y", "x"]),
        ]

        for coords_dict, expected_names in naming_variations:
            data = np.random.rand(8, 10)
            data[2:4, 3:6] = np.nan

            da = xr.DataArray(data, coords=coords_dict, dims=list(coords_dict.keys()))

            y_name, x_name, _, _ = _find_spatial_coordinates(da)
            assert y_name == expected_names[0]
            assert x_name == expected_names[1]

    def test_masked_array_input(self):
        """Test handling of masked array input."""
        lons = np.linspace(0, 360, 20, endpoint=False)
        lats = np.linspace(-80, 80, 10)

        # Create masked array
        data = np.random.rand(10, 20)
        mask = np.zeros_like(data, dtype=bool)
        mask[3:6, 8:12] = True
        masked_data = ma.array(data, mask=mask)

        da = xr.DataArray(
            masked_data,
            coords={
                "lat": ("lat", lats, {"standard_name": "latitude"}),
                "lon": ("lon", lons, {"standard_name": "longitude"}),
            },
            dims=["lat", "lon"],
        )

        filled = fill(da, eps=1e-3)

        assert isinstance(filled, xr.DataArray)
        # Check that previously masked values are now filled
        assert not np.any(ma.getmask(filled.values))

    def test_convergence_warning(self):
        """Test warning when algorithm doesn't converge."""
        # Create data that's difficult to converge
        lons = np.linspace(0, 360, 10, endpoint=False)
        lats = np.linspace(-80, 80, 8)

        data = np.random.rand(8, 10) * 1000  # Large values
        data[2:6, 3:7] = np.nan  # Large gap

        da = xr.DataArray(
            data,
            coords={
                "lat": ("lat", lats, {"standard_name": "latitude"}),
                "lon": ("lon", lons, {"standard_name": "longitude"}),
            },
            dims=["lat", "lon"],
        )

        # Use strict tolerance and low iteration limit
        with pytest.warns(UserWarning, match="gridfill did not converge"):
            filled = fill(da, eps=1e-10, itermax=5)

    def test_history_preservation(self, sample_data):
        """Test that processing history is added to attributes."""
        # Add existing history
        sample_data.attrs["history"] = "Original processing step"

        filled = fill(sample_data, eps=1e-3)

        assert "history" in filled.attrs
        assert "gridfill" in filled.attrs["history"]
        assert "Original processing step" in filled.attrs["history"]


@pytest.mark.skipif(not HAS_XARRAY, reason="xarray not available")
class TestSpecialCases:
    """Test special cases and edge conditions."""

    def test_2d_data(self):
        """Test with simple 2D data."""
        lons = np.linspace(0, 360, 12, endpoint=False)
        lats = np.linspace(-60, 60, 8)

        data = np.random.rand(8, 12)
        data[3:5, 4:7] = np.nan

        da = xr.DataArray(
            data,
            coords={
                "lat": ("lat", lats, {"standard_name": "latitude"}),
                "lon": ("lon", lons, {"standard_name": "longitude"}),
            },
            dims=["lat", "lon"],
        )

        filled = fill(da, eps=1e-3)

        assert filled.shape == da.shape
        assert not np.any(np.isnan(filled.values))

    def test_irregular_coordinates_warning(self):
        """Test warning for irregular coordinate spacing."""
        # Create irregular longitude spacing
        lons = np.array([0, 10, 30, 60, 100, 150, 220, 300])
        lats = np.linspace(-60, 60, 8)

        data = np.random.rand(8, 8)
        data[3:5, 2:5] = np.nan

        da = xr.DataArray(
            data,
            coords={
                "lat": ("lat", lats, {"standard_name": "latitude"}),
                "lon": ("lon", lons, {"standard_name": "longitude"}),
            },
            dims=["lat", "lon"],
        )

        validation = validate_grid_coverage(da)

        # Should detect irregular spacing
        irregular_messages = [
            msg for msg in validation["messages"] if "regular" in msg.lower()
        ]
        assert len(irregular_messages) > 0

    def test_single_missing_point(self):
        """Test with only a single missing point."""
        lons = np.linspace(0, 360, 10, endpoint=False)
        lats = np.linspace(-60, 60, 8)

        data = np.random.rand(8, 10)
        data[4, 5] = np.nan  # Single missing point

        da = xr.DataArray(
            data,
            coords={
                "lat": ("lat", lats, {"standard_name": "latitude"}),
                "lon": ("lon", lons, {"standard_name": "longitude"}),
            },
            dims=["lat", "lon"],
        )

        filled = fill(da, eps=1e-3)

        assert not np.any(np.isnan(filled.values))
        # Value should be reasonable interpolation
        assert 0 < filled.values[4, 5] < 1
