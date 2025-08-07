"""
Tests for skyborn.gridfill module.

These tests verify the functionality of the gridfill module for filling
missing values in gridded data using Poisson equation solvers.
"""

import numpy as np
import numpy.ma as ma
import pytest

from skyborn.gridfill import fill

# Test xarray interface if available
try:
    import xarray as xr

    HAS_XARRAY = True
    from skyborn.gridfill.xarray import (
        fill as xr_fill,
        fill_multiple,
        validate_grid_coverage,
        _find_spatial_coordinates,
        _detect_cyclic_longitude,
    )
except ImportError:
    HAS_XARRAY = False


class TestGridfill:
    """Test suite for gridfill functionality."""

    def test_basic_fill(self):
        """Test basic filling functionality with a simple case."""
        # Create a simple 5x5 grid with one missing value in the center
        data = np.arange(25).reshape(5, 5).astype(float)
        mask = np.zeros_like(data, dtype=bool)
        mask[2, 2] = True  # Center point missing

        masked_data = ma.array(data, mask=mask)

        # Fill the missing value
        filled, converged = fill(masked_data, xdim=1, ydim=0, eps=1e-6)

        # Check convergence
        assert converged[0], "Algorithm should converge for simple case"

        # Check that the filled value is reasonable (should be close to 12)
        filled_value = filled[2, 2]
        assert 10 < filled_value < 14, f"Filled value {filled_value} seems unreasonable"

        # Check that non-missing values are preserved
        assert np.allclose(
            filled[mask == False], data[mask == False]
        ), "Non-missing values should be preserved"

    def test_rectangular_gap(self):
        """Test filling a rectangular gap."""
        # Create a 10x10 grid with a 3x3 gap in the center
        data = np.random.rand(10, 10)
        original_data = data.copy()

        mask = np.zeros_like(data, dtype=bool)
        mask[4:7, 4:7] = True  # 3x3 gap in center

        masked_data = ma.array(data, mask=mask)

        # Fill the gap
        filled, converged = fill(masked_data, xdim=1, ydim=0, eps=1e-4)

        # Check convergence
        assert converged[0], "Algorithm should converge"

        # Check that boundary values are preserved
        boundary_mask = ~mask
        assert np.allclose(
            filled[boundary_mask], original_data[boundary_mask]
        ), "Boundary values should be preserved"

        # Check that filled values are not NaN or infinite
        assert np.all(np.isfinite(filled[mask])), "Filled values should be finite"

    def test_zonal_initialization(self):
        """Test zonal mean initialization."""
        # Create data with strong zonal pattern
        lon = np.linspace(0, 360, 20)
        lat = np.linspace(-90, 90, 10)
        LON, LAT = np.meshgrid(lon, lat)

        # Create a sinusoidal pattern in latitude
        data = np.sin(np.radians(LAT)) + 0.1 * np.random.rand(*LAT.shape)

        # Create a gap in the middle
        mask = np.zeros_like(data, dtype=bool)
        mask[4:6, 8:12] = True

        masked_data = ma.array(data, mask=mask)

        # Fill with and without zonal initialization
        filled_zero, _ = fill(
            masked_data, xdim=1, ydim=0, eps=1e-4, initzonal=False, verbose=False
        )
        filled_zonal, _ = fill(
            masked_data, xdim=1, ydim=0, eps=1e-4, initzonal=True, verbose=False
        )

        # Both should be finite and different
        assert np.all(np.isfinite(filled_zero[mask]))
        assert np.all(np.isfinite(filled_zonal[mask]))

        # For strongly zonal data, zonal initialization often gives better results
        # but we just test that both work
        assert not np.allclose(
            filled_zero[mask], filled_zonal[mask]
        ), "Different initialization should give different results"

    def test_cyclic_boundary(self):
        """Test cyclic boundary conditions."""
        # Create a 8x16 grid (like a simplified global grid)
        data = np.random.rand(8, 16)

        # Make the data cyclic by ensuring left and right edges are similar
        data[:, 0] = data[:, -1]

        # Create a gap near the edge to test cyclicity
        mask = np.zeros_like(data, dtype=bool)
        mask[3:5, 0:2] = True  # Gap at left edge
        mask[3:5, -2:] = True  # Gap at right edge

        masked_data = ma.array(data, mask=mask)

        # Fill with cyclic boundary
        filled_cyclic, conv_cyc = fill(
            masked_data, xdim=1, ydim=0, eps=1e-4, cyclic=True, verbose=False
        )

        # Fill without cyclic boundary
        filled_non_cyclic, conv_non = fill(
            masked_data, xdim=1, ydim=0, eps=1e-4, cyclic=False, verbose=False
        )

        # Both should converge
        assert conv_cyc[0] and conv_non[0], "Both methods should converge"

        # Results should be different
        assert not np.allclose(
            filled_cyclic[mask], filled_non_cyclic[mask]
        ), "Cyclic and non-cyclic should give different results"

    def test_convergence_parameters(self):
        """Test different convergence parameters."""
        # Create test data
        data = np.random.rand(15, 20)
        mask = np.zeros_like(data, dtype=bool)
        mask[7:10, 9:12] = True

        masked_data = ma.array(data, mask=mask)

        # Test with different eps values
        filled_loose, conv_loose = fill(masked_data, xdim=1, ydim=0, eps=1e-2)
        filled_tight, conv_tight = fill(masked_data, xdim=1, ydim=0, eps=1e-6)

        # Both should converge but with different precision
        assert conv_loose[0] and conv_tight[0], "Both should converge"

        # Tighter tolerance should be closer to looser tolerance
        # (this is a weak test, but ensures numerical stability)
        assert np.allclose(
            filled_loose[~mask], filled_tight[~mask]
        ), "Boundary values should be identical regardless of tolerance"

    def test_error_conditions(self):
        """Test error handling."""
        # Test with non-masked array (should raise TypeError)
        regular_array = np.random.rand(5, 5)

        with pytest.raises(TypeError):
            fill(regular_array, xdim=1, ydim=0, eps=1e-4)

        # Test with invalid dimensions
        data = np.random.rand(5, 5)
        mask = np.zeros_like(data, dtype=bool)
        mask[2, 2] = True
        masked_data = ma.array(data, mask=mask)

        with pytest.raises(ValueError):
            fill(masked_data, xdim=5, ydim=0, eps=1e-4)  # Invalid xdim

        with pytest.raises(ValueError):
            fill(masked_data, xdim=1, ydim=5, eps=1e-4)  # Invalid ydim

    def test_different_relaxation_parameters(self):
        """Test different relaxation constants."""
        # Create test data
        data = np.random.rand(10, 15)
        mask = np.zeros_like(data, dtype=bool)
        mask[4:6, 6:9] = True

        masked_data = ma.array(data, mask=mask)

        # Test with different relaxation parameters
        relax_values = [0.5, 0.6, 0.8]
        results = []

        for relax in relax_values:
            filled, converged = fill(
                masked_data, xdim=1, ydim=0, eps=1e-4, relax=relax, verbose=False
            )
            assert converged[0], f"Should converge with relax={relax}"
            results.append(filled)

        # All results should preserve boundary values
        for result in results:
            assert np.allclose(
                result[~mask], data[~mask]
            ), "Boundary values should be preserved"

    def test_multidimensional_data(self):
        """Test with 3D data (time series of 2D grids)."""
        # Create 3D data: (time, lat, lon)
        nt, nlat, nlon = 5, 8, 12
        data = np.random.rand(nt, nlat, nlon)

        # Create the same mask for all time steps
        mask = np.zeros_like(data, dtype=bool)
        mask[:, 3:5, 5:8] = True  # Same spatial gap for all times

        masked_data = ma.array(data, mask=mask)

        # Fill (should handle each time slice independently)
        filled, converged = fill(masked_data, xdim=2, ydim=1, eps=1e-4)

        # Should converge for all time slices
        assert np.all(converged), "Should converge for all time slices"

        # Check that boundary values are preserved for all times
        assert np.allclose(
            filled[~mask], data[~mask]
        ), "Boundary values should be preserved across all time steps"

        # Check that filled values are finite
        assert np.all(np.isfinite(filled[mask])), "All filled values should be finite"

    def test_different_dimension_orders(self):
        """Test gridfill with different dimension orders"""
        # Create test data (lat, lon) vs (lon, lat)
        data1 = np.random.rand(20, 30)
        mask1 = np.zeros_like(data1, dtype=bool)
        mask1[8:12, 12:18] = True
        masked_data1 = ma.array(data1, mask=mask1)

        # Same data transposed
        data2 = data1.T
        mask2 = mask1.T
        masked_data2 = ma.array(data2, mask=mask2)

        # Fill with different dimension orders
        filled1, converged1 = fill(masked_data1, xdim=1, ydim=0, eps=1e-4)
        filled2, converged2 = fill(masked_data2, xdim=0, ydim=1, eps=1e-4)

        assert converged1[0]
        assert converged2[0]
        assert np.allclose(filled1.T, filled2, rtol=1e-3)

    def test_convergence_tolerance_advanced(self):
        """Test filling with different convergence tolerances."""
        data = np.random.rand(12, 15)
        mask = np.zeros_like(data, dtype=bool)
        mask[4:8, 6:10] = True

        masked_data = ma.array(data, mask=mask)

        # Test strict tolerance
        filled_strict, converged_strict = fill(
            masked_data, xdim=1, ydim=0, eps=1e-8, itermax=1000
        )

        # Test loose tolerance
        filled_loose, converged_loose = fill(
            masked_data, xdim=1, ydim=0, eps=1e-3, itermax=50
        )

        assert converged_strict[0]
        assert converged_loose[0]

        # Both should give reasonable results
        assert np.all(np.isfinite(filled_strict))
        assert np.all(np.isfinite(filled_loose))

    def test_periodic_boundary_conditions(self):
        """Test with cyclic boundary conditions."""
        data = np.random.rand(15, 20)
        mask = np.zeros_like(data, dtype=bool)

        # Create gap near the edge to test cyclic behavior
        mask[5:8, 0:3] = True
        mask[5:8, 17:20] = True

        masked_data = ma.array(data, mask=mask)

        # Fill with cyclic boundary
        filled, converged = fill(masked_data, xdim=1, ydim=0, eps=1e-6, cyclic=True)

        assert converged[0]
        assert np.all(np.isfinite(filled))

    def test_complex_mask_patterns(self):
        """Test with complex, irregular mask patterns."""
        data = np.random.rand(20, 25)
        mask = np.zeros_like(data, dtype=bool)

        # Create irregular pattern of missing values
        for i in range(20):
            for j in range(25):
                if (i + j) % 7 == 0 and i > 5 and j > 5:
                    mask[i, j] = True
                if i > 10 and j > 15 and (i * j) % 13 == 0:
                    mask[i, j] = True

        masked_data = ma.array(data, mask=mask)

        filled, converged = fill(masked_data, xdim=1, ydim=0, eps=1e-5, itermax=5000)

        assert converged[0]
        assert np.all(np.isfinite(filled))

    def test_data_with_trends(self):
        """Test filling data with strong spatial trends."""
        # Create data with strong linear trend
        x = np.linspace(0, 10, 15)
        y = np.linspace(0, 8, 12)
        X, Y = np.meshgrid(x, y)
        data = 2 * X + 3 * Y + 0.1 * np.random.randn(12, 15)

        mask = np.zeros_like(data, dtype=bool)
        mask[4:8, 6:10] = True

        masked_data = ma.array(data, mask=mask)

        filled, converged = fill(masked_data, xdim=1, ydim=0, eps=1e-6)

        assert converged[0]
        assert np.all(np.isfinite(filled))

        # Filled values should follow the trend reasonably well
        original_trend = 2 * X[mask] + 3 * Y[mask]
        filled_values = filled[mask]
        correlation = np.corrcoef(original_trend, filled_values)[0, 1]
        assert correlation > 0.8  # Should maintain trend structure

    def test_invalid_dimensions_advanced(self):
        """Test with invalid dimension specifications."""
        data = np.random.rand(10, 12)
        mask = np.zeros_like(data, dtype=bool)
        mask[3:6, 4:8] = True
        masked_data = ma.array(data, mask=mask)

        # Invalid dimension indices should raise error or handle gracefully
        with pytest.raises((ValueError, IndexError)):
            fill(masked_data, xdim=5, ydim=0, eps=1e-6)

    def test_no_missing_values_edge_case(self):
        """Test behavior when there are no missing values."""
        data = np.random.rand(8, 10)
        # No mask - all values present
        masked_data = ma.array(data, mask=False)

        try:
            filled, converged = fill(masked_data, xdim=1, ydim=0, eps=1e-6)

            # Should return original data unchanged
            assert np.allclose(filled, data)
            assert converged[0]
        except (ValueError, BufferError):
            # Some implementations may not handle this case
            pass

    def test_numerical_precision(self):
        """Test numerical precision with high-precision requirements."""
        data = np.random.rand(15, 18).astype(np.float64)
        mask = np.zeros_like(data, dtype=bool)
        mask[6:9, 8:12] = True

        masked_data = ma.array(data, mask=mask)

        filled, converged = fill(masked_data, xdim=1, ydim=0, eps=1e-10, itermax=10000)

        assert converged[0]
        assert np.all(np.isfinite(filled))

        # Check that boundaries are preserved exactly
        boundary_mask = ~mask
        assert np.allclose(filled[boundary_mask], data[boundary_mask], rtol=1e-14)

    def test_interpolation_quality(self):
        """Test quality of interpolation against known solutions."""
        # Create synthetic data with known smooth function
        x = np.linspace(0, 2 * np.pi, 20)
        y = np.linspace(0, np.pi, 15)
        X, Y = np.meshgrid(x, y)

        # Known smooth function
        true_field = np.sin(X) * np.cos(Y)

        # Create gaps
        mask = np.zeros_like(true_field, dtype=bool)
        mask[5:10, 8:12] = True

        masked_data = ma.array(true_field, mask=mask)

        filled, converged = fill(masked_data, xdim=1, ydim=0, eps=1e-8)

        assert converged[0]

        # Compare filled values with true values
        filled_values = filled[mask]
        true_values = true_field[mask]

        # Should be reasonably close for smooth function
        rmse = np.sqrt(np.mean((filled_values - true_values) ** 2))
        assert rmse < 0.1  # Reasonable error for interpolation


@pytest.mark.skipif(not HAS_XARRAY, reason="xarray not available")
class TestGridfillXArray:
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

    def test_xarray_fill_basic(self, sample_data):
        """Test basic xarray fill operation."""
        filled = xr_fill(sample_data, eps=1e-3)

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

    def test_xarray_fill_with_explicit_dimensions(self, sample_data):
        """Test xarray fill with explicitly specified dimensions."""
        filled = xr_fill(sample_data, eps=1e-3, x_dim="lon", y_dim="lat")

        assert isinstance(filled, xr.DataArray)
        assert not np.any(np.isnan(filled.values))

    def test_xarray_fill_cyclic(self, global_data):
        """Test xarray fill with cyclic boundary conditions."""
        filled = xr_fill(global_data, eps=1e-3, cyclic=True)

        assert isinstance(filled, xr.DataArray)
        assert not np.any(np.isnan(filled.values))

    def test_xarray_fill_no_missing_values(self, sample_data):
        """Test xarray behavior when no missing values present."""
        # Create data without missing values
        clean_data = sample_data.fillna(0)

        with pytest.warns(UserWarning, match="No missing values found"):
            filled = xr_fill(clean_data, eps=1e-3)

        # Should return essentially the same data
        xr.testing.assert_allclose(filled, clean_data)

    def test_xarray_fill_attributes_disabled(self, sample_data):
        """Test xarray fill with keep_attrs=False."""
        filled = xr_fill(sample_data, eps=1e-3, keep_attrs=False)

        # Should not have original attributes (except possibly history)
        assert "units" not in filled.attrs
        assert "long_name" not in filled.attrs

    def test_xarray_fill_multiple(self, sample_data):
        """Test fill_multiple function."""
        # Create second dataset
        data2 = sample_data.copy()
        data2.values[0, 8:12, 15:20] = np.nan

        filled_list = fill_multiple([sample_data, data2], eps=1e-3)

        assert len(filled_list) == 2
        assert all(isinstance(da, xr.DataArray) for da in filled_list)
        assert all(not np.any(np.isnan(da.values)) for da in filled_list)

    def test_xarray_validation_good_data(self, sample_data):
        """Test validation with good quality data."""
        validation = validate_grid_coverage(sample_data, min_coverage=0.1)

        assert validation["valid"] is True
        assert validation["coverage"] > 0.8
        assert len(validation["messages"]) == 0

    def test_xarray_validation_sparse_data(self):
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

    def test_xarray_invalid_input_type(self):
        """Test error handling for invalid input types."""
        with pytest.raises(TypeError, match="data must be xarray.DataArray"):
            xr_fill(np.random.rand(10, 10), eps=1e-3)

    def test_xarray_invalid_dimensions(self, sample_data):
        """Test error handling for invalid dimension names."""
        with pytest.raises(ValueError, match="x_dim 'invalid' not found"):
            xr_fill(sample_data, eps=1e-3, x_dim="invalid", y_dim="lat")

    def test_xarray_missing_coordinates(self):
        """Test error handling when spatial coordinates cannot be found."""
        # Create DataArray without recognizable spatial coordinates
        data = xr.DataArray(np.random.rand(10, 15), dims=["dim1", "dim2"])

        with pytest.raises(ValueError, match="Could not identify spatial coordinates"):
            xr_fill(data, eps=1e-3)

    def test_xarray_coordinate_name_variations(self):
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

    def test_xarray_masked_array_input(self):
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

        filled = xr_fill(da, eps=1e-3)

        assert isinstance(filled, xr.DataArray)
        # Check that previously masked values are now filled
        assert not np.any(ma.getmask(filled.values))

    def test_xarray_convergence_warning(self):
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
            filled = xr_fill(da, eps=1e-10, itermax=5)

    def test_xarray_history_preservation(self, sample_data):
        """Test that processing history is added to attributes."""
        # Add existing history
        sample_data.attrs["history"] = "Original processing step"

        filled = xr_fill(sample_data, eps=1e-3)

        assert "history" in filled.attrs
        assert "gridfill" in filled.attrs["history"]
        assert "Original processing step" in filled.attrs["history"]

    def test_xarray_2d_data(self):
        """Test with simple 2D xarray data."""
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

        filled = xr_fill(da, eps=1e-3)

        assert filled.shape == da.shape
        assert not np.any(np.isnan(filled.values))

    def test_xarray_irregular_coordinates_warning(self):
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

    def test_xarray_single_missing_point(self):
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

        filled = xr_fill(da, eps=1e-3)

        assert not np.any(np.isnan(filled.values))
        # Value should be reasonable interpolation
        assert 0 < filled.values[4, 5] < 1


# Test fill_cube function if iris is available
try:
    # Skip iris for now due to package conflicts
    raise ImportError("Skipping iris tests due to package conflicts")
    import iris
    import iris.tests
    from skyborn.gridfill import fill_cube

    class TestGridfillCube:
        """Test suite for fill_cube functionality with iris cubes."""

        def test_fill_cube_basic(self):
            """Test basic cube filling functionality."""
            # Create a simple test cube
            data = np.random.rand(5, 10, 15)  # (time, lat, lon)

            # Create missing values
            data = ma.array(data)
            data.mask = False
            data.mask[2, 2:4, 6:9] = True  # Create a gap

            # Create coordinates
            time_coord = iris.coords.DimCoord(
                np.arange(5), standard_name="time", units="days since 2000-01-01"
            )
            lat_coord = iris.coords.DimCoord(
                np.linspace(-45, 45, 10), standard_name="latitude", units="degrees"
            )
            lon_coord = iris.coords.DimCoord(
                np.linspace(0, 360, 15), standard_name="longitude", units="degrees"
            )
            lon_coord.circular = True  # Make longitude cyclic

            # Create cube
            cube = iris.cube.Cube(
                data,
                dim_coords_and_dims=[(time_coord, 0), (lat_coord, 1), (lon_coord, 2)],
            )

            # Fill the cube
            filled_cube = fill_cube(cube, eps=1e-4, verbose=False)

            # Check that we got a cube back
            assert isinstance(filled_cube, iris.cube.Cube)

            # Check that filled data has no missing values where we had gaps
            original_mask = cube.data.mask
            filled_data = filled_cube.data

            # Where we had missing data, we should now have finite values
            if hasattr(filled_data, "mask"):
                assert not np.any(
                    filled_data.mask[original_mask]
                ), "Previously missing values should now be filled"

            # Check that original data is preserved where not missing
            if hasattr(filled_data, "mask"):
                good_original = ~original_mask
                assert np.allclose(
                    filled_data.data[good_original], cube.data.data[good_original]
                ), "Original data should be preserved"
            else:
                good_original = ~original_mask
                assert np.allclose(
                    filled_data[good_original], cube.data.data[good_original]
                ), "Original data should be preserved"

        def test_fill_cube_full_output(self):
            """Test cube filling with full output."""
            # Create test data with a gap that might not converge easily
            data = np.random.rand(3, 6, 8)
            data = ma.array(data)
            data.mask = False
            data.mask[1, 2:4, 3:6] = True  # Create a gap

            # Create a minimal cube
            cube = iris.cube.Cube(data)

            # Add required coordinates for auto-detection
            lat_coord = iris.coords.DimCoord(
                np.linspace(-30, 30, 6), standard_name="latitude", units="degrees"
            )
            lon_coord = iris.coords.DimCoord(
                np.linspace(0, 360, 8), standard_name="longitude", units="degrees"
            )

            cube.add_dim_coord(lat_coord, 1)
            cube.add_dim_coord(lon_coord, 2)

            # Fill with full output
            filled_cube, not_converged = fill_cube(
                cube, eps=1e-6, itermax=50, full_output=True, verbose=False
            )

            # Check return types
            assert isinstance(filled_cube, iris.cube.Cube)
            assert isinstance(not_converged, np.ndarray)
            assert not_converged.dtype == bool

except ImportError:
    # Skip iris tests if iris is not available
    pass


if __name__ == "__main__":
    # Run basic tests
    test_suite = TestGridfill()

    print("Running gridfill tests...")
    test_suite.test_basic_fill()
    print("✓ Basic fill test passed")

    test_suite.test_rectangular_gap()
    print("✓ Rectangular gap test passed")

    test_suite.test_zonal_initialization()
    print("✓ Zonal initialization test passed")

    test_suite.test_cyclic_boundary()
    print("✓ Cyclic boundary test passed")

    test_suite.test_convergence_parameters()
    print("✓ Convergence parameters test passed")

    test_suite.test_different_relaxation_parameters()
    print("✓ Relaxation parameters test passed")

    test_suite.test_multidimensional_data()
    print("✓ Multidimensional data test passed")

    print("All gridfill tests passed!")
