"""
Tests for skyborn.interp.regridding module.

This module tests the regridding and interpolation functionality
in the skyborn.interp.regridding module.
"""

import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_array_almost_equal, assert_array_equal
from sklearn import neighbors

import skyborn.interp.regridding as regridding_module
from skyborn.interp.regridding import (
    BilinearRegridder,
    ConservativeRegridder,
    Grid,
    NearestRegridder,
    Regridder,
    _assert_increasing,
    _conservative_latitude_weights,
    _conservative_longitude_weights,
    _detect_coordinate_names,
    _latitude_overlap,
    _longitude_overlap,
    nearest_neighbor_indices,
    regrid_dataset,
)


def _balltree_haversine_indices(source_grid: Grid, target_grid: Grid) -> np.ndarray:
    source_mesh = np.meshgrid(source_grid.lon, source_grid.lat, indexing="ij")
    target_mesh = np.meshgrid(target_grid.lon, target_grid.lat, indexing="ij")
    source_coords = np.stack([source_mesh[1].ravel(), source_mesh[0].ravel()], axis=-1)
    target_coords = np.stack([target_mesh[1].ravel(), target_mesh[0].ravel()], axis=-1)
    tree = neighbors.BallTree(source_coords, metric="haversine")
    return tree.query(target_coords, return_distance=False).squeeze(axis=-1)


class TestGrid:
    """Test Grid class functionality."""

    def test_grid_creation_from_degrees(self):
        """Test Grid creation from degree values."""
        lon_deg = np.array([0, 90, 180, 270])
        lat_deg = np.array([-90, -45, 0, 45, 90])

        grid = Grid.from_degrees(lon_deg, lat_deg)

        # Check that values are converted to radians
        expected_lon = np.deg2rad(lon_deg)
        expected_lat = np.deg2rad(lat_deg)

        assert_array_almost_equal(grid.lon, expected_lon)
        assert_array_almost_equal(grid.lat, expected_lat)

    def test_grid_creation_from_dataset(self):
        """Test Grid creation from xarray Dataset."""
        lon_vals = np.linspace(0, 360, 73)
        lat_vals = np.linspace(-90, 90, 37)

        dataset = xr.Dataset(
            {"temperature": (["lat", "lon"], np.random.randn(37, 73))},
            coords={"lat": lat_vals, "lon": lon_vals},
        )

        grid = Grid.from_dataset(dataset)

        expected_lon = np.deg2rad(lon_vals)
        expected_lat = np.deg2rad(lat_vals)

        assert_array_almost_equal(grid.lon, expected_lon)
        assert_array_almost_equal(grid.lat, expected_lat)

    def test_grid_shape_property(self):
        """Test Grid shape property."""
        lon_deg = np.linspace(0, 360, 144)
        lat_deg = np.linspace(-90, 90, 73)

        grid = Grid.from_degrees(lon_deg, lat_deg)

        assert grid.shape == (144, 73)

    def test_grid_equality_and_hashing(self):
        """Test Grid equality and hashing."""
        lon_deg = np.array([0, 90, 180])
        lat_deg = np.array([-45, 0, 45])

        grid1 = Grid.from_degrees(lon_deg, lat_deg)
        grid2 = Grid.from_degrees(lon_deg, lat_deg)
        grid3 = Grid.from_degrees(lon_deg + 1, lat_deg)

        # Test equality
        assert grid1 == grid2
        assert grid1 != grid3

        # Test hashing (grids should be hashable)
        grid_set = {grid1, grid2, grid3}
        assert len(grid_set) == 2  # grid1 and grid2 are the same

    def test_grid_from_dataset_different_coord_names(self):
        """Test Grid creation with different coordinate naming conventions."""
        lon_vals = np.linspace(0, 360, 10)
        lat_vals = np.linspace(-90, 90, 5)

        # Test with 'longitude' and 'latitude'
        dataset = xr.Dataset(
            {"temp": (["latitude", "longitude"], np.random.randn(5, 10))},
            coords={"latitude": lat_vals, "longitude": lon_vals},
        )

        grid = Grid.from_dataset(dataset)
        assert len(grid.lon) == 10
        assert len(grid.lat) == 5

        # Test with 'x' and 'y'
        dataset_xy = xr.Dataset(
            {"temp": (["y", "x"], np.random.randn(5, 10))},
            coords={"y": lat_vals, "x": lon_vals},
        )

        grid_xy = Grid.from_dataset(dataset_xy)
        assert len(grid_xy.lon) == 10
        assert len(grid_xy.lat) == 5


class TestCoordinateDetection:
    """Test coordinate name detection functionality."""

    def test_detect_coordinate_names_standard(self):
        """Test coordinate detection with standard names."""
        dataset = xr.Dataset(
            {"temp": (["lat", "lon"], np.random.randn(5, 10))},
            coords={"lat": np.linspace(-90, 90, 5), "lon": np.linspace(0, 360, 10)},
        )

        lon_name, lat_name = _detect_coordinate_names(dataset)
        assert lon_name == "lon"
        assert lat_name == "lat"

    def test_detect_coordinate_names_alternative(self):
        """Test coordinate detection with alternative names."""
        dataset = xr.Dataset(
            {"temp": (["latitude", "longitude"], np.random.randn(5, 10))},
            coords={
                "latitude": np.linspace(-90, 90, 5),
                "longitude": np.linspace(0, 360, 10),
            },
        )

        lon_name, lat_name = _detect_coordinate_names(dataset)
        assert lon_name == "longitude"
        assert lat_name == "latitude"

    def test_detect_coordinate_names_failure(self):
        """Test coordinate detection failure."""
        dataset = xr.Dataset(
            {"temp": (["unknown_x", "unknown_y"], np.random.randn(5, 10))},
            coords={
                "unknown_x": np.linspace(-90, 90, 5),
                "unknown_y": np.linspace(0, 360, 10),
            },
        )

        with pytest.raises(ValueError, match="Could not detect longitude/latitude"):
            _detect_coordinate_names(dataset)


class TestNearestRegridder:
    """Test nearest neighbor regridding functionality."""

    @pytest.fixture
    def sample_grids(self):
        """Create sample source and target grids."""
        source_lon = np.linspace(0, 360, 9, endpoint=False)  # 8 points
        source_lat = np.linspace(-60, 60, 7)  # 7 points
        source_grid = Grid.from_degrees(source_lon, source_lat)

        target_lon = np.linspace(0, 360, 5, endpoint=False)  # 4 points
        target_lat = np.linspace(-40, 40, 5)  # 5 points
        target_grid = Grid.from_degrees(target_lon, target_lat)

        return source_grid, target_grid

    def test_nearest_regridder_creation(self, sample_grids):
        """Test NearestRegridder creation."""
        source_grid, target_grid = sample_grids
        regridder = NearestRegridder(source_grid, target_grid)

        assert regridder.source == source_grid
        assert regridder.target == target_grid

    def test_nearest_regridder_indices_caching(self, sample_grids):
        """Test that indices are computed and cached correctly."""
        source_grid, target_grid = sample_grids
        regridder = NearestRegridder(source_grid, target_grid)

        # First access should compute indices
        indices1 = regridder.indices
        # Second access should return cached indices
        indices2 = regridder.indices

        assert indices1 is indices2  # Should be the same object
        assert len(indices1) == target_grid.shape[0] * target_grid.shape[1]

    def test_nearest_regridder_array(self, sample_grids):
        """Test array regridding with nearest neighbor."""
        source_grid, target_grid = sample_grids
        regridder = NearestRegridder(source_grid, target_grid)

        # Create test data
        source_data = np.random.randn(*source_grid.shape)
        result = regridder.regrid_array(source_data)

        assert result.shape == target_grid.shape
        assert np.all(np.isfinite(result))

    def test_nearest_regridder_dataset(self, sample_grids):
        """Test dataset regridding with nearest neighbor."""
        source_grid, target_grid = sample_grids
        regridder = NearestRegridder(source_grid, target_grid)

        # Create test dataset
        dataset = xr.Dataset(
            {
                "temperature": (
                    ["lat", "lon"],
                    np.random.randn(*source_grid.shape[::-1]),
                )
            },  # Note: xarray uses (lat, lon) order
            coords={
                "lat": np.rad2deg(source_grid.lat),
                "lon": np.rad2deg(source_grid.lon),
            },
        )

        result = regridder.regrid_dataset(dataset)

        assert "temperature" in result.data_vars
        # (lat, lon)
        assert result["temperature"].shape == target_grid.shape[::-1]
        assert np.all(np.isfinite(result["temperature"].values))

    def test_nearest_neighbor_indices_function(self, sample_grids):
        """Test standalone nearest neighbor indices function."""
        source_grid, target_grid = sample_grids
        indices = nearest_neighbor_indices(source_grid, target_grid)

        assert len(indices) == target_grid.shape[0] * target_grid.shape[1]
        assert np.all(indices >= 0)
        assert np.all(indices < source_grid.shape[0] * source_grid.shape[1])

    def test_native_regrid_module_imports_under_new_name(self):
        """The native backend should be importable as skyborn.interp.regrid."""
        native_regrid_module = pytest.importorskip("skyborn.interp.regrid")
        assert hasattr(native_regrid_module, "nearest_neighbor_indices")
        assert hasattr(native_regrid_module, "nearest_regrid_apply")

    def test_nearest_neighbor_indices_identity_grid(self):
        """Test that identical regular grids map to the identity permutation."""
        grid = Grid.from_degrees(
            np.array([0.0, 90.0, 180.0, 270.0]),
            np.array([-45.0, 0.0, 45.0]),
        )

        indices = nearest_neighbor_indices(grid, grid)

        assert_array_equal(indices, np.arange(grid.shape[0] * grid.shape[1]))

    def test_nearest_regridder_identity_grid(self):
        """Nearest regridding should preserve a field on the same grid."""
        grid = Grid.from_degrees(
            np.array([0.0, 90.0, 180.0, 270.0]),
            np.array([-45.0, 0.0, 45.0]),
        )
        regridder = NearestRegridder(grid, grid)

        source_data = np.arange(grid.shape[0] * grid.shape[1]).reshape(grid.shape)
        result = regridder.regrid_array(source_data)

        assert_array_equal(result, source_data)

    @pytest.mark.parametrize("shape_prefix", [(), (2,), (2, 3)])
    def test_nearest_regridder_nd_matches_flat_gather(self, sample_grids, shape_prefix):
        """Nearest regridding should apply the cached flat indices across leading dims."""
        source_grid, target_grid = sample_grids
        regridder = NearestRegridder(source_grid, target_grid)

        source_data = np.arange(
            np.prod(shape_prefix + source_grid.shape),
            dtype=float,
        ).reshape(shape_prefix + source_grid.shape)
        result = regridder.regrid_array(source_data)

        flat = source_data.reshape(-1, source_grid.shape[0] * source_grid.shape[1])
        expected = np.take(flat, regridder.indices, axis=1).reshape(
            shape_prefix + target_grid.shape
        )

        assert_array_equal(result, expected)

    def test_nearest_neighbor_indices_match_balltree_reference(self):
        """Nearest indices should agree with the BallTree haversine reference."""
        source_grid = Grid.from_degrees(
            np.linspace(0, 360, 16, endpoint=False),
            np.linspace(-80, 80, 9),
        )
        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 23, endpoint=False),
            np.linspace(-75, 75, 13),
        )

        indices = nearest_neighbor_indices(source_grid, target_grid)
        expected = _balltree_haversine_indices(source_grid, target_grid)

        assert_array_equal(indices, expected)

    @pytest.mark.parametrize("shape_prefix", [(), (2,), (2, 3)])
    def test_nearest_native_matches_python_fallback(self, monkeypatch, shape_prefix):
        """Nearest backend choice should not change 2D/3D/4D outputs."""
        source_grid = Grid.from_degrees(
            np.linspace(0, 360, 16, endpoint=False),
            np.linspace(-80, 80, 9),
        )
        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 23, endpoint=False),
            np.linspace(-75, 75, 13),
        )
        native_regridder = NearestRegridder(source_grid, target_grid)
        source_data = np.arange(
            np.prod(shape_prefix + source_grid.shape),
            dtype=float,
        ).reshape(shape_prefix + source_grid.shape)
        native_result = native_regridder.regrid_array(source_data)

        monkeypatch.setattr(
            regridding_module,
            "_native_nearest_neighbor_indices",
            None,
        )
        monkeypatch.setattr(
            regridding_module,
            "_native_nearest_regrid_apply",
            None,
        )
        fallback_regridder = NearestRegridder(source_grid, target_grid)
        fallback_result = fallback_regridder.regrid_array(source_data)

        assert_array_equal(native_result, fallback_result)

    @pytest.mark.parametrize("shape_prefix", [(), (2,), (2, 3)])
    def test_nearest_native_apply_matches_numpy_take(self, shape_prefix):
        """Native nearest gather should match the NumPy flat-take reference exactly."""
        source_grid = Grid.from_degrees(
            np.linspace(0, 360, 24, endpoint=False),
            np.linspace(-85, 85, 15),
        )
        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 37, endpoint=False),
            np.linspace(-80, 80, 17),
        )
        regridder = NearestRegridder(source_grid, target_grid)
        source_data = np.arange(
            np.prod(shape_prefix + source_grid.shape),
            dtype=np.float64,
        ).reshape(shape_prefix + source_grid.shape)

        result = regridder.regrid_array(source_data)
        flat = source_data.reshape(-1, source_grid.shape[0] * source_grid.shape[1])
        expected = np.take(flat, regridder.indices, axis=1).reshape(
            shape_prefix + target_grid.shape
        )

        assert_array_equal(result, expected)

    def test_nearest_polar_tie_outputs_stay_on_equivalent_poles(self):
        """Exact-pole source rows may choose different longitudes but stay on the same pole."""
        source_grid = Grid.from_degrees(
            np.linspace(0, 360, 16, endpoint=False),
            np.linspace(-90, 90, 9),
        )
        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 23, endpoint=False),
            np.linspace(-85, 85, 13),
        )
        source_data = np.tile(
            np.arange(source_grid.shape[1], dtype=float), (source_grid.shape[0], 1)
        )

        native_regridder = NearestRegridder(source_grid, target_grid)
        native_result = native_regridder.regrid_array(source_data)

        source_mesh = np.meshgrid(source_grid.lon, source_grid.lat, indexing="ij")
        target_mesh = np.meshgrid(target_grid.lon, target_grid.lat, indexing="ij")
        source_coords = np.stack(
            [source_mesh[1].ravel(), source_mesh[0].ravel()], axis=-1
        )
        target_coords = np.stack(
            [target_mesh[1].ravel(), target_mesh[0].ravel()], axis=-1
        )
        expected_indices = (
            neighbors.BallTree(source_coords, metric="haversine")
            .query(target_coords, return_distance=False)
            .squeeze(axis=-1)
        )
        fallback_result = source_data.ravel()[expected_indices].reshape(
            target_grid.shape
        )

        assert_array_equal(native_result, fallback_result)

    def test_nearest_neighbor_2d_wrong_shape_raises(self, sample_grids):
        """Test the private 2D kernel rejects arrays with the wrong source shape."""
        source_grid, target_grid = sample_grids
        regridder = NearestRegridder(source_grid, target_grid)

        with pytest.raises(ValueError, match="to match source.shape"):
            regridder._nearest_neighbor_2d(np.zeros((2, 2)))

    def test_nearest_neighbor_2d_matches_cached_indices(self, sample_grids):
        """The private 2D nearest kernel should apply the cached flat indices."""
        source_grid, target_grid = sample_grids
        regridder = NearestRegridder(source_grid, target_grid)

        source_data = np.arange(np.prod(source_grid.shape), dtype=float).reshape(
            source_grid.shape
        )
        result = regridder._nearest_neighbor_2d(source_data)
        expected = source_data.ravel()[regridder.indices].reshape(target_grid.shape)

        assert_array_equal(result, expected)

    def test_nearest_neighbor_indices_uses_native_helper_when_enabled(
        self, sample_grids, monkeypatch
    ):
        """The native nearest-index helper should be used only when the feature flag is enabled."""
        source_grid, target_grid = sample_grids
        expected = np.arange(target_grid.shape[0] * target_grid.shape[1], dtype=np.intp)
        calls = []

        def fake_native(src_lon, src_lat, tgt_lon, tgt_lat):
            calls.append((src_lon.shape, src_lat.shape, tgt_lon.shape, tgt_lat.shape))
            return expected

        monkeypatch.setattr(regridding_module, "_ENABLE_NATIVE_NEAREST", True)
        monkeypatch.setattr(
            regridding_module,
            "_native_nearest_neighbor_indices",
            fake_native,
        )

        result = nearest_neighbor_indices(source_grid, target_grid)

        assert_array_equal(result, expected)
        assert calls == [
            (
                source_grid.lon.shape,
                source_grid.lat.shape,
                target_grid.lon.shape,
                target_grid.lat.shape,
            )
        ]

    @pytest.mark.parametrize("shape_prefix", [(), (2,), (2, 3)])
    def test_nearest_regridder_nd_falls_back_to_flat_gather(
        self, sample_grids, monkeypatch, shape_prefix
    ):
        """Nearest ND regridding should fall back to flat gather when native apply is absent."""
        source_grid, target_grid = sample_grids
        regridder = NearestRegridder(source_grid, target_grid)

        monkeypatch.setattr(
            regridding_module,
            "_native_nearest_regrid_apply",
            None,
        )

        source_data = np.arange(
            np.prod(shape_prefix + source_grid.shape),
            dtype=float,
        ).reshape(shape_prefix + source_grid.shape)
        result = regridder.regrid_array(source_data)
        expected = regridding_module._gather_flat_spatial(
            source_data,
            regridder.indices,
            source_grid.shape,
            target_grid.shape,
        )

        assert_array_equal(result, expected)

    def test_nearest_regridder_wrong_shape_raises(self, sample_grids):
        """Nearest ND regridding should reject arrays with the wrong spatial shape."""
        source_grid, target_grid = sample_grids
        regridder = NearestRegridder(source_grid, target_grid)

        wrong = np.zeros((source_grid.shape[0], source_grid.shape[1] + 1))

        with pytest.raises(ValueError, match="Expected field shape"):
            regridder.regrid_array(wrong)


def test_gather_flat_spatial_applies_indices_across_leading_dimensions():
    """The flat gather helper should preserve leading dimensions while remapping space."""
    field = np.arange(2 * 3 * 4 * 5, dtype=float).reshape(2, 3, 4, 5)
    indices = np.array([0, 6, 11, 19], dtype=np.intp)

    result = regridding_module._gather_flat_spatial(field, indices, (4, 5), (2, 2))
    expected = np.take(field.reshape(-1, 20), indices, axis=1).reshape(2, 3, 2, 2)

    assert_array_equal(result, expected)


class TestBilinearRegridder:
    """Test bilinear regridding functionality."""

    @pytest.fixture
    def sample_grids(self):
        """Create sample grids for bilinear interpolation."""
        source_lon = np.linspace(0, 360, 17, endpoint=False)
        source_lat = np.linspace(-80, 80, 9)
        source_grid = Grid.from_degrees(source_lon, source_lat)

        target_lon = np.linspace(0, 360, 9, endpoint=False)
        target_lat = np.linspace(-60, 60, 7)
        target_grid = Grid.from_degrees(target_lon, target_lat)

        return source_grid, target_grid

    def test_bilinear_regridder_creation(self, sample_grids):
        """Test BilinearRegridder creation."""
        source_grid, target_grid = sample_grids
        regridder = BilinearRegridder(source_grid, target_grid)

        assert regridder.source == source_grid
        assert regridder.target == target_grid

    def test_bilinear_regridder_array(self, sample_grids):
        """Test array regridding with bilinear interpolation."""
        source_grid, target_grid = sample_grids
        regridder = BilinearRegridder(source_grid, target_grid)

        # Create smooth test data (sine function)
        lon_2d, lat_2d = np.meshgrid(source_grid.lon, source_grid.lat, indexing="ij")
        source_data = np.sin(2 * lon_2d) * np.cos(3 * lat_2d)

        result = regridder.regrid_array(source_data)

        assert result.shape == target_grid.shape
        assert np.all(np.isfinite(result))

    def test_bilinear_regridder_uses_native_regular_grid_helper(
        self, sample_grids, monkeypatch
    ):
        """Test that regular monotone grids dispatch to the native helper."""
        source_grid, target_grid = sample_grids
        regridder = BilinearRegridder(source_grid, target_grid)

        calls = []

        def fake_native(field, src_lon, src_lat, tgt_lon, tgt_lat):
            calls.append(
                (
                    field.shape,
                    src_lon.shape,
                    src_lat.shape,
                    tgt_lon.shape,
                    tgt_lat.shape,
                )
            )
            return np.zeros((len(tgt_lon), len(tgt_lat)))

        monkeypatch.setattr(
            regridding_module,
            "_native_bilinear_regrid",
            fake_native,
        )

        source_data = np.random.randn(*source_grid.shape)
        result = regridder.regrid_array(source_data)

        assert calls == [
            (
                source_data.shape,
                source_grid.lon.shape,
                source_grid.lat.shape,
                target_grid.lon.shape,
                target_grid.lat.shape,
            )
        ]
        assert result.shape == target_grid.shape

    @pytest.mark.parametrize("shape_prefix", [(2,), (2, 3)])
    def test_bilinear_regridder_nd_uses_native_regular_grid_helper(
        self, sample_grids, monkeypatch, shape_prefix
    ):
        """Multidimensional bilinear regridding should use the ND native helper."""
        source_grid, target_grid = sample_grids
        regridder = BilinearRegridder(source_grid, target_grid)

        calls = []

        def fake_native(field, src_lon, src_lat, tgt_lon, tgt_lat):
            calls.append(
                (
                    field.shape,
                    src_lon.shape,
                    src_lat.shape,
                    tgt_lon.shape,
                    tgt_lat.shape,
                )
            )
            return np.zeros(field.shape[:-2] + (len(tgt_lon), len(tgt_lat)))

        monkeypatch.setattr(
            regridding_module,
            "_native_bilinear_regrid_nd",
            fake_native,
        )

        source_data = np.random.randn(*shape_prefix, *source_grid.shape)
        result = regridder.regrid_array(source_data)

        assert calls == [
            (
                source_data.shape,
                source_grid.lon.shape,
                source_grid.lat.shape,
                target_grid.lon.shape,
                target_grid.lat.shape,
            )
        ]
        assert result.shape == shape_prefix + target_grid.shape

    def test_bilinear_regridder_linear_function(self):
        """Test bilinear interpolation with linear function (should be exact)."""
        # Create regular grids
        source_lon = np.linspace(0, 2, 5)
        source_lat = np.linspace(0, 2, 5)
        source_grid = Grid.from_degrees(source_lon, source_lat)

        target_lon = np.linspace(0, 2, 3)
        target_lat = np.linspace(0, 2, 3)
        target_grid = Grid.from_degrees(target_lon, target_lat)

        regridder = BilinearRegridder(source_grid, target_grid)

        # Create linear function
        lon_2d, lat_2d = np.meshgrid(source_grid.lon, source_grid.lat, indexing="ij")
        source_data = 2 * lon_2d + 3 * lat_2d  # Linear function

        result = regridder.regrid_array(source_data)

        # For bilinear interpolation of linear functions, result should be very close
        lon_target_2d, lat_target_2d = np.meshgrid(
            target_grid.lon, target_grid.lat, indexing="ij"
        )
        expected = 2 * lon_target_2d + 3 * lat_target_2d

        assert_array_almost_equal(result, expected, decimal=10)

    def test_bilinear_regridder_falls_back_when_native_helper_is_missing(
        self, sample_grids, monkeypatch
    ):
        """Test the Python fallback stays available when the native helper is absent."""
        source_grid, target_grid = sample_grids
        regridder = BilinearRegridder(source_grid, target_grid)

        called = False

        def fake_native(*args, **kwargs):
            nonlocal called
            called = True
            return np.zeros((len(target_grid.lon), len(target_grid.lat)))

        monkeypatch.setattr(
            regridding_module,
            "_native_bilinear_regrid",
            None,
        )

        source_data = np.random.randn(*source_grid.shape)
        result = regridder.regrid_array(source_data)

        assert not called
        assert result.shape == target_grid.shape

    @pytest.mark.parametrize("shape_prefix", [(2,), (2, 3)])
    def test_bilinear_regridder_nd_falls_back_to_python_loop(
        self, sample_grids, monkeypatch, shape_prefix
    ):
        """Multidimensional bilinear regridding should keep the Python loop fallback."""
        source_grid, target_grid = sample_grids
        regridder = BilinearRegridder(source_grid, target_grid)

        monkeypatch.setattr(
            regridding_module,
            "_native_bilinear_regrid_nd",
            None,
        )
        monkeypatch.setattr(
            regridding_module,
            "_native_bilinear_regrid",
            None,
        )

        source_data = np.random.randn(*shape_prefix, *source_grid.shape)
        result = regridder.regrid_array(source_data)
        expected = np.empty(shape_prefix + target_grid.shape, dtype=np.float64)
        for index in np.ndindex(shape_prefix):
            expected[index] = regridder._bilinear_2d(source_data[index])

        assert_array_almost_equal(result, expected, decimal=12)

    def test_bilinear_regridder_wrong_shape(self, sample_grids):
        """Test bilinear regridder with wrong input shape."""
        source_grid, target_grid = sample_grids
        regridder = BilinearRegridder(source_grid, target_grid)

        # Create data with wrong shape
        wrong_data = np.random.randn(5, 5)  # Wrong shape

        with pytest.raises(ValueError, match="Expected field shape"):
            regridder.regrid_array(wrong_data)

    def test_bilinear_2d_wrong_shape_raises(self, sample_grids):
        """The private 2D bilinear kernel should reject arrays with the wrong source shape."""
        source_grid, target_grid = sample_grids
        regridder = BilinearRegridder(source_grid, target_grid)

        with pytest.raises(ValueError, match="Expected field shape"):
            regridder._bilinear_2d(np.zeros((2, 2)))

    @pytest.mark.parametrize("shape_prefix", [(2,), (2, 3)])
    def test_bilinear_regridder_multidimensional_matches_stacked_2d(
        self, sample_grids, shape_prefix
    ):
        """Bilinear regridding should act independently on leading dimensions."""
        source_grid, target_grid = sample_grids
        regridder = BilinearRegridder(source_grid, target_grid)

        lon_2d, lat_2d = np.meshgrid(source_grid.lon, source_grid.lat, indexing="ij")
        base = np.sin(2 * lon_2d) * np.cos(3 * lat_2d)
        source_data = np.stack(
            [base + float(i) for i in range(np.prod(shape_prefix))],
            axis=0,
        ).reshape(shape_prefix + source_grid.shape)

        result = regridder.regrid_array(source_data)
        expected = np.stack(
            [
                regridder._bilinear_2d(source_data[index])
                for index in np.ndindex(shape_prefix)
            ],
            axis=0,
        ).reshape(shape_prefix + target_grid.shape)

        assert_array_almost_equal(result, expected, decimal=12)


class TestConservativeRegridder:
    """Test conservative regridding functionality."""

    @pytest.fixture
    def sample_grids(self):
        """Create sample grids for conservative interpolation."""
        source_lon = np.linspace(0, 360, 13, endpoint=False)
        source_lat = np.linspace(-60, 60, 7)
        source_grid = Grid.from_degrees(source_lon, source_lat)

        target_lon = np.linspace(0, 360, 7, endpoint=False)
        target_lat = np.linspace(-40, 40, 5)
        target_grid = Grid.from_degrees(target_lon, target_lat)

        return source_grid, target_grid

    def test_conservative_regridder_creation(self, sample_grids):
        """Test ConservativeRegridder creation."""
        source_grid, target_grid = sample_grids
        regridder = ConservativeRegridder(source_grid, target_grid)

        assert regridder.source == source_grid
        assert regridder.target == target_grid

    def test_conservative_regridder_weights_caching(self, sample_grids):
        """Test that weights are computed and cached correctly."""
        source_grid, target_grid = sample_grids
        regridder = ConservativeRegridder(source_grid, target_grid)

        # First access should compute weights
        lon_weights1 = regridder.lon_weights
        lat_weights1 = regridder.lat_weights

        # Second access should return cached weights
        lon_weights2 = regridder.lon_weights
        lat_weights2 = regridder.lat_weights

        assert lon_weights1 is lon_weights2
        assert lat_weights1 is lat_weights2

    def test_conservative_regridder_array(self, sample_grids):
        """Test array regridding with conservative interpolation."""
        source_grid, target_grid = sample_grids
        regridder = ConservativeRegridder(source_grid, target_grid)

        # Create test data
        source_data = np.random.randn(*source_grid.shape)
        result = regridder.regrid_array(source_data)

        assert result.shape == target_grid.shape
        assert np.all(np.isfinite(result))

    def test_conservative_regridder_conservation(self, sample_grids):
        """Test that conservative regridding conserves total mass."""
        source_grid, target_grid = sample_grids
        regridder = ConservativeRegridder(source_grid, target_grid)

        # Create uniform field
        source_data = np.ones(source_grid.shape)
        result = regridder.regrid_array(source_data)

        # Conservative interpolation of constant field should remain constant
        assert_array_almost_equal(result, np.ones(target_grid.shape), decimal=10)

    def test_conservative_regridder_with_nans(self, sample_grids):
        """Test conservative regridding with NaN values."""
        source_grid, target_grid = sample_grids
        regridder = ConservativeRegridder(source_grid, target_grid)

        # Create data with some NaN values
        source_data = np.random.randn(*source_grid.shape)
        source_data[0, 0] = np.nan
        source_data[1, 1] = np.nan

        result = regridder.regrid_array(source_data)

        assert result.shape == target_grid.shape
        # Result may contain NaN but should not be all NaN
        assert not np.all(np.isnan(result))

    @pytest.mark.parametrize("shape_prefix", [(2,), (2, 3)])
    def test_conservative_regridder_nd_uses_native_helper(
        self, sample_grids, monkeypatch, shape_prefix
    ):
        """Multidimensional conservative regridding should use the native helper."""
        source_grid, target_grid = sample_grids
        regridder = ConservativeRegridder(source_grid, target_grid)

        calls = []

        def fake_native(field, lon_weights, lat_weights):
            calls.append((field.shape, lon_weights.shape, lat_weights.shape))
            return np.zeros(field.shape[:-2] + target_grid.shape)

        monkeypatch.setattr(
            regridding_module,
            "_native_conservative_regrid",
            fake_native,
        )

        source_data = np.random.randn(*shape_prefix, *source_grid.shape)
        result = regridder.regrid_array(source_data)

        assert calls == [
            (
                source_data.shape,
                (len(target_grid.lon), len(source_grid.lon)),
                (len(target_grid.lat), len(source_grid.lat)),
            )
        ]
        assert result.shape == shape_prefix + target_grid.shape

    @pytest.mark.parametrize("shape_prefix", [(), (2,), (2, 3)])
    def test_conservative_native_helper_matches_python_nanmean(
        self, sample_grids, shape_prefix
    ):
        """The optional native conservative helper should match the Python nanmean path."""
        if regridding_module._native_conservative_regrid is None:
            pytest.skip("native conservative helper is unavailable")

        source_grid, target_grid = sample_grids
        regridder = ConservativeRegridder(source_grid, target_grid)

        source_data = np.arange(
            np.prod(shape_prefix + source_grid.shape),
            dtype=float,
        ).reshape(shape_prefix + source_grid.shape)
        source_data = source_data / np.nanmax(source_data)
        source_data = source_data.copy()
        source_data[..., 0, 0] = np.nan
        source_data[..., -1, -1] = np.nan

        native_result = regridding_module._native_conservative_regrid(
            np.asarray(source_data, dtype=np.float64),
            regridder.lon_weights,
            regridder.lat_weights,
        )
        fallback_result = regridder._python_nanmean(source_data)

        assert_array_almost_equal(native_result, fallback_result, decimal=12)

    @pytest.mark.parametrize("shape_prefix", [(), (2,), (2, 3)])
    def test_conservative_regridder_falls_back_when_native_helper_is_missing(
        self, sample_grids, monkeypatch, shape_prefix
    ):
        """Conservative regridding should keep the Python fallback when native code is absent."""
        source_grid, target_grid = sample_grids
        regridder = ConservativeRegridder(source_grid, target_grid)

        monkeypatch.setattr(
            regridding_module,
            "_native_conservative_regrid",
            None,
        )

        source_data = np.random.randn(*shape_prefix, *source_grid.shape)
        source_data = source_data.copy()
        source_data[..., 0, 0] = np.nan

        result = regridder.regrid_array(source_data)
        expected = regridder._python_nanmean(source_data)

        assert result.shape == shape_prefix + target_grid.shape
        assert_array_almost_equal(result, expected, decimal=12)

    @pytest.mark.parametrize("shape_prefix", [(2,), (2, 3)])
    def test_conservative_regridder_multidimensional_constant_field(
        self, sample_grids, shape_prefix
    ):
        """Conservative regridding should preserve constant stacked fields."""
        source_grid, target_grid = sample_grids
        regridder = ConservativeRegridder(source_grid, target_grid)

        source_data = np.ones(shape_prefix + source_grid.shape, dtype=float)
        result = regridder.regrid_array(source_data)

        assert result.shape == shape_prefix + target_grid.shape
        assert_array_almost_equal(
            result,
            np.ones(shape_prefix + target_grid.shape),
            decimal=10,
        )

    def test_conservative_regridder_wrong_shape(self, sample_grids):
        """Conservative regridding should reject arrays whose spatial shape is wrong."""
        source_grid, target_grid = sample_grids
        regridder = ConservativeRegridder(source_grid, target_grid)

        wrong = np.random.randn(len(source_grid.lon), len(source_grid.lat) + 1)

        with pytest.raises(ValueError, match="Expected field shape"):
            regridder.regrid_array(wrong)

    def test_conservative_2d_wrong_shape_raises(self, sample_grids):
        """The private 2D conservative kernel should reject arrays with the wrong source shape."""
        source_grid, target_grid = sample_grids
        regridder = ConservativeRegridder(source_grid, target_grid)

        with pytest.raises(ValueError, match="Expected field shape"):
            regridder._conservative_2d(np.zeros((2, 2)))


class TestConservativeWeights:
    """Test conservative interpolation weight calculation functions."""

    def test_latitude_weights_basic(self):
        """Test basic latitude weight calculation."""
        source_lat = np.deg2rad(np.array([-60, -30, 0, 30, 60]))
        target_lat = np.deg2rad(np.array([-45, 0, 45]))

        weights = _conservative_latitude_weights(source_lat, target_lat)

        assert weights.shape == (len(target_lat), len(source_lat))
        # Each row should sum to approximately 1
        row_sums = np.sum(weights, axis=1)
        assert_array_almost_equal(row_sums, np.ones(len(target_lat)), decimal=5)

    def test_latitude_weights_python_path_normalizes_nonzero_rows(self, monkeypatch):
        """The Python latitude-weight path should normalize ordinary overlap rows."""
        source_lat = np.deg2rad(np.array([-60, -30, 0, 30, 60]))
        target_lat = np.deg2rad(np.array([-45, 0, 45]))

        monkeypatch.setattr(
            regridding_module,
            "_native_conservative_latitude_weights",
            None,
        )

        weights = _conservative_latitude_weights(source_lat, target_lat)

        assert np.all(weights >= 0)
        assert np.any(weights > 0)
        assert_array_almost_equal(
            weights.sum(axis=1), np.ones(len(target_lat)), decimal=12
        )

    def test_longitude_weights_basic(self):
        """Test basic longitude weight calculation."""
        source_lon = np.deg2rad(np.array([0, 90, 180, 270]))
        target_lon = np.deg2rad(np.array([45, 135, 225, 315]))

        weights = _conservative_longitude_weights(source_lon, target_lon)

        assert weights.shape == (len(target_lon), len(source_lon))
        # Each row should sum to approximately 1
        row_sums = np.sum(weights, axis=1)
        assert_array_almost_equal(row_sums, np.ones(len(target_lon)), decimal=5)

    def test_longitude_weights_python_path_normalizes_nonzero_rows(self, monkeypatch):
        """The Python longitude-weight path should normalize ordinary overlap rows."""
        source_lon = np.deg2rad(np.array([0, 90, 180, 270]))
        target_lon = np.deg2rad(np.array([45, 135, 225, 315]))

        monkeypatch.setattr(
            regridding_module,
            "_native_conservative_longitude_weights",
            None,
        )

        weights = _conservative_longitude_weights(source_lon, target_lon)

        assert np.all(weights >= 0)
        assert np.any(weights > 0)
        assert_array_almost_equal(
            weights.sum(axis=1), np.ones(len(target_lon)), decimal=12
        )

    def test_latitude_weights_validate_before_native_dispatch(self, monkeypatch):
        """Invalid latitude coordinates should fail before calling the native helper."""
        called = False

        def fake_native(source_points, target_points):
            nonlocal called
            called = True
            return np.empty((len(target_points), len(source_points)))

        monkeypatch.setattr(
            regridding_module,
            "_native_conservative_latitude_weights",
            fake_native,
        )

        with pytest.raises(ValueError, match="Array is not increasing"):
            _conservative_latitude_weights(
                np.deg2rad(np.array([0, -30, 30])),
                np.deg2rad(np.array([-15, 15])),
            )

        assert not called

    def test_longitude_weights_validate_before_native_dispatch(self, monkeypatch):
        """Invalid longitude coordinates should fail before calling the native helper."""
        called = False

        def fake_native(source_points, target_points):
            nonlocal called
            called = True
            return np.empty((len(target_points), len(source_points)))

        monkeypatch.setattr(
            regridding_module,
            "_native_conservative_longitude_weights",
            fake_native,
        )

        with pytest.raises(ValueError, match="Array is not increasing"):
            _conservative_longitude_weights(
                np.deg2rad(np.array([0, 180, 90])),
                np.deg2rad(np.array([45, 135])),
            )

        assert not called

    def test_latitude_weights_zero_sum_rows_fall_back_to_uniform(self, monkeypatch):
        """Test zero-overlap latitude rows use equal weights instead of dividing by zero."""
        source_lat = np.deg2rad(np.array([-60, 0, 60]))
        target_lat = np.deg2rad(np.array([-30, 30]))

        monkeypatch.setattr(
            regridding_module,
            "_native_conservative_latitude_weights",
            None,
        )
        monkeypatch.setattr(
            regridding_module,
            "_latitude_overlap",
            lambda source_points, target_points: np.zeros(
                (len(target_points), len(source_points))
            ),
        )

        weights = _conservative_latitude_weights(source_lat, target_lat)
        expected = np.full((len(target_lat), len(source_lat)), 1.0 / len(source_lat))

        assert_array_almost_equal(weights, expected)

    def test_longitude_weights_zero_sum_rows_fall_back_to_uniform(self, monkeypatch):
        """Test zero-overlap longitude rows use equal weights instead of dividing by zero."""
        source_lon = np.deg2rad(np.array([0, 90, 180, 270]))
        target_lon = np.deg2rad(np.array([45, 135]))

        monkeypatch.setattr(
            regridding_module,
            "_native_conservative_longitude_weights",
            None,
        )
        monkeypatch.setattr(
            regridding_module,
            "_longitude_overlap",
            lambda first_points, second_points: np.zeros(
                (len(first_points), len(second_points))
            ),
        )

        weights = _conservative_longitude_weights(source_lon, target_lon)
        expected = np.full((len(target_lon), len(source_lon)), 1.0 / len(source_lon))

        assert_array_almost_equal(weights, expected)

    def test_latitude_overlap_function(self):
        """Test latitude overlap calculation."""
        source_lat = np.deg2rad(np.array([-30, 0, 30]))
        target_lat = np.deg2rad(np.array([-15, 15]))

        overlap = _latitude_overlap(source_lat, target_lat)

        assert overlap.shape == (len(target_lat), len(source_lat))
        assert np.all(overlap >= 0)  # Overlaps should be non-negative

    def test_longitude_overlap_function(self):
        """Test longitude overlap calculation."""
        source_lon = np.deg2rad(np.array([0, 90, 180, 270]))
        target_lon = np.deg2rad(np.array([45, 135]))

        overlap = _longitude_overlap(target_lon, source_lon)

        assert overlap.shape == (len(target_lon), len(source_lon))
        assert np.all(overlap >= 0)  # Overlaps should be non-negative

    def test_assert_increasing_function(self):
        """Test _assert_increasing utility function."""
        # Should not raise for increasing array
        increasing_array = np.array([1, 2, 3, 4, 5])
        _assert_increasing(increasing_array)  # Should not raise

        # Should raise for non-increasing array
        non_increasing_array = np.array([1, 3, 2, 4, 5])
        with pytest.raises(ValueError, match="Array is not increasing"):
            _assert_increasing(non_increasing_array)


class TestRegridDataset:
    """Test the convenience regrid_dataset function."""

    @pytest.fixture
    def sample_dataset(self):
        """Create sample dataset for testing."""
        lon_vals = np.linspace(0, 360, 73, endpoint=False)
        lat_vals = np.linspace(-90, 90, 37)
        time_vals = np.arange(10)

        data = np.random.randn(len(time_vals), len(lat_vals), len(lon_vals))

        dataset = xr.Dataset(
            {
                "temperature": (["time", "lat", "lon"], data),
                "pressure": (
                    ["lat", "lon"],
                    np.random.randn(len(lat_vals), len(lon_vals)),
                ),
            },
            coords={"time": time_vals, "lat": lat_vals, "lon": lon_vals},
        )

        return dataset

    @pytest.fixture
    def target_grid(self):
        """Create target grid."""
        target_lon = np.linspace(0, 360, 37, endpoint=False)
        target_lat = np.linspace(-45, 45, 19)
        return Grid.from_degrees(target_lon, target_lat)

    def test_regrid_dataset_nearest(self, sample_dataset, target_grid):
        """Test regrid_dataset with nearest neighbor method."""
        result = regrid_dataset(sample_dataset, target_grid, method="nearest")

        assert "temperature" in result.data_vars
        assert "pressure" in result.data_vars
        assert result["temperature"].shape == (10, 19, 37)  # (time, lat, lon)
        assert result["pressure"].shape == (19, 37)  # (lat, lon)

    def test_regrid_dataset_bilinear(self, sample_dataset, target_grid):
        """Test regrid_dataset with bilinear method."""
        result = regrid_dataset(sample_dataset, target_grid, method="bilinear")

        assert "temperature" in result.data_vars
        assert "pressure" in result.data_vars
        assert result["temperature"].shape == (10, 19, 37)
        assert result["pressure"].shape == (19, 37)

    def test_regrid_dataset_conservative(self, sample_dataset, target_grid):
        """Test regrid_dataset with conservative method."""
        result = regrid_dataset(sample_dataset, target_grid, method="conservative")

        assert "temperature" in result.data_vars
        assert "pressure" in result.data_vars
        assert result["temperature"].shape == (10, 19, 37)
        assert result["pressure"].shape == (19, 37)

    def test_regrid_dataset_dimension_preservation(self, sample_dataset, target_grid):
        """Test that dimension order is preserved."""
        original_temp_dims = list(sample_dataset["temperature"].dims)
        original_pressure_dims = list(sample_dataset["pressure"].dims)

        result = regrid_dataset(sample_dataset, target_grid, method="nearest")

        assert list(result["temperature"].dims) == original_temp_dims
        assert list(result["pressure"].dims) == original_pressure_dims

    def test_regrid_dataset_coordinate_preservation(self, sample_dataset, target_grid):
        """Test that non-spatial coordinates are preserved."""
        result = regrid_dataset(sample_dataset, target_grid, method="nearest")

        # Time coordinate should be preserved
        assert_array_equal(result["time"].values, sample_dataset["time"].values)

        # Spatial coordinates should be updated
        assert len(result["lat"]) == 19
        assert len(result["lon"]) == 37

    def test_regrid_dataset_attrs_preservation(self, sample_dataset, target_grid):
        """Test that attributes are preserved."""
        sample_dataset.attrs["title"] = "Test Dataset"
        sample_dataset["temperature"].attrs["units"] = "K"

        result = regrid_dataset(sample_dataset, target_grid, method="nearest")

        assert result.attrs["title"] == "Test Dataset"
        assert result["temperature"].attrs["units"] == "K"

    def test_regrid_dataset_invalid_method(self, sample_dataset, target_grid):
        """Test regrid_dataset with invalid method."""
        with pytest.raises(ValueError, match="Unknown method"):
            regrid_dataset(sample_dataset, target_grid, method="invalid_method")

    def test_regrid_dataset_explicit_dimensions(self, sample_dataset, target_grid):
        """Test regrid_dataset with explicitly specified dimension names."""
        result = regrid_dataset(
            sample_dataset, target_grid, method="nearest", lon_dim="lon", lat_dim="lat"
        )

        assert "temperature" in result.data_vars
        assert result["temperature"].shape == (10, 19, 37)

    def test_regrid_dataset_subset_variables(self, sample_dataset, target_grid):
        """Test regridding subset of variables."""
        temp_only = sample_dataset[["temperature"]]
        result = regrid_dataset(temp_only, target_grid, method="bilinear")

        assert "temperature" in result.data_vars
        assert "pressure" not in result.data_vars
        assert result["temperature"].shape == (10, 19, 37)


class TestRegridderEdgeCases:
    """Test edge cases and error conditions for regridders."""

    def test_regridder_abstract_base_class(self):
        """Test that Regridder is abstract and cannot be instantiated directly."""
        source_grid = Grid.from_degrees([0, 1], [0, 1])
        target_grid = Grid.from_degrees([0, 1], [0, 1])

        regridder = Regridder(source_grid, target_grid)

        # regrid_array should raise NotImplementedError
        with pytest.raises(NotImplementedError):
            regridder.regrid_array(np.array([[1, 2], [3, 4]]))

    def test_regrid_dataset_with_non_spatial_variables(self):
        """Test regridding dataset with variables that don't have spatial dimensions."""
        lon_vals = np.linspace(0, 360, 10, endpoint=False)
        lat_vals = np.linspace(-45, 45, 5)
        time_vals = np.arange(3)

        dataset = xr.Dataset(
            {
                "temperature": (["time", "lat", "lon"], np.random.randn(3, 5, 10)),
                # No spatial dims
                "time_bounds": (["time", "bnds"], np.random.randn(3, 2)),
                # Spatial dims
                "area": (["lat", "lon"], np.random.randn(5, 10)),
            },
            coords={"time": time_vals, "lat": lat_vals, "lon": lon_vals},
        )

        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 5, endpoint=False), np.linspace(-30, 30, 3)
        )

        result = regrid_dataset(dataset, target_grid, method="nearest")

        # Variables with spatial dimensions should be regridded
        assert result["temperature"].shape == (3, 3, 5)
        assert result["area"].shape == (3, 5)

        # Variables without spatial dimensions should remain unchanged
        assert result["time_bounds"].shape == (3, 2)
        assert_array_equal(result["time_bounds"].values, dataset["time_bounds"].values)

    @pytest.mark.parametrize("method", ["nearest", "bilinear", "conservative"])
    def test_regrid_with_descending_latitude_matches_ascending(self, method):
        """Descending latitude inputs should match the equivalent ascending dataset."""
        lon_vals = np.linspace(0, 360, 8, endpoint=False)
        lat_vals = np.linspace(-60, 60, 7)
        lon_2d, lat_2d = np.meshgrid(
            np.deg2rad(lon_vals),
            np.deg2rad(lat_vals),
            indexing="xy",
        )
        values = np.sin(lon_2d) + np.cos(lat_2d)

        ascending = xr.Dataset(
            {"temp": (["lat", "lon"], values)},
            coords={"lat": lat_vals, "lon": lon_vals},
        )
        descending = ascending.isel(lat=slice(None, None, -1))

        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 4, endpoint=False), np.linspace(-30, 30, 3)
        )

        expected = regrid_dataset(ascending, target_grid, method=method)
        result = regrid_dataset(descending, target_grid, method=method)

        assert result["temp"].shape == (3, 4)
        assert_array_almost_equal(result["temp"].values, expected["temp"].values)

    def test_regrid_dataset_rejects_partial_spatial_data_variables(self):
        """Variables with only lon or only lat should fail before target-coordinate alignment."""
        lon_vals = np.linspace(0, 360, 8, endpoint=False)
        lat_vals = np.linspace(-60, 60, 7)
        dataset = xr.Dataset(
            {
                "temp": (["lat", "lon"], np.random.randn(7, 8)),
                "lat_bounds": (
                    ["lat", "bnds"],
                    np.column_stack((lat_vals - 5.0, lat_vals + 5.0)),
                ),
            },
            coords={"lat": lat_vals, "lon": lon_vals},
        )
        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 4, endpoint=False), np.linspace(-30, 30, 3)
        )

        with pytest.raises(ValueError, match="lat_bounds.*only one horizontal"):
            regrid_dataset(dataset, target_grid, method="nearest")

    def test_regrid_dataset_rejects_partial_spatial_coordinates(self):
        """Coordinates with only lon or only lat should fail before target-coordinate alignment."""
        lon_vals = np.linspace(0, 360, 8, endpoint=False)
        lat_vals = np.linspace(-60, 60, 7)
        dataset = xr.Dataset(
            {"temp": (["lat", "lon"], np.random.randn(7, 8))},
            coords={
                "lat": lat_vals,
                "lon": lon_vals,
                "lat_bounds": (
                    ["lat", "bnds"],
                    np.column_stack((lat_vals - 5.0, lat_vals + 5.0)),
                ),
            },
        )
        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 4, endpoint=False), np.linspace(-30, 30, 3)
        )

        with pytest.raises(ValueError, match="lat_bounds.*only one horizontal"):
            regrid_dataset(dataset, target_grid, method="nearest")

    def test_regrid_dataset_rejects_nonmonotonic_latitude(self):
        """Latitude coordinates must be strictly ascending or descending."""
        lon_vals = np.linspace(0, 360, 8, endpoint=False)
        lat_vals = np.array([-60.0, 0.0, -30.0, 30.0])
        dataset = xr.Dataset(
            {"temp": (["lat", "lon"], np.random.randn(4, 8))},
            coords={"lat": lat_vals, "lon": lon_vals},
        )
        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 4, endpoint=False), np.linspace(-30, 30, 3)
        )

        with pytest.raises(ValueError, match="strictly monotonic"):
            regrid_dataset(dataset, target_grid, method="nearest")

    def test_regridder_rejects_source_grid_mismatch(self):
        """A direct Regridder should fail clearly when its source grid differs from the dataset."""
        lon_vals = np.linspace(0, 360, 8, endpoint=False)
        lat_vals = np.linspace(-60, 60, 7)
        dataset = xr.Dataset(
            {"temp": (["lat", "lon"], np.random.randn(7, 8))},
            coords={"lat": lat_vals, "lon": lon_vals},
        )
        source_grid = Grid.from_degrees(lon_vals + 1.0, lat_vals)
        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 4, endpoint=False), np.linspace(-30, 30, 3)
        )
        regridder = NearestRegridder(source_grid, target_grid)

        with pytest.raises(ValueError, match="source grid does not match"):
            regridder.regrid_dataset(dataset)

    def test_regrid_single_point_grids(self):
        """Test regridding with single point grids."""
        # Single point source and target
        source_grid = Grid.from_degrees([180], [0])
        target_grid = Grid.from_degrees([0], [0])

        regridder = NearestRegridder(source_grid, target_grid)

        source_data = np.array([[5.0]])
        result = regridder.regrid_array(source_data)

        assert result.shape == (1, 1)
        assert result[0, 0] == 5.0

    def test_conservative_regridder_identical_grids(self):
        """Test conservative regridding with identical source and target grids."""
        grid = Grid.from_degrees(
            np.linspace(0, 360, 8, endpoint=False), np.linspace(-45, 45, 5)
        )

        regridder = ConservativeRegridder(grid, grid)

        source_data = np.random.randn(*grid.shape)
        result = regridder.regrid_array(source_data)

        # Should be very close to identity
        assert_array_almost_equal(result, source_data, decimal=10)


class TestPerformanceOptimizations:
    """Test performance optimization features."""

    def test_weight_caching_conservative(self):
        """Test that conservative regridder caches weights for performance."""
        source_grid = Grid.from_degrees(
            np.linspace(0, 360, 20, endpoint=False), np.linspace(-60, 60, 10)
        )
        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 10, endpoint=False), np.linspace(-30, 30, 5)
        )

        regridder = ConservativeRegridder(source_grid, target_grid)

        # Initially weights should not be computed
        assert regridder._lon_weights is None
        assert regridder._lat_weights is None

        # First access computes weights
        lon_weights = regridder.lon_weights
        lat_weights = regridder.lat_weights

        # Now they should be cached
        assert regridder._lon_weights is not None
        assert regridder._lat_weights is not None

        # Second access returns cached values
        assert regridder.lon_weights is lon_weights
        assert regridder.lat_weights is lat_weights

    def test_indices_caching_nearest(self):
        """Test that nearest regridder caches indices for performance."""
        source_grid = Grid.from_degrees(
            np.linspace(0, 360, 20, endpoint=False), np.linspace(-60, 60, 10)
        )
        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 10, endpoint=False), np.linspace(-30, 30, 5)
        )

        regridder = NearestRegridder(source_grid, target_grid)

        # Initially indices should not be computed
        assert regridder._indices is None

        # First access computes indices
        indices = regridder.indices

        # Now they should be cached
        assert regridder._indices is not None

        # Second access returns cached values
        assert regridder.indices is indices


# Integration tests
class TestRegridIntegration:
    """Integration tests with realistic climate data scenarios."""

    def test_native_regrid_module_imports_under_new_name(self):
        """The compiled accelerator should be importable as skyborn.interp.regrid."""
        native_regrid_module = pytest.importorskip("skyborn.interp.regrid")
        assert hasattr(native_regrid_module, "nearest_neighbor_indices")
        assert hasattr(native_regrid_module, "nearest_regrid_apply")

    def test_climate_data_regridding_workflow(self):
        """Test complete workflow similar to climate data processing."""
        # Simulate CMIP-like data
        time = np.arange(12)  # 12 months
        source_lon = np.linspace(0, 360, 144, endpoint=False)  # 2.5 degree
        source_lat = np.linspace(-90, 90, 73)  # 2.5 degree

        # Create realistic-looking precipitation data
        np.random.seed(42)
        pr_data = np.abs(np.random.gamma(2, 2, (12, 73, 144)))  # Always positive

        dataset = xr.Dataset(
            {"pr": (["time", "lat", "lon"], pr_data)},
            coords={
                "time": time,
                "lat": source_lat,
                "lon": source_lon,
            },
            attrs={"title": "CMIP6 precipitation data"},
        )

        # Target high-resolution grid
        target_lon = np.arange(0, 361, 1)  # 1 degree
        target_lat = np.arange(-90, 91, 1)  # 1 degree
        target_grid = Grid.from_degrees(target_lon, target_lat)

        # Regrid using conservative method (best for precipitation)
        regridded = regrid_dataset(dataset, target_grid, method="conservative")

        # Verify results
        assert regridded["pr"].shape == (12, 181, 361)  # (time, lat, lon)
        assert list(regridded["pr"].dims) == ["time", "lat", "lon"]
        # Precipitation should remain positive
        assert np.all(regridded["pr"].values >= 0)
        assert regridded.attrs["title"] == "CMIP6 precipitation data"

        # Test subsetting (like SPCZ region selection)
        spcz_region = regridded.sel(lat=slice(-25, 0), lon=slice(130, 220))
        assert spcz_region["pr"].shape[1:] == (26, 91)  # Correct spatial subset


# Performance tests (marked as slow)
@pytest.mark.slow
class TestRegridPerformance:
    """Performance tests for regridding operations."""

    def test_large_dataset_regridding(self):
        """Test regridding with large datasets."""
        # Large dataset simulation
        time = np.arange(100)  # 100 time steps
        source_lon = np.linspace(0, 360, 288, endpoint=False)  # 1.25 degree
        source_lat = np.linspace(-90, 90, 145)  # 1.25 degree

        # Create large dataset
        large_data = np.random.randn(100, 145, 288)

        dataset = xr.Dataset(
            {"temperature": (["time", "lat", "lon"], large_data)},
            coords={"time": time, "lat": source_lat, "lon": source_lon},
        )

        # Target grid
        target_lon = np.linspace(0, 360, 144, endpoint=False)
        target_lat = np.linspace(-60, 60, 73)
        target_grid = Grid.from_degrees(target_lon, target_lat)

        # Should complete without memory issues
        result = regrid_dataset(dataset, target_grid, method="conservative")

        assert result["temperature"].shape == (100, 73, 144)
        assert np.all(np.isfinite(result["temperature"].values))

    def test_repeated_regridding_performance(self):
        """Test that repeated regridding benefits from caching."""
        source_grid = Grid.from_degrees(
            np.linspace(0, 360, 72, endpoint=False), np.linspace(-60, 60, 37)
        )
        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 36, endpoint=False), np.linspace(-30, 30, 19)
        )

        regridder = ConservativeRegridder(source_grid, target_grid)

        # First regridding computes weights
        data1 = np.random.randn(*source_grid.shape)
        result1 = regridder.regrid_array(data1)

        # Second regridding should be faster (uses cached weights)
        data2 = np.random.randn(*source_grid.shape)
        result2 = regridder.regrid_array(data2)

        assert result1.shape == result2.shape == target_grid.shape
        # Both should produce valid results
        assert np.all(np.isfinite(result1))
        assert np.all(np.isfinite(result2))


if __name__ == "__main__":
    # Quick test runner for development
    pytest.main([__file__, "-v"])
