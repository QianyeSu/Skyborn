"""
Comprehensive tests for skyborn.calc.geostrophic.interface module.

This test suite provides >95% code coverage for the geostrophic wind interface module,
testing all functions, methods, and edge cases for geostrophic wind calculations.
"""

import sys
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

# Mock the geostrophic module before importing with specific function behavior
mock_geostrophic = MagicMock()


def mock_z2geouv(z, zmsg=-999.0, glon=None, glat=None, iopt=1):
    """Mock z2geouv function."""
    nlat, nlon = z.shape
    ug = np.random.rand(nlat, nlon).astype(np.float32) * 10
    vg = np.random.rand(nlat, nlon).astype(np.float32) * 10
    return ug, vg


def mock_z2geouv_3d(z, zmsg=-999.0, glon=None, glat=None, iopt=1):
    """Mock z2geouv_3d function."""
    nlat, nlon, nother = z.shape
    ug = np.random.rand(nlat, nlon, nother).astype(np.float32) * 10
    vg = np.random.rand(nlat, nlon, nother).astype(np.float32) * 10
    return ug, vg


mock_geostrophic.z2geouv = mock_z2geouv
mock_geostrophic.z2geouv_3d = mock_z2geouv_3d
sys.modules["geostrophicwind"] = mock_geostrophic


from skyborn.calc.geostrophic.interface import (
    GeostrophicWind,
    _ensure_south_to_north,
    _is_longitude_cyclic,
    geostrophic_speed,
    geostrophic_uv,
    geostrophic_wind,
)


class TestIsLongitudeCyclic:
    """Test the _is_longitude_cyclic function."""

    def test_standard_global_grid_cyclic(self):
        """Test standard global longitude grid (0 to 357.5)."""
        lon = np.arange(0, 360, 2.5)  # 0, 2.5, ..., 357.5
        assert _is_longitude_cyclic(lon) == True

    def test_full_360_grid_cyclic(self):
        """Test longitude grid that includes both 0 and 360."""
        lon = np.arange(0, 361, 10)  # 0, 10, ..., 360
        assert _is_longitude_cyclic(lon) == True

    def test_regional_grid_not_cyclic(self):
        """Test regional longitude grid (not global)."""
        lon = np.arange(-30, 30, 2.5)  # -30 to 27.5
        assert _is_longitude_cyclic(lon) == False

    def test_coarse_grid_cyclic(self):
        """Test coarse resolution global grid."""
        lon = np.array([0, 90, 180, 270])  # 90-degree spacing
        assert _is_longitude_cyclic(lon) == True

    def test_fine_grid_cyclic(self):
        """Test fine resolution global grid."""
        lon = np.arange(0, 360, 0.25)  # 0.25-degree spacing
        assert _is_longitude_cyclic(lon) == True

    def test_irregular_spacing_not_cyclic(self):
        """Test irregular longitude spacing."""
        lon = np.array([0, 10, 30, 100, 200])
        assert _is_longitude_cyclic(lon) == False

    def test_too_few_points_not_cyclic(self):
        """Test that grids with <3 points are not cyclic."""
        lon = np.array([0, 180])  # Only 2 points
        assert _is_longitude_cyclic(lon) == False

        lon = np.array([0])  # Only 1 point
        assert _is_longitude_cyclic(lon) == False

    def test_negative_longitude_range_cyclic(self):
        """Test longitude grid with negative values (e.g., -180 to 180)."""
        lon = np.arange(-180, 180, 2.5)  # -180 to 177.5
        assert _is_longitude_cyclic(lon) == True

    def test_shifted_global_grid_cyclic(self):
        """Test globally-spanning grid with different start point."""
        lon = np.arange(30, 390, 2.5)  # 30 to 387.5 degrees
        assert _is_longitude_cyclic(lon) == True

    def test_custom_tolerance(self):
        """Test custom tolerance parameter."""
        # Grid that's almost but not quite global
        lon = np.arange(0, 355, 5)  # Missing last 5 degrees

        # With default tolerance (1.0), should not be cyclic
        assert _is_longitude_cyclic(lon, tolerance=1.0) == False

        # With larger tolerance, should be cyclic
        assert _is_longitude_cyclic(lon, tolerance=10.0) == True


class TestEnsureSouthToNorth:
    """Test the _ensure_south_to_north function."""

    def setup_method(self):
        """Set up test data."""
        self.nlat, self.nlon = 5, 8
        self.z_2d = np.random.randn(self.nlat, self.nlon)
        self.z_3d = np.random.randn(10, self.nlat, self.nlon)  # (time, lat, lon)

    def test_already_south_to_north_2d(self):
        """Test data already ordered south-to-north."""
        lat = np.array([-60, -30, 0, 30, 60])  # Already south-to-north
        dim_order = "yx"

        z_ordered, lat_ordered = _ensure_south_to_north(self.z_2d, lat, dim_order)

        # Should be unchanged
        np.testing.assert_array_equal(z_ordered, self.z_2d)
        np.testing.assert_array_equal(lat_ordered, lat)

    def test_north_to_south_reversal_2d(self):
        """Test reversal of north-to-south data for 2D case."""
        lat = np.array([60, 30, 0, -30, -60])  # North-to-south
        dim_order = "yx"

        z_ordered, lat_ordered = _ensure_south_to_north(self.z_2d, lat, dim_order)

        # Latitude should be reversed
        expected_lat = np.array([-60, -30, 0, 30, 60])
        np.testing.assert_array_equal(lat_ordered, expected_lat)

        # Data should be flipped along latitude axis (axis 0 in 'yx')
        expected_z = np.flip(self.z_2d, axis=0)
        np.testing.assert_array_equal(z_ordered, expected_z)

    def test_north_to_south_reversal_3d(self):
        """Test reversal for 3D data (time, lat, lon)."""
        lat = np.array([60, 30, 0, -30, -60])  # North-to-south
        dim_order = "tyx"  # time, lat, lon

        z_ordered, lat_ordered = _ensure_south_to_north(self.z_3d, lat, dim_order)

        # Latitude should be reversed
        expected_lat = np.array([-60, -30, 0, 30, 60])
        np.testing.assert_array_equal(lat_ordered, expected_lat)

        # Data should be flipped along latitude axis (axis 1 in 'tyx')
        expected_z = np.flip(self.z_3d, axis=1)
        np.testing.assert_array_equal(z_ordered, expected_z)

    def test_same_latitude_values(self):
        """Test handling of same latitude values (edge case)."""
        lat = np.array([30, 30, 30, 30, 30])  # All same
        dim_order = "yx"

        z_ordered, lat_ordered = _ensure_south_to_north(self.z_2d, lat, dim_order)

        # Should be unchanged since lat[0] == lat[-1]
        np.testing.assert_array_equal(z_ordered, self.z_2d)
        np.testing.assert_array_equal(lat_ordered, lat)

    def test_invalid_dimension_order(self):
        """Test error for invalid dimension order (no 'y')."""
        lat = np.array([60, 30, 0, -30, -60])
        dim_order = "tx"  # No 'y' for latitude

        with pytest.raises(ValueError, match="Latitude dimension 'y' not found"):
            _ensure_south_to_north(self.z_2d, lat, dim_order)

    def test_different_dimension_orders(self):
        """Test different dimension orders."""
        lat = np.array([60, 30, 0, -30, -60])  # North-to-south

        # Test 'xy' order (lon, lat)
        z_xy = np.random.randn(self.nlon, self.nlat)
        z_ordered, lat_ordered = _ensure_south_to_north(z_xy, lat, "xy")
        expected_z = np.flip(z_xy, axis=1)  # Flip along lat axis (axis 1)
        np.testing.assert_array_equal(z_ordered, expected_z)

        # Test 'zyx' order (level, lat, lon)
        z_zyx = np.random.randn(7, self.nlat, self.nlon)
        z_ordered, lat_ordered = _ensure_south_to_north(z_zyx, lat, "zyx")
        expected_z = np.flip(z_zyx, axis=1)  # Flip along lat axis (axis 1)
        np.testing.assert_array_equal(z_ordered, expected_z)

    def test_complex_4d_case(self):
        """Test 4D case with multiple dimensions."""
        lat = np.array([60, 30, 0, -30, -60])
        z_4d = np.random.randn(12, 17, self.nlat, self.nlon)  # (time, level, lat, lon)
        dim_order = "tzyx"

        z_ordered, lat_ordered = _ensure_south_to_north(z_4d, lat, dim_order)

        # Should flip along axis 2 ('y' is at position 2 in 'tzyx')
        expected_z = np.flip(z_4d, axis=2)
        np.testing.assert_array_equal(z_ordered, expected_z)


class TestGeostrophicWindFunction:
    """Test the main geostrophic_wind function."""

    def setup_method(self):
        """Set up test data."""
        self.nlat, self.nlon = 73, 144
        self.lat = np.linspace(-90, 90, self.nlat)
        self.lon = np.linspace(0, 357.5, self.nlon)
        self.z_2d = np.random.randn(self.nlat, self.nlon) * 100 + 5500

        # 3D data (time, lat, lon)
        self.ntime = 12
        self.z_3d = np.random.randn(self.ntime, self.nlat, self.nlon) * 100 + 5500

        # 4D data (time, level, lat, lon)
        self.nlev = 17
        self.z_4d = (
            np.random.randn(self.ntime, self.nlev, self.nlat, self.nlon) * 100 + 5500
        )

    def test_2d_geostrophic_wind_cyclic(self):
        """Test 2D geostrophic wind calculation with cyclic longitude."""
        ug, vg = geostrophic_wind(self.z_2d, self.lon, self.lat, "yx")

        # Check return values have correct shape
        assert ug.shape == (self.nlat, self.nlon)
        assert vg.shape == (self.nlat, self.nlon)
        assert ug.dtype == np.float32
        assert vg.dtype == np.float32

    def test_2d_geostrophic_wind_non_cyclic(self):
        """Test 2D geostrophic wind calculation with non-cyclic longitude."""
        # Regional longitude grid
        regional_lon = np.linspace(-30, 30, 60)
        z_regional = np.random.randn(self.nlat, 60) * 100 + 5500

        ug, vg = geostrophic_wind(z_regional, regional_lon, self.lat, "yx")

        # Check return values have correct shape
        assert ug.shape == (self.nlat, 60)
        assert vg.shape == (self.nlat, 60)

    @patch("skyborn.windspharm.tools.prep_data")
    @patch("skyborn.windspharm.tools.recover_data")
    @patch("geostrophicwind.z2geouv_3d")
    def test_3d_geostrophic_wind(self, mock_fortran, mock_recover, mock_prep):
        """Test 3D geostrophic wind calculation."""
        # Mock prep_data return (nlat, nlon, combined_other)
        prepared_z = np.random.randn(self.nlat, self.nlon, self.ntime).astype(
            np.float32
        )
        recovery_info = {"test": "recovery_info"}
        mock_prep.return_value = (prepared_z, recovery_info)

        # Mock Fortran return
        mock_ug = np.random.randn(self.nlat, self.nlon, self.ntime).astype(np.float32)
        mock_vg = np.random.randn(self.nlat, self.nlon, self.ntime).astype(np.float32)
        mock_fortran.return_value = (mock_ug, mock_vg)

        # Mock recover_data
        final_ug = np.random.randn(self.ntime, self.nlat, self.nlon)
        final_vg = np.random.randn(self.ntime, self.nlat, self.nlon)
        mock_recover.side_effect = [final_ug, final_vg]

        ug, vg = geostrophic_wind(self.z_3d, self.lon, self.lat, "tyx")

        # Check that prep_data was called
        mock_prep.assert_called_once_with(self.z_3d, "tyx")

        # Check that Fortran function was called with prepared data
        mock_fortran.assert_called_once()
        call_args = mock_fortran.call_args[1]
        np.testing.assert_array_equal(call_args["z"], prepared_z)

        # Check that recover_data was called
        assert mock_recover.call_count == 2

        # Check return values
        np.testing.assert_array_equal(ug, final_ug)
        np.testing.assert_array_equal(vg, final_vg)

    def test_custom_missing_value(self):
        """Test custom missing value parameter."""
        with patch("geostrophicwind.z2geouv") as mock_fortran:
            mock_fortran.return_value = (
                np.random.randn(self.nlat, self.nlon),
                np.random.randn(self.nlat, self.nlon),
            )

            geostrophic_wind(self.z_2d, self.lon, self.lat, "yx", missing_value=-888.0)

            call_args = mock_fortran.call_args[1]
            assert call_args["zmsg"] == -888.0

    def test_latitude_ordering_correction(self):
        """Test that north-to-south latitude is corrected."""
        # North-to-south latitude
        lat_n2s = self.lat[::-1]
        z_n2s = self.z_2d[::-1, :]  # Flip data to match

        with patch("geostrophicwind.z2geouv") as mock_fortran:
            mock_fortran.return_value = (
                np.random.randn(self.nlat, self.nlon),
                np.random.randn(self.nlat, self.nlon),
            )

            geostrophic_wind(z_n2s, self.lon, lat_n2s, "yx")

            # Check that latitude was corrected to south-to-north
            call_args = mock_fortran.call_args[1]
            passed_lat = call_args["glat"]
            assert passed_lat[0] < passed_lat[-1]  # Should be south-to-north

    @patch("skyborn.windspharm.tools.prep_data")
    def test_coordinate_validation_errors(self, mock_prep):
        """Test coordinate validation errors."""
        mock_prep.return_value = (
            np.random.randn(self.nlat + 1, self.nlon, self.ntime).astype(np.float32),
            {},
        )

        with pytest.raises(ValueError, match="Latitude array length"):
            geostrophic_wind(self.z_3d, self.lon, self.lat, "tyx")

        mock_prep.return_value = (
            np.random.randn(self.nlat, self.nlon + 1, self.ntime).astype(np.float32),
            {},
        )

        with pytest.raises(ValueError, match="Longitude array length"):
            geostrophic_wind(self.z_3d, self.lon, self.lat, "tyx")

    def test_data_type_conversion(self):
        """Test that input data is converted to float32."""
        z_int = self.z_2d.astype(np.int32)

        with patch("geostrophicwind.z2geouv") as mock_fortran:
            mock_fortran.return_value = (
                np.random.randn(self.nlat, self.nlon),
                np.random.randn(self.nlat, self.nlon),
            )

            geostrophic_wind(z_int, self.lon, self.lat, "yx")

            # Check that z was converted to float32
            call_args = mock_fortran.call_args[0]
            z_passed = call_args[0]
            assert z_passed.dtype == np.float32


class TestGeostrophicWindClass:
    """Test the GeostrophicWind class."""

    def setup_method(self):
        """Set up test data."""
        self.nlat, self.nlon = 73, 144
        self.lat = np.linspace(-90, 90, self.nlat)
        self.lon = np.linspace(0, 357.5, self.nlon)
        self.z = np.random.randn(self.nlat, self.nlon) * 100 + 5500

    def test_class_initialization(self):
        """Test class initialization and property access."""
        with patch("skyborn.calc.geostrophic.interface.geostrophic_wind") as mock_func:
            mock_ug = np.random.randn(self.nlat, self.nlon) * 10
            mock_vg = np.random.randn(self.nlat, self.nlon) * 10
            mock_func.return_value = (mock_ug, mock_vg)

            gw = GeostrophicWind(self.z, self.lon, self.lat, "yx")

            # Check properties
            np.testing.assert_array_equal(gw.geopotential_height, self.z)
            np.testing.assert_array_equal(gw.longitude, self.lon)
            np.testing.assert_array_equal(gw.latitude, self.lat)

            # Check that geostrophic_wind was called during initialization
            mock_func.assert_called_once_with(
                self.z, self.lon, self.lat, "yx", missing_value=-999.0
            )

    def test_uv_components_method(self):
        """Test uv_components method."""
        with patch("skyborn.calc.geostrophic.interface.geostrophic_wind") as mock_func:
            mock_ug = np.random.randn(self.nlat, self.nlon) * 10
            mock_vg = np.random.randn(self.nlat, self.nlon) * 10
            mock_func.return_value = (mock_ug, mock_vg)

            gw = GeostrophicWind(self.z, self.lon, self.lat, "yx")
            ug, vg = gw.uv_components()

            np.testing.assert_array_equal(ug, mock_ug)
            np.testing.assert_array_equal(vg, mock_vg)

    def test_speed_calculation(self):
        """Test wind speed calculation."""
        with patch("skyborn.calc.geostrophic.interface.geostrophic_wind") as mock_func:
            mock_ug = np.array([[3.0, 4.0], [0.0, -5.0]])
            mock_vg = np.array([[4.0, 3.0], [12.0, 12.0]])
            mock_func.return_value = (mock_ug, mock_vg)

            gw = GeostrophicWind(
                np.random.randn(2, 2), np.array([0, 180]), np.array([-45, 45]), "yx"
            )
            speed = gw.speed()

            # Calculate expected speed: sqrt(ug^2 + vg^2)
            expected_speed = np.array([[5.0, 5.0], [12.0, 13.0]])
            np.testing.assert_array_almost_equal(speed, expected_speed)

    def test_speed_with_missing_values(self):
        """Test speed calculation with missing values."""
        with patch("skyborn.calc.geostrophic.interface.geostrophic_wind") as mock_func:
            mock_ug = np.array([[3.0, -999.0], [0.0, 4.0]])
            mock_vg = np.array([[4.0, 5.0], [-999.0, 3.0]])
            mock_func.return_value = (mock_ug, mock_vg)

            gw = GeostrophicWind(
                np.random.randn(2, 2),
                np.array([0, 180]),
                np.array([-45, 45]),
                "yx",
                missing_value=-999.0,
            )
            speed = gw.speed()

            expected_speed = np.array([[5.0, -999.0], [-999.0, 5.0]])
            np.testing.assert_array_equal(speed, expected_speed)

    def test_custom_missing_value(self):
        """Test class with custom missing value."""
        with patch("skyborn.calc.geostrophic.interface.geostrophic_wind") as mock_func:
            mock_ug = np.random.randn(self.nlat, self.nlon) * 10
            mock_vg = np.random.randn(self.nlat, self.nlon) * 10
            mock_func.return_value = (mock_ug, mock_vg)

            gw = GeostrophicWind(self.z, self.lon, self.lat, "yx", missing_value=-888.0)

            # Check that custom missing value was passed
            mock_func.assert_called_once_with(
                self.z, self.lon, self.lat, "yx", missing_value=-888.0
            )


class TestConvenienceFunctions:
    """Test convenience functions geostrophic_uv and geostrophic_speed."""

    def setup_method(self):
        """Set up test data."""
        self.nlat, self.nlon = 73, 144
        self.lat = np.linspace(-90, 90, self.nlat)
        self.lon = np.linspace(0, 357.5, self.nlon)
        self.z = np.random.randn(self.nlat, self.nlon) * 100 + 5500

    @patch("skyborn.calc.geostrophic.interface.GeostrophicWind")
    def test_geostrophic_uv_function(self, mock_class):
        """Test geostrophic_uv convenience function."""
        # Mock GeostrophicWind instance
        mock_instance = MagicMock()
        mock_ug = np.random.randn(self.nlat, self.nlon) * 10
        mock_vg = np.random.randn(self.nlat, self.nlon) * 10
        mock_instance.uv_components.return_value = (mock_ug, mock_vg)
        mock_class.return_value = mock_instance

        ug, vg = geostrophic_uv(self.z, self.lon, self.lat, "yx", missing_value=-888.0)

        # Check that GeostrophicWind was created with correct parameters
        mock_class.assert_called_once_with(
            self.z, self.lon, self.lat, "yx", missing_value=-888.0
        )

        # Check return values
        np.testing.assert_array_equal(ug, mock_ug)
        np.testing.assert_array_equal(vg, mock_vg)

    @patch("skyborn.calc.geostrophic.interface.GeostrophicWind")
    def test_geostrophic_speed_function(self, mock_class):
        """Test geostrophic_speed convenience function."""
        # Mock GeostrophicWind instance
        mock_instance = MagicMock()
        mock_speed = np.random.randn(self.nlat, self.nlon) * 10
        mock_instance.speed.return_value = mock_speed
        mock_class.return_value = mock_instance

        speed = geostrophic_speed(
            self.z, self.lon, self.lat, "yx", missing_value=-777.0
        )

        # Check that GeostrophicWind was created with correct parameters
        mock_class.assert_called_once_with(
            self.z, self.lon, self.lat, "yx", missing_value=-777.0
        )

        # Check return value
        np.testing.assert_array_equal(speed, mock_speed)


class TestEdgeCasesAndErrorHandling:
    """Test edge cases and error handling."""

    def setup_method(self):
        """Set up test data."""
        self.lat = np.array([-45, 0, 45])
        self.lon = np.array([0, 120, 240])
        self.z_small = np.random.randn(3, 3) * 100 + 5500

    def test_minimum_grid_size(self):
        """Test with minimum viable grid size."""
        with patch("geostrophicwind.z2geouv") as mock_fortran:
            mock_fortran.return_value = (
                np.random.randn(3, 3),
                np.random.randn(3, 3),
            )

            ug, vg = geostrophic_wind(self.z_small, self.lon, self.lat, "yx")

            # Should work without errors
            assert ug.shape == (3, 3)
            assert vg.shape == (3, 3)

    def test_single_point_longitude(self):
        """Test handling of very small longitude arrays."""
        single_lon = np.array([0])
        assert _is_longitude_cyclic(single_lon) == False

        two_point_lon = np.array([0, 180])
        assert _is_longitude_cyclic(two_point_lon) == False

    def test_large_tolerance_cyclic_detection(self):
        """Test cyclic detection with very large tolerance."""
        # Grid that spans only 270 degrees
        partial_lon = np.arange(0, 270, 10)

        # With very large tolerance, should be considered cyclic
        assert _is_longitude_cyclic(partial_lon, tolerance=100.0) == True

    def test_extremely_fine_grid(self):
        """Test cyclic detection with very fine resolution."""
        fine_lon = np.arange(0, 360, 0.1)  # 0.1-degree resolution
        assert _is_longitude_cyclic(fine_lon) == True

    def test_non_monotonic_latitude_handling(self):
        """Test handling of non-monotonic latitude arrays."""
        non_monotonic_lat = np.array([0, -30, 30, 60, -60])  # Non-monotonic
        z_test = np.random.randn(5, 3) * 100 + 5500

        # Should not crash, but behavior depends on implementation
        # The function checks only first and last elements
        z_ordered, lat_ordered = _ensure_south_to_north(z_test, non_monotonic_lat, "yx")

        # Should work (only checks first vs last element)
        assert z_ordered.shape == z_test.shape
        assert lat_ordered.shape == non_monotonic_lat.shape


class TestDataTypeHandling:
    """Test different data types and precision handling."""

    def setup_method(self):
        """Set up test data."""
        self.nlat, self.nlon = 10, 15
        self.lat = np.linspace(-90, 90, self.nlat)
        self.lon = np.linspace(0, 357.5, self.nlon)

    def test_float64_input_conversion(self):
        """Test that float64 input is converted to float32."""
        z_float64 = (
            np.random.randn(self.nlat, self.nlon).astype(np.float64) * 100 + 5500
        )

        with patch("geostrophicwind.z2geouv") as mock_fortran:
            mock_fortran.return_value = (
                np.random.randn(self.nlat, self.nlon).astype(np.float32),
                np.random.randn(self.nlat, self.nlon).astype(np.float32),
            )

            geostrophic_wind(z_float64, self.lon, self.lat, "yx")

            # Check that data was converted to float32
            call_args = mock_fortran.call_args[0]
            z_passed = call_args[0]
            assert z_passed.dtype == np.float32

    def test_integer_input_conversion(self):
        """Test that integer input is converted to float32."""
        z_int = (np.random.randn(self.nlat, self.nlon) * 100 + 5500).astype(np.int32)

        with patch("geostrophicwind.z2geouv") as mock_fortran:
            mock_fortran.return_value = (
                np.random.randn(self.nlat, self.nlon).astype(np.float32),
                np.random.randn(self.nlat, self.nlon).astype(np.float32),
            )

            geostrophic_wind(z_int, self.lon, self.lat, "yx")

            # Check that data was converted to float32
            call_args = mock_fortran.call_args[0]
            z_passed = call_args[0]
            assert z_passed.dtype == np.float32

    def test_float32_input_preserved(self):
        """Test that float32 input is preserved."""
        z_float32 = (
            np.random.randn(self.nlat, self.nlon).astype(np.float32) * 100 + 5500
        )

        with patch("geostrophicwind.z2geouv") as mock_fortran:
            mock_fortran.return_value = (
                np.random.randn(self.nlat, self.nlon).astype(np.float32),
                np.random.randn(self.nlat, self.nlon).astype(np.float32),
            )

            geostrophic_wind(z_float32, self.lon, self.lat, "yx")

            # Check that data type was preserved
            call_args = mock_fortran.call_args[0]
            z_passed = call_args[0]
            assert z_passed.dtype == np.float32


class TestComplexDimensionOrders:
    """Test various dimension orders and arrangements."""

    def setup_method(self):
        """Set up test data for different dimension orders."""
        self.nlat, self.nlon, self.ntime, self.nlev = 18, 36, 12, 17
        self.lat = np.linspace(-90, 90, self.nlat)
        self.lon = np.linspace(0, 357.5, self.nlon)

    @patch("skyborn.windspharm.tools.prep_data")
    @patch("skyborn.windspharm.tools.recover_data")
    @patch("geostrophicwind.z2geouv_3d")
    def test_tzyx_dimension_order(self, mock_fortran, mock_recover, mock_prep):
        """Test (time, level, lat, lon) dimension order."""
        z_4d = np.random.randn(self.ntime, self.nlev, self.nlat, self.nlon) * 100 + 5500

        # Mock prep_data return
        prepared_z = np.random.randn(
            self.nlat, self.nlon, self.ntime * self.nlev
        ).astype(np.float32)
        mock_prep.return_value = (prepared_z, {"test": "info"})

        # Mock Fortran return
        mock_ug = np.random.randn(self.nlat, self.nlon, self.ntime * self.nlev).astype(
            np.float32
        )
        mock_vg = np.random.randn(self.nlat, self.nlon, self.ntime * self.nlev).astype(
            np.float32
        )
        mock_fortran.return_value = (mock_ug, mock_vg)

        # Mock recover_data
        final_ug = np.random.randn(self.ntime, self.nlev, self.nlat, self.nlon)
        final_vg = np.random.randn(self.ntime, self.nlev, self.nlat, self.nlon)
        mock_recover.side_effect = [final_ug, final_vg]

        ug, vg = geostrophic_wind(z_4d, self.lon, self.lat, "tzyx")

        # Check that prep_data was called with correct dimension order
        mock_prep.assert_called_once_with(z_4d, "tzyx")

        # Check return shapes
        assert ug.shape == (self.ntime, self.nlev, self.nlat, self.nlon)
        assert vg.shape == (self.ntime, self.nlev, self.nlat, self.nlon)

    @patch("skyborn.windspharm.tools.prep_data")
    @patch("skyborn.windspharm.tools.recover_data")
    @patch("geostrophicwind.z2geouv_3d")
    def test_yxzt_dimension_order(self, mock_fortran, mock_recover, mock_prep):
        """Test (lat, lon, level, time) dimension order."""
        z_4d = np.random.randn(self.nlat, self.nlon, self.nlev, self.ntime) * 100 + 5500

        # Mock prep_data return
        prepared_z = np.random.randn(
            self.nlat, self.nlon, self.nlev * self.ntime
        ).astype(np.float32)
        mock_prep.return_value = (prepared_z, {"test": "info"})

        # Mock Fortran return
        mock_ug = np.random.randn(self.nlat, self.nlon, self.nlev * self.ntime).astype(
            np.float32
        )
        mock_vg = np.random.randn(self.nlat, self.nlon, self.nlev * self.ntime).astype(
            np.float32
        )
        mock_fortran.return_value = (mock_ug, mock_vg)

        # Mock recover_data
        final_ug = np.random.randn(self.nlat, self.nlon, self.nlev, self.ntime)
        final_vg = np.random.randn(self.nlat, self.nlon, self.nlev, self.ntime)
        mock_recover.side_effect = [final_ug, final_vg]

        ug, vg = geostrophic_wind(z_4d, self.lon, self.lat, "yxzt")

        # Check that prep_data was called with correct dimension order
        mock_prep.assert_called_once_with(z_4d, "yxzt")

        # Check return shapes
        assert ug.shape == (self.nlat, self.nlon, self.nlev, self.ntime)
        assert vg.shape == (self.nlat, self.nlon, self.nlev, self.ntime)


class TestRegressionTests:
    """Test specific scenarios that might have caused issues in the past."""

    def setup_method(self):
        """Set up regression test data."""
        self.lat = np.array([-90, -45, 0, 45, 90])
        self.lon = np.array([0, 90, 180, 270])  # Exactly 90-degree spacing
        self.z = np.random.randn(5, 4) * 100 + 5500

    def test_90_degree_longitude_spacing(self):
        """Test longitude cyclicity with 90-degree spacing."""
        # This is a coarse but global grid
        assert _is_longitude_cyclic(self.lon) == True

    def test_poles_included_latitude(self):
        """Test handling when poles are included in latitude array."""
        # Includes -90 and +90
        with patch("geostrophicwind.z2geouv") as mock_fortran:
            mock_fortran.return_value = (
                np.random.randn(5, 4),
                np.random.randn(5, 4),
            )

            # Should work without issues
            ug, vg = geostrophic_wind(self.z, self.lon, self.lat, "yx")

            assert ug.shape == (5, 4)
            assert vg.shape == (5, 4)

    def test_zero_geopotential_height(self):
        """Test handling of zero or very small geopotential heights."""
        z_zero = np.zeros((5, 4))

        with patch("geostrophicwind.z2geouv") as mock_fortran:
            mock_fortran.return_value = (
                np.random.randn(5, 4),
                np.random.randn(5, 4),
            )

            # Should handle zero heights
            ug, vg = geostrophic_wind(z_zero, self.lon, self.lat, "yx")

            # Check that zero array was passed
            call_args = mock_fortran.call_args[0]
            z_passed = call_args[0]
            np.testing.assert_array_equal(z_passed, z_zero.astype(np.float32))

    def test_very_large_geopotential_height(self):
        """Test handling of very large geopotential heights."""
        z_large = np.ones((5, 4)) * 1e6  # Very large values

        with patch("geostrophicwind.z2geouv") as mock_fortran:
            mock_fortran.return_value = (
                np.random.randn(5, 4),
                np.random.randn(5, 4),
            )

            # Should handle large values
            ug, vg = geostrophic_wind(z_large, self.lon, self.lat, "yx")

            # Check that large values were passed correctly
            call_args = mock_fortran.call_args[0]
            z_passed = call_args[0]
            assert np.allclose(z_passed, z_large.astype(np.float32))


class TestExtraCodeCoverageTargets:
    """Additional tests to target specific code lines for >95% coverage."""

    def setup_method(self):
        """Set up test data."""
        self.lat = np.linspace(-90, 90, 10)
        self.lon = np.linspace(0, 357.5, 15)
        self.z = np.random.randn(10, 15) * 100 + 5500

    def test_all_class_properties_accessed(self):
        """Test that all class properties are accessible."""
        with patch("skyborn.calc.geostrophic.interface.geostrophic_wind") as mock_func:
            mock_ug = np.random.randn(10, 15) * 10
            mock_vg = np.random.randn(10, 15) * 10
            mock_func.return_value = (mock_ug, mock_vg)

            gw = GeostrophicWind(self.z, self.lon, self.lat, "yx")

            # Access all properties to ensure coverage
            _ = gw.geopotential_height
            _ = gw.longitude
            _ = gw.latitude

            # Test methods
            ug, vg = gw.uv_components()
            speed = gw.speed()

            # Verify all work correctly
            assert ug.shape == (10, 15)
            assert vg.shape == (10, 15)
            assert speed.shape == (10, 15)

    def test_speed_all_missing_values(self):
        """Test speed calculation when all values are missing."""
        with patch("skyborn.calc.geostrophic.interface.geostrophic_wind") as mock_func:
            mock_ug = np.full((2, 2), -999.0)
            mock_vg = np.full((2, 2), -999.0)
            mock_func.return_value = (mock_ug, mock_vg)

            gw = GeostrophicWind(
                np.random.randn(2, 2),
                np.array([0, 180]),
                np.array([-45, 45]),
                "yx",
                missing_value=-999.0,
            )
            speed = gw.speed()

            # All values should be missing
            expected_speed = np.full((2, 2), -999.0)
            np.testing.assert_array_equal(speed, expected_speed)

    def test_mixed_missing_values_in_speed(self):
        """Test speed calculation with mixed missing and valid values."""
        with patch("skyborn.calc.geostrophic.interface.geostrophic_wind") as mock_func:
            # Mixed scenario: some valid, some missing in different components
            mock_ug = np.array([[3.0, -999.0], [0.0, 4.0]])
            mock_vg = np.array([[4.0, 3.0], [-999.0, -999.0]])
            mock_func.return_value = (mock_ug, mock_vg)

            gw = GeostrophicWind(
                np.random.randn(2, 2),
                np.array([0, 180]),
                np.array([-45, 45]),
                "yx",
                missing_value=-999.0,
            )
            speed = gw.speed()

            # Expected: valid where both components valid, missing otherwise
            expected_speed = np.array([[5.0, -999.0], [-999.0, -999.0]])
            np.testing.assert_array_equal(speed, expected_speed)

    def test_longitude_range_edge_case_values(self):
        """Test longitude cyclicity detection with edge case values."""
        # Test longitude range that's very close to 360 but not quite
        almost_global = np.arange(0, 359.9, 1.0)  # 359 degrees
        assert _is_longitude_cyclic(almost_global, tolerance=0.5) == False
        assert _is_longitude_cyclic(almost_global, tolerance=1.5) == True

        # Test longitude range that exceeds 360 slightly
        over_global = np.arange(0, 361, 1.0)  # 361 degrees
        assert _is_longitude_cyclic(over_global) == True


class TestCoverageCompletion:
    """Tests specifically to achieve 100% coverage of interface.py."""

    def test_coordinate_length_mismatch_latitude(self):
        """Test ValueError when latitude array length doesn't match data shape."""
        z_3d = np.random.rand(12, 10, 15) * 100 + 5500
        glon = np.linspace(0, 360, 15)
        glat_wrong = np.linspace(-90, 90, 8)  # Wrong length - should be 10

        with pytest.raises(
            ValueError, match="Latitude array length.*doesn't match data"
        ):
            geostrophic_wind(z_3d, glon, glat_wrong, "tyx")

    def test_coordinate_length_mismatch_longitude(self):
        """Test ValueError when longitude array length doesn't match data shape."""
        z_3d = np.random.rand(12, 10, 15) * 100 + 5500
        glon_wrong = np.linspace(0, 360, 12)  # Wrong length - should be 15
        glat = np.linspace(-90, 90, 10)

        with pytest.raises(
            ValueError, match="Longitude array length.*doesn't match data"
        ):
            geostrophic_wind(z_3d, glon_wrong, glat, "tyx")

    def test_multidimensional_data_recovery(self):
        """Test multidimensional data processing to cover recovery code paths."""
        # Test 3D data to ensure recover_data is called
        ntime, nlat, nlon = 6, 8, 12
        z_3d = np.random.rand(ntime, nlat, nlon) * 100 + 5500
        glon = np.linspace(0, 360, nlon)
        glat = np.linspace(-90, 90, nlat)

        ug, vg = geostrophic_wind(z_3d, glon, glat, "tyx")

        # Check that data recovery worked correctly
        assert ug.shape == (ntime, nlat, nlon)
        assert vg.shape == (ntime, nlat, nlon)
        assert ug.dtype == np.float32
        assert vg.dtype == np.float32

    def test_4d_data_recovery(self):
        """Test 4D data processing to ensure all recovery paths are covered."""
        ntime, nlev, nlat, nlon = 3, 5, 8, 12
        z_4d = np.random.rand(ntime, nlev, nlat, nlon) * 100 + 5500
        glon = np.linspace(0, 360, nlon)
        glat = np.linspace(-90, 90, nlat)

        ug, vg = geostrophic_wind(z_4d, glon, glat, "tzyx")

        # Check that 4D data recovery worked correctly
        assert ug.shape == (ntime, nlev, nlat, nlon)
        assert vg.shape == (ntime, nlev, nlat, nlon)
        assert ug.dtype == np.float32
        assert vg.dtype == np.float32

    def test_geostrophic_wind_class_basic(self):
        """Test GeostrophicWind class basic functionality."""
        nlat, nlon = 8, 12
        z = np.random.rand(nlat, nlon) * 100 + 5500
        glon = np.linspace(0, 360, nlon)
        glat = np.linspace(-90, 90, nlat)

        gw = GeostrophicWind(z, glon, glat, "yx")

        # Test properties
        np.testing.assert_array_equal(gw.geopotential_height, z)
        np.testing.assert_array_equal(gw.longitude, glon)
        np.testing.assert_array_equal(gw.latitude, glat)

        # Test methods
        ug, vg = gw.uv_components()
        assert ug.shape == (nlat, nlon)
        assert vg.shape == (nlat, nlon)

        speed = gw.speed()
        assert speed.shape == (nlat, nlon)

    def test_convenience_functions(self):
        """Test geostrophic_uv and geostrophic_speed convenience functions."""
        nlat, nlon = 8, 12
        z = np.random.rand(nlat, nlon) * 100 + 5500
        glon = np.linspace(0, 360, nlon)
        glat = np.linspace(-90, 90, nlat)

        # Test geostrophic_uv
        ug, vg = geostrophic_uv(z, glon, glat, "yx")
        assert ug.shape == (nlat, nlon)
        assert vg.shape == (nlat, nlon)

        # Test geostrophic_speed
        speed = geostrophic_speed(z, glon, glat, "yx")
        assert speed.shape == (nlat, nlon)

    def test_2d_direct_calculation(self):
        """Test 2D direct calculation path."""
        nlat, nlon = 8, 12
        z = np.random.rand(nlat, nlon) * 100 + 5500
        glon = np.linspace(0, 360, nlon)
        glat = np.linspace(-90, 90, nlat)

        ug, vg = geostrophic_wind(z, glon, glat, "yx")

        assert ug.shape == (nlat, nlon)
        assert vg.shape == (nlat, nlon)

    def test_longitude_too_few_points_coverage(self):
        """Test _is_longitude_cyclic with < 3 points to cover line 47."""
        # Test with 0 points
        assert _is_longitude_cyclic(np.array([])) == False

        # Test with 1 point
        assert _is_longitude_cyclic(np.array([0])) == False

        # Test with 2 points
        assert _is_longitude_cyclic(np.array([0, 180])) == False

    def test_north_to_south_latitude_coverage(self):
        """Test _ensure_south_to_north with north-to-south data to cover lines 90-98."""
        nlat, nlon = 5, 8
        z_test = np.random.rand(nlat, nlon) * 100 + 5500

        # Create north-to-south latitude array
        lat_n2s = np.array([60, 30, 0, -30, -60])  # North to south

        # Test the function directly
        z_ordered, lat_ordered = _ensure_south_to_north(z_test, lat_n2s, "yx")

        # Check that latitude was reversed to south-to-north
        expected_lat = np.array([-60, -30, 0, 30, 60])
        np.testing.assert_array_equal(lat_ordered, expected_lat)

        # Check that data was flipped along latitude axis
        expected_z = np.flip(z_test, axis=0)
        np.testing.assert_array_equal(z_ordered, expected_z)

    def test_north_to_south_3d_data_coverage(self):
        """Test _ensure_south_to_north with 3D north-to-south data."""
        ntime, nlat, nlon = 3, 5, 8
        z_3d = np.random.rand(ntime, nlat, nlon) * 100 + 5500

        # North-to-south latitude
        lat_n2s = np.array([60, 30, 0, -30, -60])

        z_ordered, lat_ordered = _ensure_south_to_north(z_3d, lat_n2s, "tyx")

        # Check latitude reversal
        expected_lat = np.array([-60, -30, 0, 30, 60])
        np.testing.assert_array_equal(lat_ordered, expected_lat)

        # Check data flip along latitude axis (axis 1 in 'tyx')
        expected_z = np.flip(z_3d, axis=1)
        np.testing.assert_array_equal(z_ordered, expected_z)

    def test_invalid_dimension_order_coverage(self):
        """Test _ensure_south_to_north with invalid dimension order to cover line 95."""
        z_test = np.random.rand(5, 8)
        lat_n2s = np.array([60, 30, 0, -30, -60])  # North to south

        # Use dimension order without 'y' - should raise error
        with pytest.raises(ValueError, match="Latitude dimension 'y' not found"):
            _ensure_south_to_north(z_test, lat_n2s, "tx")  # No 'y' dimension

    def test_multidimensional_recovery_lines_212_215(self):
        """Test to specifically cover lines 212-215: recover_data calls."""
        # This test specifically ensures we go through the multidimensional path
        # that calls recover_data (lines 212-213) and returns the result (line 215)
        ntime, nlat, nlon = 2, 4, 6
        z_3d = np.random.rand(ntime, nlat, nlon) * 100 + 5500
        glon = np.linspace(0, 360, nlon)
        glat = np.linspace(-90, 90, nlat)

        # Use 3D data to force the multidimensional path
        ug, vg = geostrophic_wind(z_3d, glon, glat, "tyx")

        # Verify the recovery worked
        assert ug.shape == (ntime, nlat, nlon)
        assert vg.shape == (ntime, nlat, nlon)

        # This should hit lines 212-215 in the _geostrophic_wind_multidim function
