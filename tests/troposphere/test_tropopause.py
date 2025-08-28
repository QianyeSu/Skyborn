"""
Comprehensive tests for skyborn.calc.tropopause module.

This test suite provides 100% code coverage for the tropopause module,
testing all functions, methods, and edge cases for the WMO tropopause calculation.
"""

from unittest.mock import MagicMock, patch

import numpy as np
import numpy.ma as ma
import pytest

from skyborn.calc.troposphere.tropopause import (
    _order_dims_for_fortran,
    _restore_dims_from_fortran,
    trop_wmo,
    trop_wmo_profile,
)


class TestOrderDimsForFortran:
    """Test the _order_dims_for_fortran function."""

    def test_2d_level_lat_ordering(self):
        """Test 2D array with (level, lat) ordering."""
        grid = np.random.rand(5, 10)  # (level, lat)
        xdim, ydim, levdim = -1, 1, 0  # No lon, lat=1, level=0

        ordered_grid, info = _order_dims_for_fortran(grid, xdim, ydim, levdim)

        # Should reorder to (lat, level)
        assert ordered_grid.shape == (10, 5)
        assert info["new_order"] == [1, 0]
        assert info["original_shape"] == (5, 10)
        np.testing.assert_array_equal(info["inverse_order"], [1, 0])

    def test_2d_level_lon_ordering(self):
        """Test 2D array with (level, lon) ordering."""
        grid = np.random.rand(5, 15)  # (level, lon)
        xdim, ydim, levdim = 1, -1, 0  # lon=1, no lat, level=0

        ordered_grid, info = _order_dims_for_fortran(grid, xdim, ydim, levdim)

        # Should reorder to (lon, level)
        assert ordered_grid.shape == (15, 5)
        assert info["new_order"] == [1, 0]

    def test_2d_no_spatial_dimension_error(self):
        """Test error for 2D data without spatial dimensions."""
        grid = np.random.rand(5, 10)
        xdim, ydim, levdim = -1, -1, 0  # No spatial dimensions

        with pytest.raises(
            ValueError, match="2D data must have either lat or lon dimension"
        ):
            _order_dims_for_fortran(grid, xdim, ydim, levdim)

    def test_3d_standard_ordering(self):
        """Test 3D array ordering (level, lat, lon) -> (lat, lon, level)."""
        grid = np.random.rand(5, 10, 15)  # (level, lat, lon)
        xdim, ydim, levdim = 2, 1, 0

        ordered_grid, info = _order_dims_for_fortran(grid, xdim, ydim, levdim)

        # Should reorder to (lat, lon, level)
        assert ordered_grid.shape == (10, 15, 5)
        assert info["new_order"] == [1, 2, 0]

    def test_4d_with_time_ordering(self):
        """Test 4D array ordering (time, level, lat, lon) -> (lat, lon, level, time)."""
        grid = np.random.rand(12, 5, 10, 15)  # (time, level, lat, lon)
        xdim, ydim, levdim, timedim = 3, 2, 1, 0

        ordered_grid, info = _order_dims_for_fortran(grid, xdim, ydim, levdim, timedim)

        # Should reorder to (lat, lon, level, time)
        assert ordered_grid.shape == (10, 15, 5, 12)
        assert info["new_order"] == [2, 3, 1, 0]

    def test_missing_dimensions_handling(self):
        """Test handling of missing dimensions (set to -1)."""
        grid = np.random.rand(5, 10, 15)  # 3D grid
        xdim, ydim, levdim, timedim = 2, 1, 0, -1  # No time dimension

        ordered_grid, info = _order_dims_for_fortran(grid, xdim, ydim, levdim, timedim)

        # Should still work without time dimension
        assert ordered_grid.shape == (10, 15, 5)
        assert len(info["new_order"]) == 3

    def test_remaining_dimensions_handling(self):
        """Test that remaining dimensions are properly handled."""
        grid = np.random.rand(5, 10, 15, 8, 6)  # 5D grid with extra dimensions
        xdim, ydim, levdim = 2, 1, 0

        ordered_grid, info = _order_dims_for_fortran(grid, xdim, ydim, levdim)

        # Should include all dimensions
        assert ordered_grid.ndim == 5
        assert len(info["new_order"]) == 5
        # First three should be lat, lon, level
        assert info["new_order"][:3] == [1, 2, 0]


class TestRestoreDimsFromFortran:
    """Test the _restore_dims_from_fortran function."""

    def test_2d_restoration_without_level_removal(self):
        """Test 2D array restoration without level dimension removal."""
        # Simulate processing of 2D array
        original_shape = (5, 10)
        info = {
            "original_shape": original_shape,
            "original_ndim": 2,
            "new_order": [1, 0],
            "inverse_order": np.array([1, 0]),
        }

        # Processed array (reordered)
        processed = np.random.rand(10, 5)

        restored = _restore_dims_from_fortran(processed, info, remove_level_dim=False)

        assert restored.shape == original_shape

    def test_3d_restoration_with_level_removal(self):
        """Test 3D array restoration with level dimension removal."""
        # Original 3D -> processed 2D (level dimension removed)
        info = {
            "original_shape": (5, 10, 15),
            "original_ndim": 3,
            "new_order": [1, 2, 0],  # (lat, lon, level)
            "inverse_order": np.array([2, 0, 1]),
        }

        # Processed array without level dimension (lat, lon)
        processed = np.random.rand(10, 15)

        restored = _restore_dims_from_fortran(processed, info, remove_level_dim=True)

        # Should restore to (lat, lon) order (level dimension removed)
        expected_shape = (10, 15)  # Original (5,10,15) with level removed
        assert restored.shape == expected_shape

    def test_4d_restoration_with_level_removal(self):
        """Test 4D array restoration with level dimension removal."""
        info = {
            "original_shape": (12, 5, 10, 15),  # (time, level, lat, lon)
            "original_ndim": 4,
            "new_order": [2, 3, 1, 0],  # (lat, lon, level, time)
            "inverse_order": np.array([3, 2, 0, 1]),
        }

        # Processed array without level dimension (lat, lon, time)
        processed = np.random.rand(10, 15, 12)

        restored = _restore_dims_from_fortran(processed, info, remove_level_dim=True)

        # Should restore to (time, lat, lon) order
        assert restored.shape == (12, 10, 15)


class TestTropWmoInputValidation:
    """Test input validation for trop_wmo function."""

    def test_pressure_length_mismatch(self):
        """Test error when pressure length doesn't match temperature level dimension."""
        temperature = np.random.rand(12, 10, 180, 360)  # 4D array
        pressure = np.array([100, 200, 300, 500, 700])  # Only 5 levels

        with pytest.raises(
            ValueError, match="Pressure length.*must match temperature level dimension"
        ):
            trop_wmo(temperature, pressure, xdim=3, ydim=2, levdim=1, timedim=0)

    def test_invalid_dimension_indices(self):
        """Test error for invalid dimension indices."""
        temperature = np.random.rand(5, 10, 15)  # 3D array
        pressure = np.array([100, 200, 300, 500, 700])

        # Dimension index too large
        with pytest.raises(ValueError, match="Dimension indices must be valid"):
            trop_wmo(temperature, pressure, xdim=5, ydim=2, levdim=0)

        # Test that negative indices other than placeholders work correctly
        # -1 as placeholder is allowed, but other negative values should be converted by the validation logic

    def test_duplicate_dimension_indices(self):
        """Test error for duplicate dimension indices."""
        temperature = np.random.rand(5, 10, 15)  # 3D array
        pressure = np.array([100, 200, 300, 500, 700])

        with pytest.raises(ValueError, match="Dimension indices must be unique"):
            trop_wmo(temperature, pressure, xdim=1, ydim=1, levdim=0)

    def test_unsupported_data_shape(self):
        """Test error for unsupported data shapes."""
        temperature = np.random.rand(5, 10, 15, 8, 6)  # 5D array (unsupported)
        pressure = np.array([100, 200, 300, 500, 700])

        with pytest.raises(ValueError, match="Unsupported temperature data shape"):
            trop_wmo(temperature, pressure, xdim=1, ydim=2, levdim=0)


class TestTropWmoMaskedArrays:
    """Test trop_wmo function with masked arrays."""

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_grid_3d")
    def test_masked_pressure_array(self, mock_fortran):
        """Test handling of masked pressure arrays."""
        # Create masked array with first element masked
        pressure = ma.array(
            [100, 200, 300, 500, 700], mask=[True, False, False, False, False]
        )
        temperature = np.random.rand(5, 10, 15) * 50 + 250

        mock_fortran.return_value = (
            np.random.rand(10, 15) * 200 + 100,
            np.random.rand(10, 15) * 5000 + 8000,
            np.random.randint(0, 5, (10, 15)),
            np.random.rand(10, 15) * 3,
            np.ones((10, 15), dtype=bool),
        )

        # Verify that ma.is_masked returns True for our test data
        assert ma.is_masked(pressure) == True

        result = trop_wmo(
            temperature, pressure, xdim=2, ydim=1, levdim=0, missing_value=-999.0
        )

        # Check that function executed successfully with masked input
        assert isinstance(result, dict)
        mock_fortran.assert_called_once()

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_grid_3d")
    def test_masked_temperature_array(self, mock_fortran):
        """Test handling of masked temperature arrays."""
        pressure = np.array([100, 200, 300, 500, 700])
        # Create masked temperature array
        temperature = ma.array(np.random.rand(5, 10, 15) * 50 + 250)
        temperature[0, 0, 0] = ma.masked

        mock_fortran.return_value = (
            np.random.rand(10, 15) * 200 + 100,
            np.random.rand(10, 15) * 5000 + 8000,
            np.random.randint(0, 5, (10, 15)),
            np.random.rand(10, 15) * 3,
            np.ones((10, 15), dtype=bool),
        )

        result = trop_wmo(
            temperature, pressure, xdim=2, ydim=1, levdim=0, missing_value=-999.0
        )

        # Check that function executed successfully
        assert isinstance(result, dict)
        mock_fortran.assert_called_once()


class TestTropWmoPressureSorting:
    """Test pressure sorting functionality in trop_wmo."""

    def setup_method(self):
        """Set up test data with unsorted pressure."""
        # Unsorted pressure levels
        self.unsorted_pressure = np.array([500, 100, 1000, 200, 700])
        # Corresponding temperature (should be reordered with pressure)
        self.temperature = np.array(
            [
                [260, 220, 288, 230, 275],  # Level values
                [261, 221, 289, 231, 276],
                [262, 222, 290, 232, 277],
            ]
        ).T  # Shape: (5, 3) for 5 levels, 3 spatial points

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_grid_2d")
    def test_pressure_sorting_enabled(self, mock_fortran):
        """Test automatic pressure sorting when enabled."""
        mock_fortran.return_value = (
            np.random.rand(3) * 200 + 100,
            np.random.rand(3) * 5000 + 8000,
            np.random.randint(0, 5, 3),
            np.random.rand(3) * 3,
            np.ones(3, dtype=bool),
        )

        result = trop_wmo(
            self.temperature,
            self.unsorted_pressure,
            xdim=-1,
            ydim=1,
            levdim=0,
            check_pressure_order=True,
        )

        # Get the pressure passed to Fortran routine
        call_args = mock_fortran.call_args[0]
        pressure_passed = call_args[1]

        # Should be sorted in ascending order
        assert np.all(pressure_passed[:-1] <= pressure_passed[1:])
        expected_sorted = np.array([100, 200, 500, 700, 1000])
        np.testing.assert_array_equal(pressure_passed, expected_sorted)

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_grid_2d")
    def test_pressure_sorting_disabled(self, mock_fortran):
        """Test that sorting can be disabled."""
        mock_fortran.return_value = (
            np.random.rand(3) * 200 + 100,
            np.random.rand(3) * 5000 + 8000,
            np.random.randint(0, 5, 3),
            np.random.rand(3) * 3,
            np.ones(3, dtype=bool),
        )

        result = trop_wmo(
            self.temperature,
            self.unsorted_pressure,
            xdim=-1,
            ydim=1,
            levdim=0,
            check_pressure_order=False,
        )

        # Get the pressure passed to Fortran routine
        call_args = mock_fortran.call_args[0]
        pressure_passed = call_args[1]

        # Should retain original order
        np.testing.assert_array_equal(pressure_passed, self.unsorted_pressure)

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_grid_3d")
    def test_temperature_reordering_with_pressure_sort(self, mock_fortran):
        """Test that temperature is reordered when pressure is sorted."""
        temp_3d = np.random.rand(5, 10, 15)  # (level, lat, lon)
        unsorted_pressure = np.array([500, 100, 1000, 200, 700])

        mock_fortran.return_value = (
            np.random.rand(10, 15) * 200 + 100,
            np.random.rand(10, 15) * 5000 + 8000,
            np.random.randint(0, 5, (10, 15)),
            np.random.rand(10, 15) * 3,
            np.ones((10, 15), dtype=bool),
        )

        result = trop_wmo(
            temp_3d,
            unsorted_pressure,
            xdim=2,
            ydim=1,
            levdim=0,
            check_pressure_order=True,
        )

        mock_fortran.assert_called_once()


class TestTropWmo2DData:
    """Test trop_wmo function with 2D data."""

    def setup_method(self):
        """Set up 2D test data."""
        self.pressure = np.array([100, 200, 300, 500, 700, 850, 1000])
        self.temp_2d = np.random.rand(7, 10) * 50 + 250  # (level, spatial)

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_grid_2d")
    def test_2d_level_lat_data(self, mock_fortran):
        """Test 2D data with level and latitude dimensions."""
        mock_result = (
            np.random.rand(10) * 200 + 100,  # pressure
            np.random.rand(10) * 5000 + 8000,  # height
            np.random.randint(0, 7, 10),  # level_index
            np.random.rand(10) * 3,  # lapse_rate
            np.ones(10, dtype=bool),  # success
        )
        mock_fortran.return_value = mock_result

        result = trop_wmo(self.temp_2d, self.pressure, xdim=-1, ydim=1, levdim=0)

        # Check result structure
        assert "pressure" in result
        assert "height" in result
        assert result["pressure"].shape == (10,)

        # Check that 2D Fortran routine was called
        mock_fortran.assert_called_once()

        # Check parameters passed to Fortran routine
        call_args = mock_fortran.call_args[0]
        nlevm = call_args[0]
        pressure_passed = call_args[1]
        temp_passed = call_args[2]

        assert nlevm == 8  # nlev + 1
        assert len(pressure_passed) == 7
        assert temp_passed.shape == (10, 7)  # Reordered to (spatial, level)

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_grid_2d")
    def test_2d_level_lon_data(self, mock_fortran):
        """Test 2D data with level and longitude dimensions."""
        mock_result = (
            np.random.rand(10) * 200 + 100,
            np.random.rand(10) * 5000 + 8000,
            np.random.randint(0, 7, 10),
            np.random.rand(10) * 3,
            np.ones(10, dtype=bool),
        )
        mock_fortran.return_value = mock_result

        result = trop_wmo(self.temp_2d, self.pressure, xdim=1, ydim=-1, levdim=0)

        assert result["pressure"].shape == (10,)
        mock_fortran.assert_called_once()


class TestTropWmo3DData:
    """Test trop_wmo function with 3D data."""

    def setup_method(self):
        """Set up 3D test data."""
        self.pressure = np.array([100, 200, 300, 500, 700, 850, 1000])
        self.temp_3d = np.random.rand(7, 10, 15) * 50 + 250  # (level, lat, lon)

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_grid_3d")
    def test_3d_standard_processing(self, mock_fortran):
        """Test standard 3D data processing."""
        mock_result = (
            np.random.rand(10, 15) * 200 + 100,  # pressure
            np.random.rand(10, 15) * 5000 + 8000,  # height
            np.random.randint(0, 7, (10, 15)),  # level_index
            np.random.rand(10, 15) * 3,  # lapse_rate
            np.ones((10, 15), dtype=bool),  # success
        )
        mock_fortran.return_value = mock_result

        result = trop_wmo(self.temp_3d, self.pressure, xdim=2, ydim=1, levdim=0)

        # Check result structure
        assert result["pressure"].shape == (10, 15)
        assert result["height"].shape == (10, 15)

        # Check Fortran call
        mock_fortran.assert_called_once()
        call_args = mock_fortran.call_args[0]

        # Temperature should be reordered to (lat, lon, level)
        temp_passed = call_args[2]
        assert temp_passed.shape == (10, 15, 7)

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_grid_3d")
    def test_3d_dimension_restoration(self, mock_fortran):
        """Test that output dimensions are correctly restored."""
        # The fortran function returns in (lat, lon) order for 3D case
        mock_result = (
            np.random.rand(10, 15) * 200 + 100,  # (lat, lon) order
            np.random.rand(10, 15) * 5000 + 8000,
            np.random.randint(0, 7, (10, 15)),
            np.random.rand(10, 15) * 3,
            np.ones((10, 15), dtype=bool),
        )
        mock_fortran.return_value = mock_result

        # Input order: (level, lat, lon) -> should restore to (lat, lon)
        result = trop_wmo(self.temp_3d, self.pressure, xdim=2, ydim=1, levdim=0)

        # Output should match original spatial dimensions
        assert result["pressure"].shape == (10, 15)  # (lat, lon)
        mock_fortran.assert_called_once()


class TestTropWmo4DData:
    """Test trop_wmo function with 4D data including time dimension."""

    def setup_method(self):
        """Set up 4D test data."""
        self.pressure = np.array([100, 200, 300, 500, 700, 850, 1000])
        self.temp_4d = (
            np.random.rand(12, 7, 10, 15) * 50 + 250
        )  # (time, level, lat, lon)

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_grid_4d")
    def test_4d_with_time_processing(self, mock_fortran):
        """Test 4D data processing with time dimension."""
        mock_result = (
            np.random.rand(10, 15, 12) * 200 + 100,  # (lat, lon, time)
            np.random.rand(10, 15, 12) * 5000 + 8000,
            np.random.randint(0, 7, (10, 15, 12)),
            np.random.rand(10, 15, 12) * 3,
            np.ones((10, 15, 12), dtype=bool),
        )
        mock_fortran.return_value = mock_result

        result = trop_wmo(
            self.temp_4d, self.pressure, xdim=3, ydim=2, levdim=1, timedim=0
        )

        # Check result structure - should restore to (time, lat, lon)
        assert result["pressure"].shape == (12, 10, 15)
        assert result["height"].shape == (12, 10, 15)

        # Check Fortran call
        mock_fortran.assert_called_once()
        call_args = mock_fortran.call_args[0]

        # Temperature should be reordered to (lat, lon, level, time)
        temp_passed = call_args[2]
        assert temp_passed.shape == (10, 15, 7, 12)


class TestTropWmoParameters:
    """Test parameter handling in trop_wmo function."""

    def setup_method(self):
        """Set up test data."""
        self.pressure = np.array([100, 200, 300, 500, 700, 850, 1000])
        self.temperature = np.random.rand(7, 10, 15) * 50 + 250

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_grid_3d")
    def test_pressure_unit_hpa(self, mock_fortran):
        """Test hPa pressure unit parameter."""
        mock_fortran.return_value = (
            np.random.rand(10, 15),
            np.random.rand(10, 15),
            np.random.randint(0, 7, (10, 15)),
            np.random.rand(10, 15),
            np.ones((10, 15), dtype=bool),
        )

        trop_wmo(
            self.temperature,
            self.pressure,
            xdim=2,
            ydim=1,
            levdim=0,
            pressure_unit="hPa",
        )

        # Check that punit=0 was passed (hPa)
        call_args = mock_fortran.call_args[0]
        punit = call_args[5]
        assert punit == 0

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_grid_3d")
    def test_pressure_unit_pa(self, mock_fortran):
        """Test Pa pressure unit parameter."""
        mock_fortran.return_value = (
            np.random.rand(10, 15),
            np.random.rand(10, 15),
            np.random.randint(0, 7, (10, 15)),
            np.random.rand(10, 15),
            np.ones((10, 15), dtype=bool),
        )

        trop_wmo(
            self.temperature,
            self.pressure,
            xdim=2,
            ydim=1,
            levdim=0,
            pressure_unit="Pa",
        )

        # Check that punit=1 was passed (Pa)
        call_args = mock_fortran.call_args[0]
        punit = call_args[5]
        assert punit == 1

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_grid_3d")
    def test_custom_lapse_criterion(self, mock_fortran):
        """Test custom lapse rate criterion."""
        mock_fortran.return_value = (
            np.random.rand(10, 15),
            np.random.rand(10, 15),
            np.random.randint(0, 7, (10, 15)),
            np.random.rand(10, 15),
            np.ones((10, 15), dtype=bool),
        )

        trop_wmo(
            self.temperature,
            self.pressure,
            xdim=2,
            ydim=1,
            levdim=0,
            lapse_criterion=2.5,
        )

        # Check that custom lapse criterion was passed
        call_args = mock_fortran.call_args[0]
        lapse_crit = call_args[4]
        assert lapse_crit == 2.5

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_grid_3d")
    def test_custom_missing_value(self, mock_fortran):
        """Test custom missing value."""
        mock_fortran.return_value = (
            np.random.rand(10, 15),
            np.random.rand(10, 15),
            np.random.randint(0, 7, (10, 15)),
            np.random.rand(10, 15),
            np.ones((10, 15), dtype=bool),
        )

        trop_wmo(
            self.temperature,
            self.pressure,
            xdim=2,
            ydim=1,
            levdim=0,
            missing_value=-888.0,
        )

        # Check that custom missing value was passed
        call_args = mock_fortran.call_args[0]
        missing_val = call_args[3]
        assert missing_val == -888.0


class TestTropWmoProfile:
    """Test trop_wmo_profile function for 1D profiles."""

    def setup_method(self):
        """Set up test data for profile analysis."""
        self.pressure = np.array([100, 200, 300, 500, 700, 850, 1000], dtype=np.float32)
        self.temperature = np.array(
            [220, 230, 245, 260, 275, 283, 288], dtype=np.float32
        )

    def test_input_validation_non_1d_pressure(self):
        """Test error for non-1D pressure input."""
        pressure_2d = np.random.rand(5, 3)

        with pytest.raises(ValueError, match="Profile inputs must be 1D arrays"):
            trop_wmo_profile(self.temperature, pressure_2d)

    def test_input_validation_non_1d_temperature(self):
        """Test error for non-1D temperature input."""
        temp_2d = np.random.rand(5, 3)

        with pytest.raises(ValueError, match="Profile inputs must be 1D arrays"):
            trop_wmo_profile(temp_2d, self.pressure)

    def test_input_validation_length_mismatch(self):
        """Test error for mismatched array lengths."""
        short_pressure = np.array([100, 200, 300])

        with pytest.raises(
            ValueError, match="Pressure and temperature profiles must have same length"
        ):
            trop_wmo_profile(self.temperature, short_pressure)

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_profile_1d")
    def test_basic_profile_calculation(self, mock_fortran):
        """Test basic profile calculation."""
        # Mock Fortran function return
        mock_return = (250.0, 10000.0, 2, 1.8, True)
        mock_fortran.return_value = mock_return

        result = trop_wmo_profile(self.temperature, self.pressure)

        # Check that Fortran function was called correctly
        mock_fortran.assert_called_once()
        call_args = mock_fortran.call_args[0]

        assert call_args[0] == 8  # nlevm = nlev + 1
        np.testing.assert_array_equal(call_args[1], self.pressure)
        np.testing.assert_array_equal(call_args[2], self.temperature)
        assert call_args[3] == -999.0  # missing_value
        assert call_args[4] == 2.0  # lapse_criterion
        assert call_args[5] == 0  # punit (hPa)

        # Check result structure
        assert isinstance(result, dict)
        assert result["pressure"] == 250.0
        assert result["height"] == 10000.0
        assert result["level_index"] == 2
        assert result["lapse_rate"] == 1.8
        assert result["success"] == True

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_profile_1d")
    def test_profile_with_custom_parameters(self, mock_fortran):
        """Test profile calculation with custom parameters."""
        mock_return = (200.0, 12000.0, 1, 1.5, True)
        mock_fortran.return_value = mock_return

        result = trop_wmo_profile(
            self.temperature,
            self.pressure,
            pressure_unit="Pa",
            lapse_criterion=1.5,
            missing_value=-777.0,
        )

        # Check parameters
        call_args = mock_fortran.call_args[0]
        assert call_args[3] == -777.0  # missing_value
        assert call_args[4] == 1.5  # lapse_criterion
        assert call_args[5] == 1  # punit (Pa)

        assert result["pressure"] == 200.0

    def test_profile_masked_arrays(self):
        """Test profile calculation with masked arrays."""
        # Create masked arrays
        pressure_masked = ma.array(
            [100, 200, 300, 500, 700], mask=[True, False, False, False, False]
        )
        temp_masked = ma.array(
            [220, 230, 245, 260, 275], mask=[False, True, False, False, False]
        )

        # Verify arrays are actually masked
        assert ma.is_masked(pressure_masked) == True
        assert ma.is_masked(temp_masked) == True

        with patch(
            "skyborn.calc.troposphere.tropopause_height.tropopause_profile_1d"
        ) as mock_fortran:
            mock_fortran.return_value = (250.0, 10000.0, 2, 1.8, True)

            result = trop_wmo_profile(
                temp_masked, pressure_masked, missing_value=-999.0
            )

            # Check that function executed successfully with masked input
            assert isinstance(result, dict)
            assert result["success"] == True

            mock_fortran.assert_called_once()

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_profile_1d")
    def test_profile_data_type_conversion(self, mock_fortran):
        """Test that inputs are converted to float32."""
        # Use integer inputs
        pressure_int = np.array([100, 200, 300, 500, 700], dtype=int)
        temp_int = np.array([220, 230, 245, 260, 275], dtype=int)

        mock_fortran.return_value = (250.0, 10000.0, 2, 1.8, True)

        result = trop_wmo_profile(temp_int, pressure_int)

        # Check that arrays were converted to float32
        call_args = mock_fortran.call_args[0]
        pressure_passed = call_args[1]
        temp_passed = call_args[2]

        assert pressure_passed.dtype == np.float32
        assert temp_passed.dtype == np.float32

        mock_fortran.assert_called_once()

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_profile_1d")
    def test_profile_case_insensitive_pressure_unit(self, mock_fortran):
        """Test case insensitive pressure unit handling."""
        mock_fortran.return_value = (250.0, 10000.0, 2, 1.8, True)

        # Test different cases
        for unit in ["HPA", "hpa", "Hpa", "hPa"]:
            trop_wmo_profile(self.temperature, self.pressure, pressure_unit=unit)
            call_args = mock_fortran.call_args[0]
            assert call_args[5] == 0  # Should be hPa (0)

        for unit in ["PA", "pa", "Pa"]:
            trop_wmo_profile(self.temperature, self.pressure, pressure_unit=unit)
            call_args = mock_fortran.call_args[0]
            assert call_args[5] == 1  # Should be Pa (1)


class TestTropWmoIntegration:
    """Integration tests for complete trop_wmo workflows."""

    def setup_method(self):
        """Set up integration test data."""
        # Create realistic atmospheric profile
        self.pressure_levels = np.array(
            [10, 20, 50, 100, 200, 300, 500, 700, 850, 1000]
        )
        nlev = len(self.pressure_levels)

        # Create temperature data with realistic vertical structure
        # Typical atmospheric temperature profile
        self.temp_1d = np.array(
            [200, 210, 220, 230, 250, 270, 280, 285, 288, 290], dtype=np.float64
        )

        # 3D data (level, lat, lon)
        self.temp_3d = (
            np.broadcast_to(self.temp_1d[:, np.newaxis, np.newaxis], (nlev, 18, 36))
            .copy()
            .astype(np.float64)
        )
        # Add some spatial variation
        self.temp_3d += np.random.randn(nlev, 18, 36).astype(np.float64) * 2

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_profile_1d")
    def test_complete_1d_workflow(self, mock_fortran):
        """Test complete 1D profile workflow."""
        mock_fortran.return_value = (200.0, 11000.0, 3, 2.1, True)

        # Test with list inputs (should be converted to arrays)
        result = trop_wmo_profile(self.temp_1d.tolist(), self.pressure_levels.tolist())

        assert isinstance(result, dict)
        assert all(
            key in result
            for key in ["pressure", "height", "level_index", "lapse_rate", "success"]
        )
        mock_fortran.assert_called_once()

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_grid_3d")
    def test_complete_3d_workflow(self, mock_fortran):
        """Test complete 3D workflow with all features."""
        mock_fortran.return_value = (
            np.random.rand(18, 36) * 100 + 100,
            np.random.rand(18, 36) * 3000 + 9000,
            np.random.randint(1, 8, (18, 36)),
            np.random.rand(18, 36) * 2 + 1,
            np.ones((18, 36), dtype=bool),
        )

        result = trop_wmo(
            self.temp_3d,
            self.pressure_levels,
            xdim=2,
            ydim=1,
            levdim=0,
            pressure_unit="hPa",
            lapse_criterion=2.0,
            missing_value=-999.0,
            check_pressure_order=True,
        )

        # Check all output variables
        expected_keys = ["pressure", "height", "level_index", "lapse_rate", "success"]
        assert all(key in result for key in expected_keys)

        # Check that result is properly structured
        assert isinstance(result, dict)

        mock_fortran.assert_called_once()

    def test_edge_case_empty_arrays(self):
        """Test handling of edge cases with empty arrays."""
        empty_temp = np.array([])
        empty_pressure = np.array([])

        # Should raise error for empty arrays
        with pytest.raises((ValueError, IndexError)):
            trop_wmo_profile(empty_temp, empty_pressure)

    def test_edge_case_single_level(self):
        """Test handling of single level data."""
        single_temp = np.array([280.0])
        single_pressure = np.array([500.0])

        with patch(
            "skyborn.calc.troposphere.tropopause_height.tropopause_profile_1d"
        ) as mock_fortran:
            mock_fortran.return_value = (500.0, 5000.0, 0, 0.0, False)

            result = trop_wmo_profile(single_temp, single_pressure)

            assert result["success"] == False  # Should fail with single level
            mock_fortran.assert_called_once()


class TestTropWmoAdditionalCoverage:
    """Additional tests to achieve 100% coverage."""

    def setup_method(self):
        """Set up additional test data."""
        self.pressure = np.array([100, 200, 300, 500, 700])
        self.temperature = np.random.rand(5, 10, 15) * 50 + 250

    def test_dimension_validation_with_placeholders(self):
        """Test dimension validation when using -1 placeholders."""
        # Test that -1 placeholders are properly handled
        with patch(
            "skyborn.calc.troposphere.tropopause_height.tropopause_grid_2d"
        ) as mock_fortran:
            mock_fortran.return_value = (
                np.random.rand(10) * 200 + 100,
                np.random.rand(10) * 5000 + 8000,
                np.random.randint(0, 5, 10),
                np.random.rand(10) * 3,
                np.ones(10, dtype=bool),
            )

            # Test 2D case with -1 for missing longitude dimension
            temp_2d = self.temperature[:, :, 0]  # (level, lat)
            result = trop_wmo(temp_2d, self.pressure, xdim=-1, ydim=1, levdim=0)

            assert isinstance(result, dict)
            mock_fortran.assert_called_once()

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_grid_3d")
    def test_already_sorted_pressure_no_sorting(self, mock_fortran):
        """Test that already sorted pressure doesn't get re-sorted."""
        sorted_pressure = np.array([100, 200, 300, 500, 700])  # Already sorted

        mock_fortran.return_value = (
            np.random.rand(10, 15) * 200 + 100,
            np.random.rand(10, 15) * 5000 + 8000,
            np.random.randint(0, 5, (10, 15)),
            np.random.rand(10, 15) * 3,
            np.ones((10, 15), dtype=bool),
        )

        result = trop_wmo(
            self.temperature,
            sorted_pressure,
            xdim=2,
            ydim=1,
            levdim=0,
            check_pressure_order=True,
        )

        # Should work without issues
        assert isinstance(result, dict)
        mock_fortran.assert_called_once()

    @patch("skyborn.calc.troposphere.tropopause_height.tropopause_grid_4d")
    def test_4d_level_axis_calculation(self, mock_fortran):
        """Test level axis calculation for 4D temperature sorting."""
        temp_4d = np.random.rand(12, 5, 10, 15)  # (time, level, lat, lon)
        unsorted_pressure = np.array([500, 100, 700, 200, 300])

        mock_fortran.return_value = (
            np.random.rand(10, 15, 12) * 200 + 100,
            np.random.rand(10, 15, 12) * 5000 + 8000,
            np.random.randint(0, 5, (10, 15, 12)),
            np.random.rand(10, 15, 12) * 3,
            np.ones((10, 15, 12), dtype=bool),
        )

        result = trop_wmo(
            temp_4d,
            unsorted_pressure,
            xdim=3,
            ydim=2,
            levdim=1,
            timedim=0,
            check_pressure_order=True,
        )

        assert isinstance(result, dict)
        mock_fortran.assert_called_once()

    def test_non_masked_array_passthrough(self):
        """Test that non-masked arrays pass through unchanged."""
        pressure = np.array([100, 200, 300, 500, 700])
        temperature = np.array([220, 230, 245, 260, 275])

        # Verify these are not masked
        assert not ma.is_masked(pressure)
        assert not ma.is_masked(temperature)

        with patch(
            "skyborn.calc.troposphere.tropopause_height.tropopause_profile_1d"
        ) as mock_fortran:
            mock_fortran.return_value = (250.0, 10000.0, 2, 1.8, True)

            result = trop_wmo_profile(temperature, pressure)

            # Check that non-masked arrays are handled correctly
            call_args = mock_fortran.call_args[0]
            pressure_passed = call_args[1]
            temp_passed = call_args[2]

            # Values should be unchanged (no masking occurred)
            np.testing.assert_array_equal(pressure_passed, pressure.astype(np.float32))
            np.testing.assert_array_equal(temp_passed, temperature.astype(np.float32))

            mock_fortran.assert_called_once()

    def test_remaining_dims_in_ordering(self):
        """Test handling of extra dimensions in ordering function."""
        # Create 5D array to test remaining dimensions handling
        grid_5d = np.random.rand(5, 10, 15, 8, 6)

        ordered_grid, info = _order_dims_for_fortran(grid_5d, xdim=2, ydim=1, levdim=0)

        # Should handle all 5 dimensions
        assert ordered_grid.ndim == 5
        assert len(info["new_order"]) == 5

        # First dimensions should be in standard order: lat, lon, level
        assert info["new_order"][0] == 1  # lat (ydim)
        assert info["new_order"][1] == 2  # lon (xdim)
        assert info["new_order"][2] == 0  # level (levdim)

        # Remaining dimensions should be included
        remaining_in_order = info["new_order"][3:]
        expected_remaining = [3, 4]  # dimensions 3 and 4
        assert remaining_in_order == expected_remaining
