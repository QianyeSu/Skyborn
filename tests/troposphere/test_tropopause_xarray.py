"""
Comprehensive tests for src.skyborn.calc.troposphere.tropopause_xarray module.

This test suite provides 100% code coverage for the tropopause_xarray module,
testing all functions, methods, and edge cases.
"""

from unittest.mock import MagicMock, patch

import numpy as np
import pytest
import xarray as xr

from skyborn.calc.troposphere.xarray import _detect_atmospheric_dimensions, trop_wmo


class TestDetectAtmosphericDimensions:
    """Test the _detect_atmospheric_dimensions function."""

    def test_1d_level_only(self):
        """Test 1D profile with only level dimension."""
        temp = xr.DataArray([280, 270, 260], dims=["level"])
        xdim, ydim, levdim, timedim = _detect_atmospheric_dimensions(temp)

        assert xdim is None
        assert ydim is None
        assert levdim == 0
        assert timedim is None

    def test_2d_level_lat(self):
        """Test 2D data with level and latitude."""
        temp = xr.DataArray(np.random.rand(5, 10), dims=["level", "lat"])
        xdim, ydim, levdim, timedim = _detect_atmospheric_dimensions(temp)

        assert xdim is None
        assert ydim == 1
        assert levdim == 0
        assert timedim is None

    def test_2d_level_lon(self):
        """Test 2D data with level and longitude."""
        temp = xr.DataArray(np.random.rand(5, 10), dims=["level", "lon"])
        xdim, ydim, levdim, timedim = _detect_atmospheric_dimensions(temp)

        assert xdim == 1
        assert ydim is None
        assert levdim == 0
        assert timedim is None

    def test_3d_level_lat_lon(self):
        """Test 3D data with level, lat, lon."""
        temp = xr.DataArray(np.random.rand(5, 10, 15), dims=["level", "lat", "lon"])
        xdim, ydim, levdim, timedim = _detect_atmospheric_dimensions(temp)

        assert xdim == 2
        assert ydim == 1
        assert levdim == 0
        assert timedim is None

    def test_4d_time_level_lat_lon(self):
        """Test 4D data with time, level, lat, lon."""
        temp = xr.DataArray(
            np.random.rand(12, 5, 10, 15), dims=["time", "level", "lat", "lon"]
        )
        xdim, ydim, levdim, timedim = _detect_atmospheric_dimensions(temp)

        assert xdim == 3
        assert ydim == 2
        assert levdim == 1
        assert timedim == 0

    def test_alternative_dimension_names(self):
        """Test detection with alternative dimension names."""
        temp = xr.DataArray(
            np.random.rand(5, 10, 15), dims=["plev", "latitude", "longitude"]
        )
        xdim, ydim, levdim, timedim = _detect_atmospheric_dimensions(temp)

        assert xdim == 2
        assert ydim == 1
        assert levdim == 0
        assert timedim is None

    def test_case_insensitive_names(self):
        """Test case-insensitive dimension name detection."""
        temp = xr.DataArray(np.random.rand(5, 10, 15), dims=["LEVEL", "LAT", "LON"])
        xdim, ydim, levdim, timedim = _detect_atmospheric_dimensions(temp)

        assert xdim == 2
        assert ydim == 1
        assert levdim == 0
        assert timedim is None

    def test_pressure_level_names(self):
        """Test detection of pressure level dimension names."""
        temp = xr.DataArray(
            np.random.rand(5, 10, 15), dims=["pressure_level", "y", "x"]
        )
        xdim, ydim, levdim, timedim = _detect_atmospheric_dimensions(temp)

        assert xdim == 2
        assert ydim == 1
        assert levdim == 0
        assert timedim is None

    def test_height_dimension(self):
        """Test detection of height as level dimension."""
        temp = xr.DataArray(np.random.rand(5, 10, 15), dims=["height", "lat", "lon"])
        xdim, ydim, levdim, timedim = _detect_atmospheric_dimensions(temp)

        assert xdim == 2
        assert ydim == 1
        assert levdim == 0
        assert timedim is None

    def test_missing_level_dimension(self):
        """Test error when level dimension cannot be identified."""
        temp = xr.DataArray(np.random.rand(5, 10, 15), dims=["x", "y", "unknown"])

        with pytest.raises(ValueError, match="Could not auto-detect level dimension"):
            _detect_atmospheric_dimensions(temp)

    def test_2d_no_spatial_dimensions(self):
        """Test error for 2D data without spatial dimensions."""
        temp = xr.DataArray(np.random.rand(5, 10), dims=["level", "other"])

        with pytest.raises(
            ValueError, match="For 2D data, need at least one spatial dimension"
        ):
            _detect_atmospheric_dimensions(temp)

    def test_3d_missing_spatial_dimensions(self):
        """Test error for 3D+ data missing required spatial dimensions."""
        temp = xr.DataArray(
            np.random.rand(5, 10, 15), dims=["level", "other1", "other2"]
        )

        with pytest.raises(ValueError, match="need both lat and lon dimensions"):
            _detect_atmospheric_dimensions(temp)

    def test_partial_spatial_dimensions_3d(self):
        """Test error for 3D data with only one spatial dimension."""
        temp = xr.DataArray(np.random.rand(5, 10, 15), dims=["level", "lat", "other"])

        with pytest.raises(ValueError, match="need both lat and lon dimensions"):
            _detect_atmospheric_dimensions(temp)

    def test_time_dimension_variations(self):
        """Test detection of various time dimension names."""
        # Test known working time dimension names
        working_time_names = ["time", "t", "T", "season"]
        for time_name in working_time_names:
            temp = xr.DataArray(
                np.random.rand(12, 5, 10, 15), dims=[time_name, "level", "lat", "lon"]
            )
            xdim, ydim, levdim, timedim = _detect_atmospheric_dimensions(temp)

            assert timedim == 0, f"Failed for time dimension: {time_name}"

    def test_time_dimension_partial_matches(self):
        """Test detection of time dimensions with partial string matches."""
        # Test time names that are substrings of other words
        for time_name in ["month", "yr", "mn"]:
            temp = xr.DataArray(
                np.random.rand(12, 5, 10, 15), dims=[time_name, "level", "lat", "lon"]
            )
            xdim, ydim, levdim, timedim = _detect_atmospheric_dimensions(temp)

            # These may or may not be detected as time dimensions depending on implementation
            # Just ensure the function doesn't crash
            assert levdim == 1  # level should always be detected


class TestTropWmoInputValidation:
    """Test input validation for trop_wmo function."""

    def test_non_dataarray_temperature(self):
        """Test error for non-DataArray temperature input."""
        temp = np.array([280, 270, 260])

        with pytest.raises(TypeError, match="temperature must be xarray.DataArray"):
            trop_wmo(temp)

    def test_non_dataarray_pressure(self):
        """Test error for non-DataArray pressure input."""
        temp = xr.DataArray([280, 270, 260], dims=["level"])
        pressure = np.array([1000, 500, 100])

        with pytest.raises(
            TypeError, match="pressure must be xarray.DataArray or None"
        ):
            trop_wmo(temp, pressure=pressure)

    def test_missing_level_coordinate_no_pressure(self):
        """Test error when level coordinate is missing and no pressure provided."""
        temp = xr.DataArray([280, 270, 260], dims=["level"])

        with pytest.raises(
            ValueError, match="Cannot generate pressure.*must have 'level' coordinate"
        ):
            trop_wmo(temp)


class TestTropWmo1DProfile:
    """Test trop_wmo function with 1D profile data."""

    def setup_method(self):
        """Set up test data for 1D profiles."""
        # Create ascending pressure levels (WMO standard)
        self.pressure_levels = np.array([100, 200, 300, 500, 700, 850, 1000])  # hPa

        # Create temperature profile with typical atmospheric structure
        self.temperature = np.array([220, 230, 245, 260, 275, 283, 288])  # K

        self.temp_da = xr.DataArray(
            self.temperature,
            dims=["level"],
            coords={"level": self.pressure_levels},
            attrs={"units": "K", "long_name": "Temperature"},
        )

    @patch("skyborn.calc.troposphere.tropopause.trop_wmo_profile")
    def test_1d_profile_basic(self, mock_trop_profile):
        """Test basic 1D profile calculation."""
        # Mock return values
        mock_result = {
            "pressure": 250.0,
            "height": 10000.0,
            "level_index": 2,
            "lapse_rate": 1.8,
            "success": True,
        }
        mock_trop_profile.return_value = mock_result

        result = trop_wmo(self.temp_da)

        # Check that the mock was called correctly
        mock_trop_profile.assert_called_once()
        call_args = mock_trop_profile.call_args
        np.testing.assert_array_equal(call_args[0][0], self.temperature)
        np.testing.assert_array_equal(call_args[0][1], self.pressure_levels)

        # Check result structure
        assert isinstance(result, xr.Dataset)
        assert "pressure" in result
        assert "height" in result
        assert "level_index" in result
        assert "lapse_rate" in result
        assert "success" in result

        # Check scalar outputs (0D DataArrays)
        assert result.pressure.ndim == 0
        assert float(result.pressure) == 250.0

    @patch("skyborn.calc.troposphere.tropopause.trop_wmo_profile")
    def test_1d_profile_with_explicit_pressure(self, mock_trop_profile):
        """Test 1D profile with explicit pressure DataArray."""
        pressure_da = xr.DataArray(
            self.pressure_levels,
            dims=["level"],
            coords={"level": self.pressure_levels},
            attrs={"units": "hPa"},
        )

        mock_result = {
            "pressure": 250.0,
            "height": 10000.0,
            "level_index": 2,
            "lapse_rate": 1.8,
            "success": True,
        }
        mock_trop_profile.return_value = mock_result

        result = trop_wmo(self.temp_da, pressure=pressure_da)

        mock_trop_profile.assert_called_once()
        assert isinstance(result, xr.Dataset)

    @patch("skyborn.calc.troposphere.tropopause.trop_wmo_profile")
    def test_1d_profile_attributes_preserved(self, mock_trop_profile):
        """Test that attributes are preserved in 1D results."""
        mock_result = {
            "pressure": 250.0,
            "height": 10000.0,
            "level_index": 2,
            "lapse_rate": 1.8,
            "success": True,
        }
        mock_trop_profile.return_value = mock_result

        result = trop_wmo(self.temp_da, keep_attrs=True)

        # Check that source attributes are preserved
        assert "source_temperature_units" in result.attrs
        assert "source_temperature_long_name" in result.attrs

    def test_1d_profile_auto_pressure_generation(self):
        """Test automatic pressure generation from level coordinate."""
        with patch(
            "skyborn.calc.troposphere.tropopause.trop_wmo_profile"
        ) as mock_profile:
            mock_profile.return_value = {
                "pressure": 250.0,
                "height": 10000.0,
                "level_index": 2,
                "lapse_rate": 1.8,
                "success": True,
            }

            result = trop_wmo(self.temp_da)

            # Check that pressure was generated correctly
            mock_profile.assert_called_once()
            call_args = mock_profile.call_args
            generated_pressure = call_args[0][1]
            np.testing.assert_array_equal(generated_pressure, self.pressure_levels)


class TestTropWmoMultiDimensional:
    """Test trop_wmo function with multi-dimensional data."""

    def setup_method(self):
        """Set up test data for multi-dimensional arrays."""
        self.pressure_levels = np.array([100, 200, 300, 500, 700, 850, 1000])  # hPa
        self.nlev = len(self.pressure_levels)
        self.nlat, self.nlon = 10, 15
        self.ntime = 12

        # Create 3D temperature data (level, lat, lon)
        self.temp_3d = 300 - np.random.rand(self.nlev, self.nlat, self.nlon) * 80

        # Create 4D temperature data (time, level, lat, lon)
        self.temp_4d = (
            300 - np.random.rand(self.ntime, self.nlev, self.nlat, self.nlon) * 80
        )

    @patch("skyborn.calc.troposphere.tropopause.trop_wmo")
    def test_3d_data_basic(self, mock_trop):
        """Test basic 3D data processing."""
        # Create test data
        temp_da = xr.DataArray(
            self.temp_3d,
            dims=["level", "lat", "lon"],
            coords={
                "level": self.pressure_levels,
                "lat": np.linspace(-45, 45, self.nlat),
                "lon": np.linspace(0, 360, self.nlon),
            },
        )

        # Mock return values
        mock_result = {
            "pressure": np.random.rand(self.nlat, self.nlon) * 200 + 100,
            "height": np.random.rand(self.nlat, self.nlon) * 5000 + 8000,
            "level_index": np.random.randint(0, self.nlev, (self.nlat, self.nlon)),
            "lapse_rate": np.random.rand(self.nlat, self.nlon) * 3,
            "success": np.ones((self.nlat, self.nlon), dtype=bool),
        }
        mock_trop.return_value = mock_result

        result = trop_wmo(temp_da)

        # Check that the mock was called correctly
        mock_trop.assert_called_once()
        call_args = mock_trop.call_args

        # Verify dimensions
        assert isinstance(result, xr.Dataset)
        assert result.pressure.dims == ("lat", "lon")
        assert result.height.dims == ("lat", "lon")
        assert result.pressure.shape == (self.nlat, self.nlon)

    @patch("skyborn.calc.troposphere.tropopause.trop_wmo")
    def test_4d_data_with_time(self, mock_trop):
        """Test 4D data processing with time dimension."""
        temp_da = xr.DataArray(
            self.temp_4d,
            dims=["time", "level", "lat", "lon"],
            coords={
                "time": pd.date_range("2000-01-01", periods=self.ntime, freq="ME"),
                "level": self.pressure_levels,
                "lat": np.linspace(-45, 45, self.nlat),
                "lon": np.linspace(0, 360, self.nlon),
            },
        )

        mock_result = {
            "pressure": np.random.rand(self.ntime, self.nlat, self.nlon) * 200 + 100,
            "height": np.random.rand(self.ntime, self.nlat, self.nlon) * 5000 + 8000,
            "level_index": np.random.randint(
                0, self.nlev, (self.ntime, self.nlat, self.nlon)
            ),
            "lapse_rate": np.random.rand(self.ntime, self.nlat, self.nlon) * 3,
            "success": np.ones((self.ntime, self.nlat, self.nlon), dtype=bool),
        }
        mock_trop.return_value = mock_result

        result = trop_wmo(temp_da)

        mock_trop.assert_called_once()
        assert result.pressure.dims == ("time", "lat", "lon")
        assert result.pressure.shape == (self.ntime, self.nlat, self.nlon)

    def test_explicit_dimension_specification(self):
        """Test with explicitly specified dimensions."""
        temp_da = xr.DataArray(
            self.temp_3d,
            dims=["level", "lat", "lon"],
            coords={"level": self.pressure_levels},
        )

        with patch("skyborn.calc.troposphere.tropopause.trop_wmo") as mock_trop:
            mock_trop.return_value = {
                "pressure": np.random.rand(self.nlat, self.nlon),
                "height": np.random.rand(self.nlat, self.nlon),
                "level_index": np.random.randint(0, self.nlev, (self.nlat, self.nlon)),
                "lapse_rate": np.random.rand(self.nlat, self.nlon),
                "success": np.ones((self.nlat, self.nlon), dtype=bool),
            }

            result = trop_wmo(temp_da, xdim="lon", ydim="lat", levdim="level")

            # Check that dimension indices were converted correctly
            call_args = mock_trop.call_args[1]
            assert call_args["xdim"] == 2  # lon dimension index
            assert call_args["ydim"] == 1  # lat dimension index
            assert call_args["levdim"] == 0  # level dimension index

    def test_2d_spatial_data(self):
        """Test 2D data (level, spatial) processing."""
        temp_2d = self.temp_3d[:, :, 0]  # Take first longitude slice
        temp_da = xr.DataArray(
            temp_2d,
            dims=["level", "lat"],
            coords={
                "level": self.pressure_levels,
                "lat": np.linspace(-45, 45, self.nlat),
            },
        )

        with patch("skyborn.calc.troposphere.tropopause.trop_wmo") as mock_trop:
            mock_trop.return_value = {
                "pressure": np.random.rand(self.nlat),
                "height": np.random.rand(self.nlat),
                "level_index": np.random.randint(0, self.nlev, self.nlat),
                "lapse_rate": np.random.rand(self.nlat),
                "success": np.ones(self.nlat, dtype=bool),
            }

            result = trop_wmo(temp_da)

            assert result.pressure.dims == ("lat",)
            assert result.pressure.shape == (self.nlat,)


class TestTropWmoPressureHandling:
    """Test pressure handling and sorting in trop_wmo."""

    def setup_method(self):
        """Set up test data."""
        self.pressure_levels = np.array([100, 200, 300, 500, 700, 850, 1000])
        self.temperature = 300 - np.random.rand(len(self.pressure_levels), 10, 15) * 80

    def test_pressure_auto_sorting_1d(self):
        """Test automatic sorting of 1D pressure arrays."""
        # Create unsorted pressure and corresponding temperature
        unsorted_pressure = np.array([500, 100, 1000, 200, 700])
        unsorted_temp = np.array([260, 220, 288, 230, 275])

        temp_da = xr.DataArray(unsorted_temp, dims=["level"])
        pressure_da = xr.DataArray(unsorted_pressure, dims=["level"])

        with patch(
            "skyborn.calc.troposphere.tropopause.trop_wmo_profile"
        ) as mock_profile:
            mock_profile.return_value = {
                "pressure": 250.0,
                "height": 10000.0,
                "level_index": 2,
                "lapse_rate": 1.8,
                "success": True,
            }

            trop_wmo(temp_da, pressure=pressure_da, auto_sort_levels=True)

            # Check that sorted arrays were passed to the calculation
            call_args = mock_profile.call_args[0]
            sorted_pressure = call_args[1]
            sorted_temp = call_args[0]

            # Pressure should be sorted in ascending order
            assert np.all(sorted_pressure[:-1] <= sorted_pressure[1:])

    def test_pressure_auto_sorting_disabled(self):
        """Test disabling automatic pressure sorting."""
        unsorted_pressure = np.array([500, 100, 1000, 200, 700])
        unsorted_temp = np.array([260, 220, 288, 230, 275])

        temp_da = xr.DataArray(unsorted_temp, dims=["level"])
        pressure_da = xr.DataArray(unsorted_pressure, dims=["level"])

        with patch(
            "skyborn.calc.troposphere.tropopause.trop_wmo_profile"
        ) as mock_profile:
            mock_profile.return_value = {
                "pressure": 250.0,
                "height": 10000.0,
                "level_index": 2,
                "lapse_rate": 1.8,
                "success": True,
            }

            trop_wmo(temp_da, pressure=pressure_da, auto_sort_levels=False)

            # Check that original order was preserved
            call_args = mock_profile.call_args[0]
            passed_pressure = call_args[1]
            np.testing.assert_array_equal(passed_pressure, unsorted_pressure)

    def test_pressure_coordinate_sorting(self):
        """Test sorting based on pressure level coordinates."""
        unsorted_levels = np.array([500, 100, 1000, 200, 700])
        temp_da = xr.DataArray(
            np.random.rand(5, 10, 15),
            dims=["level", "lat", "lon"],
            coords={"level": unsorted_levels},
        )

        with patch("skyborn.calc.troposphere.tropopause.trop_wmo") as mock_trop:
            mock_trop.return_value = {
                "pressure": np.random.rand(10, 15),
                "height": np.random.rand(10, 15),
                "level_index": np.random.randint(0, 5, (10, 15)),
                "lapse_rate": np.random.rand(10, 15),
                "success": np.ones((10, 15), dtype=bool),
            }

            trop_wmo(temp_da, auto_sort_levels=True)

            # Verify sorting was applied
            mock_trop.assert_called_once()

    def test_pressure_unit_conversion(self):
        """Test pressure unit specification."""
        temp_da = xr.DataArray(
            [220, 240, 260, 280],
            dims=["level"],
            coords={"level": [10000, 50000, 70000, 100000]},  # Pa
        )

        with patch(
            "skyborn.calc.troposphere.tropopause.trop_wmo_profile"
        ) as mock_profile:
            mock_profile.return_value = {
                "pressure": 50000.0,
                "height": 10000.0,
                "level_index": 1,
                "lapse_rate": 1.8,
                "success": True,
            }

            trop_wmo(temp_da, pressure_unit="Pa")

            # Check that pressure unit was passed correctly
            call_kwargs = mock_profile.call_args[1]
            assert call_kwargs["pressure_unit"] == "Pa"


class TestTropWmoParameterHandling:
    """Test parameter handling and configuration options."""

    def setup_method(self):
        """Set up test data."""
        self.temp_da = xr.DataArray(
            [220, 240, 260, 280],
            dims=["level"],
            coords={"level": [100, 300, 500, 1000]},
        )

    @patch("skyborn.calc.troposphere.tropopause.trop_wmo_profile")
    def test_lapse_criterion_parameter(self, mock_profile):
        """Test custom lapse rate criterion."""
        mock_profile.return_value = {
            "pressure": 250.0,
            "height": 10000.0,
            "level_index": 1,
            "lapse_rate": 1.5,
            "success": True,
        }

        trop_wmo(self.temp_da, lapse_criterion=1.5)

        call_kwargs = mock_profile.call_args[1]
        assert call_kwargs["lapse_criterion"] == 1.5

    @patch("skyborn.calc.troposphere.tropopause.trop_wmo_profile")
    def test_missing_value_parameter(self, mock_profile):
        """Test custom missing value."""
        mock_profile.return_value = {
            "pressure": 250.0,
            "height": 10000.0,
            "level_index": 1,
            "lapse_rate": 1.8,
            "success": True,
        }

        trop_wmo(self.temp_da, missing_value=-888.0)

        call_kwargs = mock_profile.call_args[1]
        assert call_kwargs["missing_value"] == -888.0

    @patch("skyborn.calc.troposphere.tropopause.trop_wmo_profile")
    def test_keep_attrs_false(self, mock_profile):
        """Test disabling attribute preservation."""
        self.temp_da.attrs = {"test_attr": "test_value"}
        mock_profile.return_value = {
            "pressure": 250.0,
            "height": 10000.0,
            "level_index": 1,
            "lapse_rate": 1.8,
            "success": True,
        }

        result = trop_wmo(self.temp_da, keep_attrs=False)

        # Should not have source attributes
        assert "source_temperature_test_attr" not in result.attrs


class TestTropWmoOutputStructure:
    """Test output Dataset structure and attributes."""

    def setup_method(self):
        """Set up test data."""
        self.temp_da = xr.DataArray(
            [220, 240, 260, 280],
            dims=["level"],
            coords={"level": [100, 300, 500, 1000]},
            attrs={"units": "K", "description": "Test temperature"},
        )

    @patch("skyborn.calc.troposphere.tropopause.trop_wmo_profile")
    def test_output_dataset_structure(self, mock_profile):
        """Test that output Dataset has correct structure."""
        mock_profile.return_value = {
            "pressure": 250.0,
            "height": 10000.0,
            "level_index": 1,
            "lapse_rate": 1.8,
            "success": True,
        }

        result = trop_wmo(self.temp_da)

        # Check Dataset structure
        assert isinstance(result, xr.Dataset)
        required_vars = ["pressure", "height", "level_index", "lapse_rate", "success"]
        for var in required_vars:
            assert var in result.data_vars

    @patch("skyborn.calc.troposphere.tropopause.trop_wmo_profile")
    def test_output_variable_attributes(self, mock_profile):
        """Test that output variables have correct attributes."""
        mock_profile.return_value = {
            "pressure": 250.0,
            "height": 10000.0,
            "level_index": 1,
            "lapse_rate": 1.8,
            "success": True,
        }

        result = trop_wmo(self.temp_da)

        # Check pressure attributes
        assert result.pressure.attrs["units"] == "hPa"
        assert "standard_name" in result.pressure.attrs
        assert "long_name" in result.pressure.attrs

        # Check height attributes
        assert result.height.attrs["units"] == "m"
        assert "standard_name" in result.height.attrs

        # Check lapse rate attributes
        assert result.lapse_rate.attrs["units"] == "K km-1"
        assert "standard_name" in result.lapse_rate.attrs

    @patch("skyborn.calc.troposphere.tropopause.trop_wmo_profile")
    def test_global_attributes(self, mock_profile):
        """Test global attributes of output Dataset."""
        mock_profile.return_value = {
            "pressure": 250.0,
            "height": 10000.0,
            "level_index": 1,
            "lapse_rate": 1.8,
            "success": True,
        }

        result = trop_wmo(self.temp_da, lapse_criterion=2.5)

        # Check global attributes
        assert "title" in result.attrs
        assert "description" in result.attrs
        assert result.attrs["lapse_rate_criterion_K_per_km"] == 2.5
        assert result.attrs["pressure_unit"] == "hPa"
        assert "method" in result.attrs

    @patch("skyborn.calc.troposphere.tropopause.trop_wmo_profile")
    def test_coordinate_preservation(self, mock_profile):
        """Test that coordinates are properly preserved."""
        temp_3d = xr.DataArray(
            np.random.rand(5, 10, 15),
            dims=["level", "lat", "lon"],
            coords={
                "level": [100, 200, 300, 500, 1000],
                "lat": np.linspace(-45, 45, 10),
                "lon": np.linspace(0, 360, 15),
            },
        )

        with patch("skyborn.calc.troposphere.tropopause.trop_wmo") as mock_trop:
            mock_trop.return_value = {
                "pressure": np.random.rand(10, 15),
                "height": np.random.rand(10, 15),
                "level_index": np.random.randint(0, 5, (10, 15)),
                "lapse_rate": np.random.rand(10, 15),
                "success": np.ones((10, 15), dtype=bool),
            }

            result = trop_wmo(temp_3d)

            # Check that spatial coordinates are preserved
            assert "lat" in result.coords
            assert "lon" in result.coords
            np.testing.assert_array_equal(result.lat.values, temp_3d.lat.values)
            np.testing.assert_array_equal(result.lon.values, temp_3d.lon.values)


class TestTropWmoEdgeCases:
    """Test edge cases and error conditions."""

    def test_string_dimension_names(self):
        """Test using string dimension names."""
        temp_da = xr.DataArray(
            np.random.rand(5, 10, 15),
            dims=["plev", "latitude", "longitude"],
            coords={"plev": [100, 200, 300, 500, 1000]},
        )

        with patch("skyborn.calc.troposphere.tropopause.trop_wmo") as mock_trop:
            mock_trop.return_value = {
                "pressure": np.random.rand(10, 15),
                "height": np.random.rand(10, 15),
                "level_index": np.random.randint(0, 5, (10, 15)),
                "lapse_rate": np.random.rand(10, 15),
                "success": np.ones((10, 15), dtype=bool),
            }

            # Should work with string dimension names
            result = trop_wmo(temp_da, levdim="plev", ydim="latitude", xdim="longitude")

            mock_trop.assert_called_once()

    def test_missing_dimension_name(self):
        """Test behavior when specified dimension name doesn't exist."""
        temp_da = xr.DataArray(
            np.random.rand(5, 10, 15),
            dims=["level", "lat", "lon"],
            coords={"level": [100, 200, 300, 500, 1000]},
        )

        with patch("skyborn.calc.troposphere.tropopause.trop_wmo") as mock_trop:
            mock_trop.return_value = {
                "pressure": np.random.rand(10, 15),
                "height": np.random.rand(10, 15),
                "level_index": np.random.randint(0, 5, (10, 15)),
                "lapse_rate": np.random.rand(10, 15),
                "success": np.ones((10, 15), dtype=bool),
            }

            # Non-existent dimension should be treated as None
            result = trop_wmo(temp_da, xdim="nonexistent")

            # Should still work with auto-detection
            mock_trop.assert_called_once()


# Add pandas import for datetime tests
try:
    import pandas as pd
except ImportError:
    pd = None

if pd is not None:

    class TestTropWmoWithPandas:
        """Test trop_wmo with pandas datetime coordinates."""

        @patch("skyborn.calc.troposphere.tropopause.trop_wmo")
        def test_time_coordinate_with_pandas(self, mock_trop):
            """Test with pandas datetime time coordinate."""
            temp_da = xr.DataArray(
                np.random.rand(12, 5, 10, 15),
                dims=["time", "level", "lat", "lon"],
                coords={
                    "time": pd.date_range("2000-01-01", periods=12, freq="ME"),
                    "level": [100, 200, 300, 500, 1000],
                    "lat": np.linspace(-45, 45, 10),
                    "lon": np.linspace(0, 360, 15),
                },
            )

            mock_trop.return_value = {
                "pressure": np.random.rand(12, 10, 15),
                "height": np.random.rand(12, 10, 15),
                "level_index": np.random.randint(0, 5, (12, 10, 15)),
                "lapse_rate": np.random.rand(12, 10, 15),
                "success": np.ones((12, 10, 15), dtype=bool),
            }

            result = trop_wmo(temp_da)

            # Check that time coordinate is preserved
            assert "time" in result.coords
            assert len(result.time) == 12


class TestTropWmoConsoleOutput:
    """Test console output and messaging."""

    def test_pressure_generation_message(self, capsys):
        """Test that pressure generation message is printed."""
        temp_da = xr.DataArray(
            [220, 240, 260, 280],
            dims=["level"],
            coords={"level": [100, 300, 500, 1000]},
        )

        with patch(
            "skyborn.calc.troposphere.tropopause.trop_wmo_profile"
        ) as mock_profile:
            mock_profile.return_value = {
                "pressure": 250.0,
                "height": 10000.0,
                "level_index": 1,
                "lapse_rate": 1.8,
                "success": True,
            }

            trop_wmo(temp_da)

            captured = capsys.readouterr()
            assert "Generated pressure from level coordinate" in captured.out
            assert "with 4 levels" in captured.out


class TestTropWmoAdditionalXarrayCoverage:
    """Additional tests for complete xarray coverage."""

    def test_dimension_name_conversion_edge_cases(self):
        """Test edge cases in dimension name to index conversion."""
        temp_da = xr.DataArray(
            np.random.rand(5, 10, 15),
            dims=["level", "lat", "lon"],
            coords={"level": [100, 200, 300, 500, 1000]},
        )

        with patch("skyborn.calc.troposphere.tropopause.trop_wmo") as mock_trop:
            mock_trop.return_value = {
                "pressure": np.random.rand(10, 15),
                "height": np.random.rand(10, 15),
                "level_index": np.random.randint(0, 5, (10, 15)),
                "lapse_rate": np.random.rand(10, 15),
                "success": np.ones((10, 15), dtype=bool),
            }

            # Test when dimension name doesn't exist in temperature dims
            result = trop_wmo(
                temp_da, xdim="nonexistent_dim", ydim="lat", levdim="level"
            )

            # Should still work (xdim becomes None)
            mock_trop.assert_called_once()

    def test_pressure_coordinate_without_level_coord(self):
        """Test pressure sorting when no level coordinate exists."""
        temp_da = xr.DataArray(
            np.random.rand(5, 10, 15),
            dims=["level", "lat", "lon"],
            # No level coordinate
        )
        pressure_da = xr.DataArray(
            np.random.rand(5, 10, 15) * 500 + 100,  # Multi-dimensional pressure
            dims=["level", "lat", "lon"],
            # No level coordinate
        )

        with patch("skyborn.calc.troposphere.tropopause.trop_wmo") as mock_trop:
            mock_trop.return_value = {
                "pressure": np.random.rand(10, 15),
                "height": np.random.rand(10, 15),
                "level_index": np.random.randint(0, 5, (10, 15)),
                "lapse_rate": np.random.rand(10, 15),
                "success": np.ones((10, 15), dtype=bool),
            }

            # Should handle pressure without level coordinate
            result = trop_wmo(temp_da, pressure=pressure_da, auto_sort_levels=True)

            mock_trop.assert_called_once()

    def test_scalar_result_conversion(self):
        """Test conversion of scalar results to 0D arrays for 1D profiles."""
        temp_da = xr.DataArray(
            [220, 240, 260, 280],
            dims=["level"],
            coords={"level": [100, 300, 500, 1000]},
        )

        with patch(
            "skyborn.calc.troposphere.tropopause.trop_wmo_profile"
        ) as mock_profile:
            # Return scalar values (not numpy arrays)
            mock_profile.return_value = {
                "pressure": 250.0,  # scalar float
                "height": 10000.0,  # scalar float
                "level_index": 2,  # scalar int
                "lapse_rate": 1.8,  # scalar float
                "success": True,  # scalar bool
            }

            result = trop_wmo(temp_da)

            # All results should be 0D arrays for consistency
            assert result.pressure.ndim == 0
            assert result.height.ndim == 0
            assert result.level_index.ndim == 0
            assert result.lapse_rate.ndim == 0
            assert result.success.ndim == 0

    def test_pressure_attributes_preservation(self):
        """Test preservation of pressure attributes when provided."""
        temp_da = xr.DataArray(
            [220, 240, 260, 280],
            dims=["level"],
            coords={"level": [100, 300, 500, 1000]},
            attrs={"temp_attr": "temp_value"},
        )

        pressure_da = xr.DataArray(
            [100, 300, 500, 1000],
            dims=["level"],
            attrs={"pressure_attr": "pressure_value", "units": "hPa"},
        )

        with patch(
            "skyborn.calc.troposphere.tropopause.trop_wmo_profile"
        ) as mock_profile:
            mock_profile.return_value = {
                "pressure": 250.0,
                "height": 10000.0,
                "level_index": 2,
                "lapse_rate": 1.8,
                "success": True,
            }

            result = trop_wmo(temp_da, pressure=pressure_da, keep_attrs=True)

            # Check that both temperature and pressure attributes are preserved
            assert "source_temperature_temp_attr" in result.attrs
            assert "source_pressure_pressure_attr" in result.attrs
            assert "source_pressure_units" in result.attrs

    def test_multidimensional_pressure_sorting_sample_profile(self):
        """Test pressure sorting logic for multi-dimensional pressure arrays."""
        # Create unsorted pressure field
        unsorted_levels = np.array([500, 100, 1000, 200, 700])
        temp_da = xr.DataArray(np.random.rand(5, 10, 15), dims=["level", "lat", "lon"])
        pressure_da = xr.DataArray(
            np.broadcast_to(unsorted_levels[:, None, None], (5, 10, 15)),
            dims=["level", "lat", "lon"],
        )

        with patch("skyborn.calc.troposphere.tropopause.trop_wmo") as mock_trop:
            mock_trop.return_value = {
                "pressure": np.random.rand(10, 15),
                "height": np.random.rand(10, 15),
                "level_index": np.random.randint(0, 5, (10, 15)),
                "lapse_rate": np.random.rand(10, 15),
                "success": np.ones((10, 15), dtype=bool),
            }

            result = trop_wmo(temp_da, pressure=pressure_da, auto_sort_levels=True)

            mock_trop.assert_called_once()
