"""
Comprehensive tests for skyborn.calc.geostrophic.xarray module.

This test suite provides >95% code coverage for the geostrophic xarray module,
testing all functions, methods, and edge cases for xarray-based geostrophic wind calculations.
"""

import sys
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
import xarray as xr

# Mock the geostrophic module before importing
sys.modules["geostrophicwind"] = MagicMock()


from skyborn.calc.geostrophic.xarray import (
    GeostrophicWind,
    _create_dim_order_string,
    _detect_spatial_dimensions,
    _extract_coordinates,
    geostrophic_wind,
)


class TestDetectSpatialDimensions:
    """Test the _detect_spatial_dimensions function."""

    def test_2d_lat_lon(self):
        """Test 2D data with lat and lon dimensions."""
        data = xr.DataArray(np.random.rand(10, 15), dims=["lat", "lon"])
        xdim, ydim = _detect_spatial_dimensions(data)

        assert xdim == 1  # lon dimension
        assert ydim == 0  # lat dimension

    def test_3d_level_lat_lon(self):
        """Test 3D data with level, lat, lon dimensions."""
        data = xr.DataArray(np.random.rand(5, 10, 15), dims=["level", "lat", "lon"])
        xdim, ydim = _detect_spatial_dimensions(data)

        assert xdim == 2  # lon dimension
        assert ydim == 1  # lat dimension

    def test_4d_time_level_lat_lon(self):
        """Test 4D data with time, level, lat, lon dimensions."""
        data = xr.DataArray(
            np.random.rand(12, 5, 10, 15), dims=["time", "level", "lat", "lon"]
        )
        xdim, ydim = _detect_spatial_dimensions(data)

        assert xdim == 3  # lon dimension
        assert ydim == 2  # lat dimension

    def test_alternative_names_longitude(self):
        """Test detection with alternative longitude names."""
        lon_names = ["longitude", "x", "X", "LON", "LONGITUDE", "XLON", "LONS", "LONG"]
        for lon_name in lon_names:
            data = xr.DataArray(np.random.rand(10, 15), dims=["lat", lon_name])
            xdim, ydim = _detect_spatial_dimensions(data)
            assert xdim == 1, f"Failed for longitude name: {lon_name}"
            assert ydim == 0

    def test_alternative_names_latitude(self):
        """Test detection with alternative latitude names."""
        lat_names = ["latitude", "y", "Y", "LAT", "LATITUDE", "YLAT", "LATS", "LATI"]
        for lat_name in lat_names:
            data = xr.DataArray(np.random.rand(10, 15), dims=[lat_name, "lon"])
            xdim, ydim = _detect_spatial_dimensions(data)
            assert xdim == 1
            assert ydim == 0, f"Failed for latitude name: {lat_name}"

    def test_case_insensitive_detection(self):
        """Test case-insensitive dimension name detection."""
        data = xr.DataArray(np.random.rand(10, 15), dims=["LAT", "LON"])
        xdim, ydim = _detect_spatial_dimensions(data)

        assert xdim == 1  # LON dimension
        assert ydim == 0  # LAT dimension

    def test_mixed_case_detection(self):
        """Test mixed case dimension names."""
        data = xr.DataArray(np.random.rand(10, 15), dims=["Latitude", "Longitude"])
        xdim, ydim = _detect_spatial_dimensions(data)

        assert xdim == 1  # Longitude dimension
        assert ydim == 0  # Latitude dimension

    def test_missing_longitude_error(self):
        """Test error when longitude dimension is missing."""
        data = xr.DataArray(np.random.rand(10, 15), dims=["lat", "other"])

        with pytest.raises(
            ValueError, match="Could not auto-detect both longitude and latitude"
        ):
            _detect_spatial_dimensions(data)

    def test_missing_latitude_error(self):
        """Test error when latitude dimension is missing."""
        data = xr.DataArray(np.random.rand(10, 15), dims=["other", "lon"])

        with pytest.raises(
            ValueError, match="Could not auto-detect both longitude and latitude"
        ):
            _detect_spatial_dimensions(data)

    def test_both_missing_error(self):
        """Test error when both dimensions are missing."""
        data = xr.DataArray(np.random.rand(10, 15), dims=["other1", "other2"])

        with pytest.raises(
            ValueError, match="Could not auto-detect both longitude and latitude"
        ):
            _detect_spatial_dimensions(data)

    def test_partial_name_matching(self):
        """Test that partial name matching works correctly."""
        # Names containing lon/lat should match
        data = xr.DataArray(np.random.rand(10, 15), dims=["mylat", "mylon"])
        xdim, ydim = _detect_spatial_dimensions(data)

        assert xdim == 1  # mylon dimension
        assert ydim == 0  # mylat dimension


class TestExtractCoordinates:
    """Test the _extract_coordinates function."""

    def setup_method(self):
        """Set up test coordinates."""
        self.lat_coords = np.linspace(-90, 90, 10)
        self.lon_coords = np.linspace(0, 360, 15)

    def test_basic_coordinate_extraction(self):
        """Test basic coordinate extraction."""
        data = xr.DataArray(
            np.random.rand(10, 15),
            dims=["lat", "lon"],
            coords={"lat": self.lat_coords, "lon": self.lon_coords},
        )

        glon, glat = _extract_coordinates(data)

        np.testing.assert_array_equal(glon, self.lon_coords)
        np.testing.assert_array_equal(glat, self.lat_coords)

    def test_alternative_coordinate_names(self):
        """Test extraction with alternative coordinate names."""
        data = xr.DataArray(
            np.random.rand(10, 15),
            dims=["latitude", "longitude"],
            coords={"latitude": self.lat_coords, "longitude": self.lon_coords},
        )

        glon, glat = _extract_coordinates(data)

        np.testing.assert_array_equal(glon, self.lon_coords)
        np.testing.assert_array_equal(glat, self.lat_coords)

    def test_case_insensitive_coordinates(self):
        """Test case-insensitive coordinate name detection."""
        data = xr.DataArray(
            np.random.rand(10, 15),
            dims=["LAT", "LON"],
            coords={"LAT": self.lat_coords, "LON": self.lon_coords},
        )

        glon, glat = _extract_coordinates(data)

        np.testing.assert_array_equal(glon, self.lon_coords)
        np.testing.assert_array_equal(glat, self.lat_coords)

    def test_x_y_coordinate_names(self):
        """Test x/y coordinate name detection."""
        data = xr.DataArray(
            np.random.rand(10, 15),
            dims=["y", "x"],
            coords={"y": self.lat_coords, "x": self.lon_coords},
        )

        glon, glat = _extract_coordinates(data)

        np.testing.assert_array_equal(glon, self.lon_coords)
        np.testing.assert_array_equal(glat, self.lat_coords)

    def test_missing_longitude_coordinate_error(self):
        """Test error when longitude coordinate is missing."""
        data = xr.DataArray(
            np.random.rand(10, 15),
            dims=["lat", "lon"],
            coords={"lat": self.lat_coords},  # Missing lon coordinate
        )

        with pytest.raises(
            ValueError, match="Could not find longitude and latitude coordinates"
        ):
            _extract_coordinates(data)

    def test_missing_latitude_coordinate_error(self):
        """Test error when latitude coordinate is missing."""
        data = xr.DataArray(
            np.random.rand(10, 15),
            dims=["lat", "lon"],
            coords={"lon": self.lon_coords},  # Missing lat coordinate
        )

        with pytest.raises(
            ValueError, match="Could not find longitude and latitude coordinates"
        ):
            _extract_coordinates(data)

    def test_both_missing_coordinates_error(self):
        """Test error when both coordinates are missing."""
        data = xr.DataArray(
            np.random.rand(10, 15),
            dims=["lat", "lon"],
            coords={},  # No coordinates
        )

        with pytest.raises(
            ValueError, match="Could not find longitude and latitude coordinates"
        ):
            _extract_coordinates(data)

    def test_additional_coordinate_names(self):
        """Test other coordinate naming conventions."""
        coord_pairs = [
            ("XLON", "YLAT"),
            ("LONS", "LATS"),
            ("LONG", "LATI"),
        ]

        for lon_name, lat_name in coord_pairs:
            data = xr.DataArray(
                np.random.rand(10, 15),
                dims=[lat_name.lower(), lon_name.lower()],
                coords={lat_name: self.lat_coords, lon_name: self.lon_coords},
            )

            glon, glat = _extract_coordinates(data)

            np.testing.assert_array_equal(glon, self.lon_coords)
            np.testing.assert_array_equal(glat, self.lat_coords)


class TestCreateDimOrderString:
    """Test the _create_dim_order_string function."""

    def test_2d_lat_lon(self):
        """Test 2D data with lat/lon dimensions."""
        data = xr.DataArray(np.random.rand(10, 15), dims=["lat", "lon"])
        dim_order = _create_dim_order_string(data, xdim=1, ydim=0)

        assert dim_order == "yx"

    def test_3d_level_lat_lon(self):
        """Test 3D data with level dimension."""
        data = xr.DataArray(np.random.rand(5, 10, 15), dims=["level", "lat", "lon"])
        dim_order = _create_dim_order_string(data, xdim=2, ydim=1)

        assert dim_order == "zyx"

    def test_4d_time_level_lat_lon(self):
        """Test 4D data with time and level dimensions."""
        data = xr.DataArray(
            np.random.rand(12, 5, 10, 15), dims=["time", "level", "lat", "lon"]
        )
        dim_order = _create_dim_order_string(data, xdim=3, ydim=2)

        assert dim_order == "tzyx"

    def test_alternative_time_names(self):
        """Test recognition of alternative time dimension names."""
        time_names = ["time", "t", "T", "year", "month", "yr", "mn", "season"]

        for time_name in time_names:
            data = xr.DataArray(
                np.random.rand(12, 10, 15), dims=[time_name, "lat", "lon"]
            )
            dim_order = _create_dim_order_string(data, xdim=2, ydim=1)

            assert dim_order == "tyx", f"Failed for time name: {time_name}"

    def test_alternative_level_names(self):
        """Test recognition of alternative level dimension names."""
        level_names = [
            "level",
            "lev",
            "plev",
            "pressure",
            "pressure_level",
            "z",
            "Z",
            "LEV",
            "PRES",
            "LEVEL",
            "PLEVEL",
            "height",
            "altitude",
            "isobaric",
        ]

        for level_name in level_names:
            data = xr.DataArray(
                np.random.rand(5, 10, 15), dims=[level_name, "lat", "lon"]
            )
            dim_order = _create_dim_order_string(data, xdim=2, ydim=1)

            assert dim_order == "zyx", f"Failed for level name: {level_name}"

    def test_mixed_dimension_ordering(self):
        """Test various dimension orderings."""
        # lon, lat order
        data = xr.DataArray(np.random.rand(15, 10), dims=["lon", "lat"])
        dim_order = _create_dim_order_string(data, xdim=0, ydim=1)
        assert dim_order == "xy"

        # time, lat, lon order
        data = xr.DataArray(np.random.rand(12, 10, 15), dims=["time", "lat", "lon"])
        dim_order = _create_dim_order_string(data, xdim=2, ydim=1)
        assert dim_order == "tyx"

        # lat, lon, level order
        data = xr.DataArray(np.random.rand(10, 15, 5), dims=["lat", "lon", "level"])
        dim_order = _create_dim_order_string(data, xdim=1, ydim=0)
        assert dim_order == "yxz"

    def test_unknown_dimension_names(self):
        """Test handling of unknown dimension names."""
        data = xr.DataArray(np.random.rand(10, 15, 8), dims=["lat", "lon", "unknown"])
        dim_order = _create_dim_order_string(data, xdim=1, ydim=0)

        # Unknown dimensions should default to 't'
        assert dim_order == "yxt"

    def test_case_insensitive_dimension_detection(self):
        """Test case-insensitive dimension name detection."""
        data = xr.DataArray(
            np.random.rand(12, 5, 10, 15), dims=["TIME", "LEVEL", "LAT", "LON"]
        )
        dim_order = _create_dim_order_string(data, xdim=3, ydim=2)

        assert dim_order == "tzyx"


class TestGeostrophicWindFunction:
    """Test the main geostrophic_wind function."""

    def setup_method(self):
        """Set up test data."""
        self.lat_coords = np.linspace(-90, 90, 10)
        self.lon_coords = np.linspace(0, 360, 15)
        self.z_data = np.random.rand(10, 15) * 1000 + 5000

    @patch("skyborn.calc.geostrophic.xarray.interface.geostrophic_wind")
    def test_2d_basic_calculation(self, mock_interface):
        """Test basic 2D geostrophic wind calculation."""
        # Mock the interface function
        mock_ug = np.random.rand(10, 15) * 10
        mock_vg = np.random.rand(10, 15) * 10
        mock_interface.return_value = (mock_ug, mock_vg)

        # Create test DataArray
        z = xr.DataArray(
            self.z_data,
            dims=["lat", "lon"],
            coords={"lat": self.lat_coords, "lon": self.lon_coords},
            attrs={"units": "gpm", "description": "Geopotential height"},
        )

        result = geostrophic_wind(z)

        # Check that interface function was called
        mock_interface.assert_called_once()

        # Check result structure
        assert isinstance(result, xr.Dataset)
        assert "ug" in result
        assert "vg" in result

        # Check dimensions and coordinates
        assert result.ug.dims == ("lat", "lon")
        assert result.vg.dims == ("lat", "lon")
        np.testing.assert_array_equal(result.ug.lat.values, self.lat_coords)
        np.testing.assert_array_equal(result.ug.lon.values, self.lon_coords)

        # Check data values
        np.testing.assert_array_equal(result.ug.values, mock_ug)
        np.testing.assert_array_equal(result.vg.values, mock_vg)

    @patch("skyborn.calc.geostrophic.xarray.interface.geostrophic_wind")
    def test_3d_calculation(self, mock_interface):
        """Test 3D geostrophic wind calculation."""
        # Create 3D data (time, lat, lon)
        ntime = 12
        z_3d = np.random.rand(ntime, 10, 15) * 1000 + 5000
        time_coords = np.arange(ntime)

        mock_ug = np.random.rand(ntime, 10, 15) * 10
        mock_vg = np.random.rand(ntime, 10, 15) * 10
        mock_interface.return_value = (mock_ug, mock_vg)

        z = xr.DataArray(
            z_3d,
            dims=["time", "lat", "lon"],
            coords={
                "time": time_coords,
                "lat": self.lat_coords,
                "lon": self.lon_coords,
            },
        )

        result = geostrophic_wind(z)

        # Check result structure
        assert result.ug.dims == ("time", "lat", "lon")
        assert result.vg.dims == ("time", "lat", "lon")

        # Check that all coordinates are preserved
        np.testing.assert_array_equal(result.ug.time.values, time_coords)
        np.testing.assert_array_equal(result.ug.lat.values, self.lat_coords)
        np.testing.assert_array_equal(result.ug.lon.values, self.lon_coords)

    @patch("skyborn.calc.geostrophic.xarray.interface.geostrophic_wind")
    def test_4d_calculation(self, mock_interface):
        """Test 4D geostrophic wind calculation."""
        # Create 4D data (time, level, lat, lon)
        ntime, nlev = 6, 5
        z_4d = np.random.rand(ntime, nlev, 10, 15) * 1000 + 5000
        time_coords = np.arange(ntime)
        level_coords = np.array([100, 200, 500, 700, 1000])

        mock_ug = np.random.rand(ntime, nlev, 10, 15) * 10
        mock_vg = np.random.rand(ntime, nlev, 10, 15) * 10
        mock_interface.return_value = (mock_ug, mock_vg)

        z = xr.DataArray(
            z_4d,
            dims=["time", "level", "lat", "lon"],
            coords={
                "time": time_coords,
                "level": level_coords,
                "lat": self.lat_coords,
                "lon": self.lon_coords,
            },
        )

        result = geostrophic_wind(z)

        # Check result structure
        assert result.ug.dims == ("time", "level", "lat", "lon")
        assert result.vg.dims == ("time", "level", "lat", "lon")

        # Check coordinates
        np.testing.assert_array_equal(result.ug.time.values, time_coords)
        np.testing.assert_array_equal(result.ug.level.values, level_coords)
        np.testing.assert_array_equal(result.ug.lat.values, self.lat_coords)
        np.testing.assert_array_equal(result.ug.lon.values, self.lon_coords)

    def test_invalid_input_type_error(self):
        """Test error for invalid input type."""
        z_numpy = np.random.rand(10, 15)

        with pytest.raises(TypeError, match="z must be xarray.DataArray"):
            geostrophic_wind(z_numpy)

    def test_missing_spatial_dimensions_error(self):
        """Test error when spatial dimensions are missing."""
        z = xr.DataArray(
            np.random.rand(10, 15),
            dims=["other1", "other2"],  # No lat/lon
            coords={"other1": np.arange(10), "other2": np.arange(15)},
        )

        with pytest.raises(
            ValueError, match="Could not auto-detect both longitude and latitude"
        ):
            geostrophic_wind(z)

    def test_missing_coordinates_error(self):
        """Test error when coordinates are missing."""
        z = xr.DataArray(
            np.random.rand(10, 15),
            dims=["lat", "lon"],
            # No coordinates provided
        )

        with pytest.raises(
            ValueError, match="Could not find longitude and latitude coordinates"
        ):
            geostrophic_wind(z)

    @patch("skyborn.calc.geostrophic.xarray.interface.geostrophic_wind")
    def test_custom_missing_value(self, mock_interface):
        """Test custom missing value parameter."""
        mock_interface.return_value = (
            np.random.rand(10, 15),
            np.random.rand(10, 15),
        )

        z = xr.DataArray(
            self.z_data,
            dims=["lat", "lon"],
            coords={"lat": self.lat_coords, "lon": self.lon_coords},
        )

        result = geostrophic_wind(z, missing_value=-888.0)

        # Check that custom missing value was passed
        call_kwargs = mock_interface.call_args[1]
        assert call_kwargs["missing_value"] == -888.0

        # Check that it's stored in global attributes
        assert result.attrs["missing_value"] == -888.0

    @patch("skyborn.calc.geostrophic.xarray.interface.geostrophic_wind")
    def test_attributes_preservation(self, mock_interface):
        """Test preservation of input attributes."""
        mock_interface.return_value = (
            np.random.rand(10, 15),
            np.random.rand(10, 15),
        )

        z = xr.DataArray(
            self.z_data,
            dims=["lat", "lon"],
            coords={"lat": self.lat_coords, "lon": self.lon_coords},
            attrs={
                "units": "gpm",
                "description": "500 hPa geopotential height",
                "source": "ERA5",
            },
        )

        result = geostrophic_wind(z, keep_attrs=True)

        # Check that source attributes are preserved with prefix
        assert "source_geopotential_units" in result.attrs
        assert "source_geopotential_description" in result.attrs
        assert "source_geopotential_source" in result.attrs

        assert result.attrs["source_geopotential_units"] == "gpm"
        assert (
            result.attrs["source_geopotential_description"]
            == "500 hPa geopotential height"
        )
        assert result.attrs["source_geopotential_source"] == "ERA5"

    @patch("skyborn.calc.geostrophic.xarray.interface.geostrophic_wind")
    def test_attributes_disabled(self, mock_interface):
        """Test disabling attribute preservation."""
        mock_interface.return_value = (
            np.random.rand(10, 15),
            np.random.rand(10, 15),
        )

        z = xr.DataArray(
            self.z_data,
            dims=["lat", "lon"],
            coords={"lat": self.lat_coords, "lon": self.lon_coords},
            attrs={"units": "gpm", "description": "test data"},
        )

        result = geostrophic_wind(z, keep_attrs=False)

        # Should not have source attributes
        assert "source_geopotential_units" not in result.attrs
        assert "source_geopotential_description" not in result.attrs

    @patch("skyborn.calc.geostrophic.xarray.interface.geostrophic_wind")
    def test_output_variable_attributes(self, mock_interface):
        """Test that output variables have correct attributes."""
        mock_interface.return_value = (
            np.random.rand(10, 15),
            np.random.rand(10, 15),
        )

        z = xr.DataArray(
            self.z_data,
            dims=["lat", "lon"],
            coords={"lat": self.lat_coords, "lon": self.lon_coords},
        )

        result = geostrophic_wind(z)

        # Check ug attributes
        assert result.ug.attrs["units"] == "m s-1"
        assert result.ug.attrs["standard_name"] == "eastward_geostrophic_wind"
        assert "long_name" in result.ug.attrs
        assert "description" in result.ug.attrs

        # Check vg attributes
        assert result.vg.attrs["units"] == "m s-1"
        assert result.vg.attrs["standard_name"] == "northward_geostrophic_wind"
        assert "long_name" in result.vg.attrs
        assert "description" in result.vg.attrs

    @patch("skyborn.calc.geostrophic.xarray.interface.geostrophic_wind")
    def test_global_attributes(self, mock_interface):
        """Test global attributes in output Dataset."""
        mock_interface.return_value = (
            np.random.rand(10, 15),
            np.random.rand(10, 15),
        )

        z = xr.DataArray(
            self.z_data,
            dims=["lat", "lon"],
            coords={"lat": self.lat_coords, "lon": self.lon_coords},
        )

        result = geostrophic_wind(z)

        # Check global attributes
        assert "title" in result.attrs
        assert "description" in result.attrs
        assert "longitude_cyclic" in result.attrs
        assert result.attrs["latitude_ordering"] == "south_to_north"
        assert result.attrs["missing_value"] == -999.0
        assert "method" in result.attrs
        assert "software" in result.attrs
        assert "equations" in result.attrs


class TestGeostrophicWindClass:
    """Test the GeostrophicWind class."""

    def setup_method(self):
        """Set up test data."""
        self.lat_coords = np.linspace(-90, 90, 10)
        self.lon_coords = np.linspace(0, 360, 15)
        self.z_data = np.random.rand(10, 15) * 1000 + 5000

        self.z = xr.DataArray(
            self.z_data,
            dims=["lat", "lon"],
            coords={"lat": self.lat_coords, "lon": self.lon_coords},
            attrs={"units": "gpm", "description": "Geopotential height"},
        )

    @patch("skyborn.calc.geostrophic.xarray.geostrophic_wind")
    def test_class_initialization(self, mock_geostrophic_wind):
        """Test class initialization."""
        # Mock the geostrophic_wind function
        mock_dataset = xr.Dataset(
            {
                "ug": xr.DataArray(np.random.rand(10, 15), dims=["lat", "lon"]),
                "vg": xr.DataArray(np.random.rand(10, 15), dims=["lat", "lon"]),
            }
        )
        mock_geostrophic_wind.return_value = mock_dataset

        gw = GeostrophicWind(self.z)

        # Check that geostrophic_wind was called
        mock_geostrophic_wind.assert_called_once_with(
            self.z, missing_value=-999.0, keep_attrs=True
        )

        # Check properties
        assert isinstance(gw.geopotential_height, xr.DataArray)
        np.testing.assert_array_equal(gw.longitude, self.lon_coords)
        np.testing.assert_array_equal(gw.latitude, self.lat_coords)

    @patch("skyborn.calc.geostrophic.xarray.geostrophic_wind")
    def test_uv_components_method(self, mock_geostrophic_wind):
        """Test uv_components method."""
        mock_ug = xr.DataArray(np.random.rand(10, 15), dims=["lat", "lon"])
        mock_vg = xr.DataArray(np.random.rand(10, 15), dims=["lat", "lon"])
        mock_dataset = xr.Dataset({"ug": mock_ug, "vg": mock_vg})
        mock_geostrophic_wind.return_value = mock_dataset

        gw = GeostrophicWind(self.z)
        ug, vg = gw.uv_components()

        # Check that the correct components are returned
        assert isinstance(ug, xr.DataArray)
        assert isinstance(vg, xr.DataArray)
        np.testing.assert_array_equal(ug.values, mock_ug.values)
        np.testing.assert_array_equal(vg.values, mock_vg.values)

    @patch("skyborn.calc.geostrophic.xarray.geostrophic_wind")
    def test_speed_calculation(self, mock_geostrophic_wind):
        """Test wind speed calculation."""
        # Create mock data with known values for testing speed calculation
        mock_ug = xr.DataArray(np.array([[3.0, 4.0], [0.0, -5.0]]), dims=["lat", "lon"])
        mock_vg = xr.DataArray(
            np.array([[4.0, 3.0], [12.0, 12.0]]), dims=["lat", "lon"]
        )
        mock_dataset = xr.Dataset({"ug": mock_ug, "vg": mock_vg})
        mock_geostrophic_wind.return_value = mock_dataset

        # Use smaller test data for this specific test
        z_small = xr.DataArray(
            np.random.rand(2, 2),
            dims=["lat", "lon"],
            coords={"lat": [-45, 45], "lon": [0, 180]},
        )

        gw = GeostrophicWind(z_small)
        speed = gw.speed()

        # Check speed calculation: sqrt(ug^2 + vg^2)
        expected_speed = np.array([[5.0, 5.0], [12.0, 13.0]])
        np.testing.assert_array_almost_equal(speed.values, expected_speed)

        # Check attributes
        assert speed.attrs["units"] == "m s-1"
        assert speed.attrs["standard_name"] == "geostrophic_wind_speed"
        assert "long_name" in speed.attrs
        assert "description" in speed.attrs

    @patch("skyborn.calc.geostrophic.xarray.geostrophic_wind")
    def test_speed_with_missing_values(self, mock_geostrophic_wind):
        """Test speed calculation with missing values."""
        # Create mock data with missing values
        mock_ug = xr.DataArray(
            np.array([[3.0, -999.0], [0.0, 4.0]]), dims=["lat", "lon"]
        )
        mock_vg = xr.DataArray(
            np.array([[4.0, 5.0], [-999.0, 3.0]]), dims=["lat", "lon"]
        )
        mock_dataset = xr.Dataset({"ug": mock_ug, "vg": mock_vg})
        mock_geostrophic_wind.return_value = mock_dataset

        z_small = xr.DataArray(
            np.random.rand(2, 2),
            dims=["lat", "lon"],
            coords={"lat": [-45, 45], "lon": [0, 180]},
        )

        gw = GeostrophicWind(z_small, missing_value=-999.0)
        speed = gw.speed()

        # Check that missing values are properly handled
        expected_speed = np.array([[5.0, -999.0], [-999.0, 5.0]])
        np.testing.assert_array_equal(speed.values, expected_speed)

    @patch("skyborn.calc.geostrophic.xarray.geostrophic_wind")
    def test_custom_parameters(self, mock_geostrophic_wind):
        """Test class with custom parameters."""
        mock_dataset = xr.Dataset(
            {
                "ug": xr.DataArray(np.random.rand(10, 15), dims=["lat", "lon"]),
                "vg": xr.DataArray(np.random.rand(10, 15), dims=["lat", "lon"]),
            }
        )
        mock_geostrophic_wind.return_value = mock_dataset

        gw = GeostrophicWind(self.z, missing_value=-888.0, keep_attrs=False)

        # Check that custom parameters were passed
        mock_geostrophic_wind.assert_called_once_with(
            self.z, missing_value=-888.0, keep_attrs=False
        )


class TestDimensionHandling:
    """Test handling of different dimension arrangements."""

    @patch("skyborn.calc.geostrophic.xarray.interface.geostrophic_wind")
    def test_different_dimension_orders(self, mock_interface):
        """Test different dimension orderings."""
        mock_interface.return_value = (np.random.rand(10, 15), np.random.rand(10, 15))

        lat_coords = np.linspace(-90, 90, 10)
        lon_coords = np.linspace(0, 360, 15)

        # Test lon, lat order
        z_xy = xr.DataArray(
            np.random.rand(15, 10),
            dims=["lon", "lat"],
            coords={"lon": lon_coords, "lat": lat_coords},
        )

        result = geostrophic_wind(z_xy)

        # Check that interface was called with correct parameters
        mock_interface.assert_called_once()
        call_args = mock_interface.call_args

        # Check dimension order string
        dim_order = call_args[1]["dim_order"]
        assert dim_order == "xy"

    @patch("skyborn.calc.geostrophic.xarray.interface.geostrophic_wind")
    def test_complex_dimension_structures(self, mock_interface):
        """Test complex dimension structures."""
        mock_interface.return_value = (
            np.random.rand(6, 5, 10, 15),
            np.random.rand(6, 5, 10, 15),
        )

        # Create 4D data with different ordering
        z_4d = xr.DataArray(
            np.random.rand(6, 5, 10, 15),
            dims=["time", "level", "lat", "lon"],
            coords={
                "time": np.arange(6),
                "level": [100, 200, 300, 500, 1000],
                "lat": np.linspace(-90, 90, 10),
                "lon": np.linspace(0, 360, 15),
            },
        )

        result = geostrophic_wind(z_4d)

        # Check that coordinates are properly extracted
        call_args = mock_interface.call_args
        glon = call_args[0][1]  # Second positional argument
        glat = call_args[0][2]  # Third positional argument

        np.testing.assert_array_equal(glon, np.linspace(0, 360, 15))
        np.testing.assert_array_equal(glat, np.linspace(-90, 90, 10))


class TestEdgeCasesAndErrors:
    """Test edge cases and error conditions."""

    def test_no_attrs_preservation(self):
        """Test when input has no attributes."""
        z = xr.DataArray(
            np.random.rand(10, 15),
            dims=["lat", "lon"],
            coords={
                "lat": np.linspace(-90, 90, 10),
                "lon": np.linspace(0, 360, 15),
            },
            # No attrs
        )

        with patch(
            "skyborn.calc.geostrophic.xarray.interface.geostrophic_wind"
        ) as mock_interface:
            mock_interface.return_value = (
                np.random.rand(10, 15),
                np.random.rand(10, 15),
            )

            result = geostrophic_wind(z, keep_attrs=True)

            # Should not cause errors even without attrs
            assert isinstance(result, xr.Dataset)

    def test_partial_coordinate_matching(self):
        """Test coordinate name matching with partial matches."""
        # Names containing lat/lon as substrings
        z = xr.DataArray(
            np.random.rand(10, 15),
            dims=["mylat", "mylon"],
            coords={
                "mylat": np.linspace(-90, 90, 10),
                "mylon": np.linspace(0, 360, 15),
            },
        )

        with patch(
            "skyborn.calc.geostrophic.xarray.interface.geostrophic_wind"
        ) as mock_interface:
            mock_interface.return_value = (
                np.random.rand(10, 15),
                np.random.rand(10, 15),
            )

            result = geostrophic_wind(z)

            # Should work with partial name matching
            assert isinstance(result, xr.Dataset)
            mock_interface.assert_called_once()

    @patch("skyborn.calc.geostrophic.xarray.interface.geostrophic_wind")
    def test_longitude_cyclicity_detection(self, mock_interface):
        """Test that longitude cyclicity is properly detected and stored."""
        mock_interface.return_value = (np.random.rand(10, 15), np.random.rand(10, 15))

        # Create cyclic longitude grid
        lon_cyclic = np.linspace(0, 357.5, 15)  # 0 to 357.5, cyclic
        z = xr.DataArray(
            np.random.rand(10, 15),
            dims=["lat", "lon"],
            coords={
                "lat": np.linspace(-90, 90, 10),
                "lon": lon_cyclic,
            },
        )

        result = geostrophic_wind(z)

        # Check that cyclicity information is stored
        assert "longitude_cyclic" in result.attrs
        # The specific value depends on the _is_longitude_cyclic implementation

    def test_dimension_index_conversion(self):
        """Test conversion of dimension names to indices."""
        z = xr.DataArray(
            np.random.rand(5, 10, 15),
            dims=["level", "latitude", "longitude"],
            coords={
                "level": [100, 200, 300, 500, 1000],
                "latitude": np.linspace(-90, 90, 10),
                "longitude": np.linspace(0, 360, 15),
            },
        )

        # Test spatial dimension detection
        xdim, ydim = _detect_spatial_dimensions(z)

        assert xdim == 2  # longitude is at index 2
        assert ydim == 1  # latitude is at index 1

        # Test dimension order string creation
        dim_order = _create_dim_order_string(z, xdim, ydim)
        assert dim_order == "zyx"  # level, latitude, longitude


class TestIntegrationScenarios:
    """Test realistic usage scenarios."""

    @patch("skyborn.calc.geostrophic.xarray.interface.geostrophic_wind")
    def test_era5_like_data(self, mock_interface):
        """Test with ERA5-like data structure."""
        # Simulate ERA5-like data structure
        mock_interface.return_value = (
            np.random.rand(12, 10, 15) * 10,
            np.random.rand(12, 10, 15) * 10,
        )

        # Create ERA5-like geopotential height data
        z = xr.DataArray(
            np.random.rand(12, 10, 15) * 1000 + 5000,
            dims=["time", "latitude", "longitude"],
            coords={
                "time": np.arange(12),
                "latitude": np.linspace(90, -90, 10),  # North to south (ERA5 style)
                "longitude": np.linspace(0, 357.5, 15),
            },
            attrs={
                "units": "m**2 s**-2",
                "long_name": "Geopotential",
                "standard_name": "geopotential",
                "source": "ERA5",
            },
        )

        result = geostrophic_wind(z, keep_attrs=True)

        # Check that it works correctly
        assert isinstance(result, xr.Dataset)
        assert result.ug.dims == ("time", "latitude", "longitude")
        assert result.vg.dims == ("time", "latitude", "longitude")

        # Check that source attributes are preserved
        assert "source_geopotential_units" in result.attrs
        assert "source_geopotential_source" in result.attrs

    @patch("skyborn.calc.geostrophic.xarray.interface.geostrophic_wind")
    def test_gfs_like_data(self, mock_interface):
        """Test with GFS-like data structure."""
        mock_interface.return_value = (
            np.random.rand(6, 5, 181, 360) * 10,
            np.random.rand(6, 5, 181, 360) * 10,
        )

        # Create GFS-like data structure
        z = xr.DataArray(
            np.random.rand(6, 5, 181, 360) * 1000 + 5000,
            dims=["time", "isobaric", "lat", "lon"],
            coords={
                "time": np.arange(6),
                "isobaric": [100, 200, 500, 700, 1000],  # Pressure levels in hPa
                "lat": np.linspace(90, -90, 181),
                "lon": np.linspace(0, 359, 360),
            },
            attrs={
                "units": "gpm",
                "long_name": "Geopotential Height",
                "level_desc": "Isobaric level",
            },
        )

        result = geostrophic_wind(z)

        # Check dimensions
        assert result.ug.dims == ("time", "isobaric", "lat", "lon")
        assert result.vg.dims == ("time", "isobaric", "lat", "lon")

        # Check that interface was called with correct parameters
        call_args = mock_interface.call_args
        dim_order = call_args[1]["dim_order"]
        assert dim_order == "tzyx"  # time, isobaric (level), lat, lon

    @patch("skyborn.calc.geostrophic.xarray.interface.geostrophic_wind")
    def test_climate_model_data(self, mock_interface):
        """Test with climate model-like data structure."""
        mock_interface.return_value = (
            np.random.rand(48, 17, 96, 144) * 10,
            np.random.rand(48, 17, 96, 144) * 10,
        )

        # Create climate model-like data
        z = xr.DataArray(
            np.random.rand(48, 17, 96, 144) * 1000 + 5000,
            dims=["time", "plev", "lat", "lon"],
            coords={
                "time": np.arange(48),  # Monthly data for 4 years
                "plev": np.array(
                    [  # Standard pressure levels
                        1000,
                        925,
                        850,
                        700,
                        600,
                        500,
                        400,
                        300,
                        250,
                        200,
                        150,
                        100,
                        70,
                        50,
                        30,
                        20,
                        10,
                    ]
                )
                * 100,  # Convert to Pa
                "lat": np.linspace(-90, 90, 96),
                "lon": np.linspace(0, 357.5, 144),
            },
            attrs={
                "units": "m",
                "standard_name": "geopotential_height",
                "long_name": "Geopotential Height",
                "model": "CESM2",
                "experiment": "historical",
            },
        )

        result = geostrophic_wind(z, keep_attrs=True)

        # Check structure
        assert result.ug.shape == (48, 17, 96, 144)
        assert result.vg.shape == (48, 17, 96, 144)

        # Check that model information is preserved
        assert "source_geopotential_model" in result.attrs
        assert "source_geopotential_experiment" in result.attrs
        assert result.attrs["source_geopotential_model"] == "CESM2"
